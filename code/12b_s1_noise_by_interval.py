"""Step 12b - Sentinel-1 change noise against the number of days between the two images (6, 12, 18, 24 days).
Same measures as step 8: window-wide spatial noise on stable pairs at 100 m and 200 m, and the station-level change not
explained by the soil-moisture change. Each pair = same relative orbit (same geometry).
Outputs: data/processed/s1_noise_pairs.csv, s1_noise_by_lag.csv, figures/fig_s1_noise_vs_interval.png"""
import os, glob, importlib.util
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROC = os.path.join(ROOT, "data", "processed")
RAW = os.path.join(ROOT, "data", "raw", "s1_windows"); FIG = os.path.join(ROOT, "figures")
spec = importlib.util.spec_from_file_location("n8", os.path.join(ROOT, "code", "08_noise_vs_resolution.py"))
n8 = importlib.util.module_from_spec(spec); spec.loader.exec_module(n8)
SCALES = [5, 10, 20]                    # 10 m pixels -> 50, 100, 200 m
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})
rng = np.random.default_rng(7)

if __name__ == "__main__":
    rows = []
    for f in glob.glob(os.path.join(RAW, "*", "*.npz")):
        z = np.load(f)
        rows.append(dict(wban=int(os.path.basename(os.path.dirname(f))), path=f, rel_orbit=int(z["rel_orbit"]),
                         t=pd.Timestamp(str(z["t"]), tz="UTC"), crs=str(z["crs"]), x0=float(z["x0"]), y0=float(z["y0"]),
                         dx=float(z["dx"]), dy=float(z["dy"])))
    w = pd.DataFrame(rows).sort_values(["wban", "t"])
    w = w[(w.wban != w.wban.shift()) | ((w.t - w.t.shift()) > pd.Timedelta(hours=1))].reset_index(drop=True)
    hourly = pd.read_parquet(os.path.join(PROC, "uscrn_hourly.parquet"), columns=["wban", "time", "sm5", "rain"])
    hourly = hourly[hourly.wban.isin(w.wban.unique())].set_index(["wban", "time"]).sort_index()
    def sm_at(wban, t):
        h = t.floor("h"); g = hourly.loc[wban]
        return g.sm5.reindex([h, h + pd.Timedelta(hours=1)]).mean(), g.rain.reindex(pd.date_range(h - pd.Timedelta(hours=3), h, freq="h")).sum()
    w[["sm", "rain3"]] = [sm_at(r.wban, r.t) for r in w.itertuples()]
    out = []
    for (wban, ro), g in w.groupby(["wban", "rel_orbit"]):
        g = g.sort_values("t").reset_index(drop=True)
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                dt = (g.t[j] - g.t[i]).total_seconds() / 86400; lag = int(round(dt / 6.0)) * 6
                if lag not in (6, 12, 18, 24) or abs(dt - lag) > 0.5:
                    continue
                ra, rb = g.iloc[i], g.iloc[j]
                if ra.crs != rb.crs:
                    continue
                A, B = np.load(ra.path), np.load(rb.path)
                oj = int(round((ra.x0 - rb.x0) / ra.dx)); oi = int(round((ra.y0 - rb.y0) / ra.dy)); n = A["vv"].shape[0]
                i0, i1 = max(0, -oi), min(n, n - oi); j0, j1 = max(0, -oj), min(n, n - oj)
                if i1 - i0 < 256 or j1 - j0 < 256:
                    continue
                sa = (slice(i0, i1), slice(j0, j1)); sb = (slice(i0 + oi, i1 + oi), slice(j0 + oj, j1 + oj))
                srow, scol = n8.HALF if hasattr(n8, "HALF") else 256, 256
                srow, scol = 256 - i0, 256 - j0
                rec = dict(wban=wban, rel_orbit=ro, t1=ra.t, t2=rb.t, lag_days=lag, dsm=rb.sm - ra.sm, rain3=max(ra.rain3, rb.rain3))
                for pol in ("vv", "vh"):
                    a, b_ = A[pol][sa].astype(float), B[pol][sb].astype(float)
                    valid = np.isfinite(a) & np.isfinite(b_) & (a > 0) & (b_ > 0)
                    for bs in SCALES:
                        ma, fa = n8.block_means(a, valid, bs); mb, fb = n8.block_means(b_, valid, bs)
                        ok = (fa >= 0.9) & (fb >= 0.9); d = 10 * np.log10(mb[ok] / ma[ok])
                        if d.size >= 25:
                            rec[f"{pol}_noise_{bs * 10}m"] = 1.4826 * np.median(np.abs(d - np.median(d)))
                        kr, kc = srow // bs, scol // bs
                        if kr < ok.shape[0] and kc < ok.shape[1] and ok[kr, kc]:
                            rec[f"{pol}_station_{bs * 10}m"] = 10 * np.log10(mb[kr, kc] / ma[kr, kc])
                out.append(rec)
    pr = pd.DataFrame(out); pr.to_csv(os.path.join(PROC, "s1_noise_pairs.csv"), index=False)
    dry = pr.rain3.fillna(0) < 0.2; stable = dry & (pr.dsm.abs() < 0.01)
    print(f"pairs {len(pr)} at {pr.wban.nunique()} stations | stable by lag", pr[stable].lag_days.value_counts().sort_index().to_dict())
    res = []
    for pol in ("vv", "vh"):
        for bs in SCALES:
            sc = bs * 10
            for lag in (6, 12, 18, 24):
                m = stable & (pr.lag_days == lag)
                v = pr.loc[m, f"{pol}_noise_{sc}m"].dropna()
                # paired by station: only stations that have both 6- and 12-day stable pairs, for a fair comparison
                d = pr[dry & (pr.lag_days == lag)].dropna(subset=[f"{pol}_station_{sc}m", "dsm"])
                resid = np.nan; resid_rb = np.nan
                if len(d) >= 20:
                    b = np.polyfit(d.dsm, d[f"{pol}_station_{sc}m"], 1); r_ = d[f"{pol}_station_{sc}m"] - np.polyval(b, d.dsm)
                    resid = r_.std(); resid_rb = 1.4826 * np.median(np.abs(r_ - r_.median()))
                res.append(dict(pol=pol.upper(), scale_m=sc, lag_days=lag, n_stable=len(v), stations=pr.loc[m, "wban"].nunique(),
                                spatial_median_db=v.median() if len(v) else np.nan, n_dry=len(d), station_resid_sd_db=resid,
                                station_resid_robust_db=resid_rb))
    rs = pd.DataFrame(res)
    # station-matched 6 vs 12 day comparison of the spatial noise (stations with stable pairs at both lags)
    both = []
    for pol in ("vv", "vh"):
        for sc in (100, 200):
            c = f"{pol}_noise_{sc}m"; s6 = pr[stable & (pr.lag_days == 6)].groupby("wban")[c].median(); s12 = pr[stable & (pr.lag_days == 12)].groupby("wban")[c].median()
            k = s6.index.intersection(s12.index); r = (s12[k] / s6[k]).dropna()
            if len(r) >= 3:
                bs_ = [np.median(rng.choice(r.to_numpy(), len(r))) for _ in range(2000)]
                both.append(dict(pol=pol.upper(), scale_m=sc, stations=len(r), median_ratio_12_over_6=np.median(r),
                                 ci_lo=np.percentile(bs_, 2.5), ci_hi=np.percentile(bs_, 97.5)))
    bt = pd.DataFrame(both)
    rs.to_csv(os.path.join(PROC, "s1_noise_by_lag.csv"), index=False); bt.to_csv(os.path.join(PROC, "s1_noise_6v12_matched.csv"), index=False)
    print(rs.round(3).to_string(index=False)); print(bt.round(3).to_string(index=False))
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for pol, col in (("VV", "#1f4e79"), ("VH", "#8e44ad")):
        for sc, ls in ((100, "-"), (200, "--")):
            q = rs[(rs.pol == pol) & (rs.scale_m == sc)].sort_values("lag_days")
            ax[0].plot(q.lag_days, q.spatial_median_db, "o" + ls, color=col, label=f"{pol}, {sc} m")
            ax[1].plot(q.lag_days, q.station_resid_sd_db, "o" + ls, color=col, label=f"{pol}, {sc} m")
    ax[0].set_title("Window-wide change noise on stable pairs", fontsize=9); ax[0].set_ylabel("dB")
    ax[1].set_title("At the station: change not explained by soil moisture", fontsize=9)
    for a in ax:
        a.set_xticks([6, 12, 18, 24]); a.set_xlabel("days between the two Sentinel-1 images (same track)"); a.legend(frameon=False, fontsize=8)
        a.set_ylim(0, None)
    fig.suptitle("Sentinel-1 (C-band), USCRN stations, Apr-Oct 2026: does a shorter interval lower the noise?", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_s1_noise_vs_interval.png")); plt.close(fig)
    print("figure written")
