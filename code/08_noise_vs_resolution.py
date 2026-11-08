"""Step 8 - measured L-band change noise against product resolution and against days between passes (real NISAR data).

Question: how fine can ARMO's soil-moisture product be, and is the team's assumption right that today's 1.45 dB change
noise is speckle that more independent looks would remove (the report asks for about 15x more looks to reach 0.38 dB)?

Data: NISAR GCOV gamma0 (10 m, HH and HV) windows of 5.12 x 5.12 km around 16 USCRN stations, 17 Jun-31 Oct 2026
(07_nisar_windows.py). Pairs = two passes on the same track and pass direction (same geometry), 12, 24 or 36 days apart.
Ground truth: USCRN 5 cm soil moisture at the two overpass hours; pairs with >= 0.2 mm rain in the 3 h before either pass
are dropped (wet canopy, ponding), as in the team's method. "Stable" pairs: |dSM| < 0.01 m3/m3.

For block sizes of 10-500 m (pixels averaged in linear power):
  spatial change noise = robust spread (1.4826 x MAD) of the block-by-block dB change across the window, after removing
                         the window median (any uniform calibration or moisture shift); this is the noise a change detector
                         sees at that product resolution, including speckle, vegetation and roughness change.
  speckle-only theory  = 4.343 * sqrt(2 * trigamma(N)) dB for N independent looks per block, N = blocks pixels x ENL,
                         ENL per 10 m pixel from the file's own resolution (5.35 m along track x 3.75 m / sin(incidence)).
Sensitivity: dB change of the station's 100 m and 200 m block against dSM over all pairs (least squares, station bootstrap).
Outputs: data/processed/noise_pairs.csv, noise_vs_scale.csv, noise_sensitivity.csv; figures/fig_noise_vs_resolution.png"""
import os, glob, json
import numpy as np, pandas as pd
from scipy.special import polygamma
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROC = os.path.join(ROOT, "data", "processed")
RAW = os.path.join(ROOT, "data", "raw", "nisar_windows"); FIG = os.path.join(ROOT, "figures")
BLOCKS = [1, 2, 3, 5, 10, 15, 20, 30, 50]          # pixels of 10 m -> 10 ... 500 m
AZ_RES, RG_RES, INC_DEG = 5.35, 3.75, 40.0         # from the GCOV metadata (sceneCenterAlongTrackResolution, rangeResolution)
ENL10 = 100.0 / (AZ_RES * RG_RES / np.sin(np.radians(INC_DEG)))
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})


def speckle_db(n):
    return 10 / np.log(10) * np.sqrt(2 * polygamma(1, n))


def block_means(img, valid, b):
    n, m = img.shape[0] // b, img.shape[1] // b
    v2 = valid[:n * b, :m * b]
    x = np.where(v2, img[:n * b, :m * b], 0).reshape(n, b, m, b); v = v2.reshape(n, b, m, b)
    cnt = v.sum((1, 3)); s = x.sum((1, 3))
    with np.errstate(invalid="ignore", divide="ignore"):
        return s / cnt, cnt / (b * b)


def load_index():
    rows = []
    for f in glob.glob(os.path.join(RAW, "*", "*.npz")):
        z = np.load(f)
        rows.append(dict(wban=int(os.path.basename(os.path.dirname(f))), path=f, track=int(z["track"]), frame=int(z["frame"]),
                         direction=str(z["direction"]), t=pd.Timestamp(str(z["t"]), tz="UTC"), epsg=int(z["epsg"]),
                         x0=float(z["x0"]), y0=float(z["y0"]), dx=float(z["dx"]), dy=float(z["dy"]), sr=int(z["sr"]), sc=int(z["sc"])))
    w = pd.DataFrame(rows).sort_values(["wban", "t"])
    # the same pass can sit in two neighbouring frames: keep one window per pass
    w = w[(w.wban != w.wban.shift()) | ((w.t - w.t.shift()) > pd.Timedelta(hours=1))]
    return w.reset_index(drop=True)


def coregister(A, B, ra, rb):
    """Crop two tiles to their common ground area. Returns (a_slices, b_slices, station row, station col) or None."""
    if ra.epsg != rb.epsg or ra.dx != rb.dx or ra.dy != rb.dy:
        return None
    oj = int(round((ra.x0 - rb.x0) / ra.dx)); oi = int(round((ra.y0 - rb.y0) / ra.dy))
    n = A["hh"].shape[0]
    i0, i1 = max(0, -oi), min(n, n - oi); j0, j1 = max(0, -oj), min(n, n - oj)
    if i1 - i0 < 256 or j1 - j0 < 256:
        return None
    return (slice(i0, i1), slice(j0, j1)), (slice(i0 + oi, i1 + oi), slice(j0 + oj, j1 + oj)), ra.sr - i0, ra.sc - j0


if __name__ == "__main__":
    w = load_index()
    hourly = pd.read_parquet(os.path.join(PROC, "uscrn_hourly.parquet"), columns=["wban", "time", "sm5", "rain"])
    hourly = hourly[hourly.wban.isin(w.wban.unique())].set_index(["wban", "time"]).sort_index()
    def sm_at(wban, t):
        h = t.floor("h")
        try:
            g = hourly.loc[wban]
            sm = g.sm5.reindex([h, h + pd.Timedelta(hours=1)]).mean()
            rain3 = g.rain.reindex(pd.date_range(h - pd.Timedelta(hours=3), h, freq="h")).sum()
            return sm, rain3
        except KeyError:
            return np.nan, np.nan
    w[["sm", "rain3"]] = [sm_at(r.wban, r.t) for r in w.itertuples()]
    rows = []
    for (wban, track, direc), g in w.groupby(["wban", "track", "direction"]):
        g = g.sort_values("t").reset_index(drop=True)
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                dt = (g.t[j] - g.t[i]).total_seconds() / 86400
                lag = int(round(dt / 12.0))
                if lag not in (1, 2, 3) or abs(dt - 12 * lag) > 1.0:
                    continue
                A, B = np.load(g.path[i]), np.load(g.path[j])
                cr = coregister(A, B, g.iloc[i], g.iloc[j])
                if cr is None:
                    continue
                sa, sb, srow, scol = cr
                rec = dict(wban=wban, track=track, direction=direc, t1=g.t[i], t2=g.t[j], lag_days=12 * lag,
                           sm1=g.sm[i], sm2=g.sm[j], dsm=g.sm[j] - g.sm[i], rain3=max(g.rain3[i], g.rain3[j]))
                for pol in ("hh", "hv"):
                    a, b_ = A[pol][sa].astype(float), B[pol][sb].astype(float)
                    valid = np.isfinite(a) & np.isfinite(b_) & (a > 0) & (b_ > 0)
                    rec[f"{pol}_valid"] = valid.mean()
                    for bs in BLOCKS:
                        ma, fa = block_means(a, valid, bs); mb, fb = block_means(b_, valid, bs)
                        ok = (fa >= 0.9) & (fb >= 0.9)
                        d = 10 * np.log10(mb[ok] / ma[ok])
                        if d.size >= 25:
                            rec[f"{pol}_noise_{bs * 10}m"] = 1.4826 * np.median(np.abs(d - np.median(d)))
                            rec[f"{pol}_shift_{bs * 10}m"] = np.median(d)
                        kr, kc = srow // bs, scol // bs
                        if kr < ok.shape[0] and kc < ok.shape[1] and ok[kr, kc]:
                            rec[f"{pol}_station_{bs * 10}m"] = 10 * np.log10(mb[kr, kc] / ma[kr, kc])
                    rec[f"{pol}_level_db"] = 10 * np.log10(np.nanmedian(a[valid])) if valid.any() else np.nan
                rows.append(rec)
    pr = pd.DataFrame(rows)
    pr.to_csv(os.path.join(PROC, "noise_pairs.csv"), index=False)
    dry = pr.rain3.fillna(0) < 0.2
    stable = dry & (pr.dsm.abs() < 0.01)
    print(f"pairs {len(pr)} at {pr.wban.nunique()} stations | dry {dry.sum()} | stable {stable.sum()} | by lag",
          pr[stable].lag_days.value_counts().to_dict())

    # noise against product resolution (stable 12-day pairs) and against lag
    out = []
    for pol in ("hh", "hv"):
        for bs in BLOCKS:
            col = f"{pol}_noise_{bs * 10}m"
            for lag in (12, 24, 36):
                v = pr.loc[stable & (pr.lag_days == lag), col].dropna()
                if len(v):
                    out.append(dict(pol=pol.upper(), scale_m=bs * 10, lag_days=lag, n_pairs=len(v), median_db=v.median(),
                                    q25_db=v.quantile(0.25), q75_db=v.quantile(0.75),
                                    speckle_theory_db=speckle_db(bs * bs * ENL10)))
            # station-block temporal spread on stable 12-day pairs
            s = pr.loc[stable & (pr.lag_days == 12), f"{pol}_station_{bs * 10}m"].dropna()
            if len(s) >= 5:
                out.append(dict(pol=pol.upper(), scale_m=bs * 10, lag_days=12, n_pairs=len(s), median_db=np.nan,
                                station_block_robust_sd_db=1.4826 * np.median(np.abs(s - s.median())),
                                station_block_sd_db=s.std(), speckle_theory_db=speckle_db(bs * bs * ENL10), kind="station block"))
    ns = pd.DataFrame(out); ns.to_csv(os.path.join(PROC, "noise_vs_scale.csv"), index=False)
    print(ns[ns.kind.isna() if "kind" in ns else slice(None)].query("lag_days == 12").round(3).to_string(index=False))

    # sensitivity: station-block change against dSM (dry pairs, 12-day)
    sens = []
    rng = np.random.default_rng(1)
    for pol in ("hh", "hv"):
        for bs in (10, 20):
            col = f"{pol}_station_{bs * 10}m"
            d = pr[dry & (pr.lag_days == 12)].dropna(subset=[col, "dsm"])
            if len(d) < 10:
                continue
            slope = np.polyfit(d.dsm, d[col], 1)[0]
            u = d.wban.unique(); bss = []
            for _ in range(1000):
                pick = rng.choice(u, len(u)); dd = pd.concat([d[d.wban == x] for x in pick])
                if dd.dsm.std() > 0: bss.append(np.polyfit(dd.dsm, dd[col], 1)[0])
            r = np.corrcoef(d.dsm, d[col])[0, 1]
            sens.append(dict(pol=pol.upper(), scale_m=bs * 10, n_pairs=len(d), stations=d.wban.nunique(),
                             slope_db_per_0p1=slope / 10, ci_lo=np.percentile(bss, 2.5) / 10, ci_hi=np.percentile(bss, 97.5) / 10, r=r))
    se = pd.DataFrame(sens); se.to_csv(os.path.join(PROC, "noise_sensitivity.csv"), index=False)
    print(se.round(3).to_string(index=False))

    # effective looks per 10 m pixel, fitted where speckle dominates (10-30 m), HH and HV together
    from scipy.optimize import minimize_scalar
    fine = ns[ns.kind.isna() & (ns.lag_days == 12) & ns.scale_m.isin([10, 20, 30])]
    obj = lambda e: np.sum((np.log(fine.median_db) - np.log(speckle_db((fine.scale_m / 10) ** 2 * e))) ** 2)
    ENL_FIT = minimize_scalar(obj, bounds=(0.5, 50), method="bounded").x
    print(f"fitted looks per 10 m pixel: {ENL_FIT:.1f} (resolution-based estimate {ENL10:.1f})")
    # station-level noise: residual SD of station-block change after removing slope x dSM (dry 12-day pairs)
    resid = []
    for pol in ("hh", "hv"):
        for bs in BLOCKS:
            col = f"{pol}_station_{bs * 10}m"; d = pr[dry & (pr.lag_days == 12)].dropna(subset=[col, "dsm"])
            if len(d) >= 20:
                b = np.polyfit(d.dsm, d[col], 1); r_ = d[col] - np.polyval(b, d.dsm)
                resid.append(dict(pol=pol.upper(), scale_m=bs * 10, n_pairs=len(d), resid_sd_db=r_.std(),
                                  resid_robust_db=1.4826 * np.median(np.abs(r_ - r_.median()))))
    rs = pd.DataFrame(resid); rs.to_csv(os.path.join(PROC, "noise_station_residual.csv"), index=False)
    print(rs.round(3).to_string(index=False))
    json.dump(dict(enl_fit=ENL_FIT, enl_resolution=ENL10), open(os.path.join(PROC, "noise_enl.json"), "w"))
    # speckle-only references: NISAR with the fitted looks; ARMO burst modes from sar_tradeoff.py (feas_modes.csv) section 2
    sc = ns.scale_m.astype(float)
    ns["speckle_nisar_fit_db"] = speckle_db((sc / 10) ** 2 * ENL_FIT)
    ns["speckle_armo_8b_40mhz_db"] = np.where(sc >= 48, speckle_db(np.maximum((sc / 48.0) * (sc / 7.1), 1)), np.nan)
    ns["speckle_armo_11b_20mhz_db"] = np.where(sc >= 66, speckle_db(np.maximum((sc / 66.0) * (sc / 14.3), 1)), np.nan)
    ns.to_csv(os.path.join(PROC, "noise_vs_scale.csv"), index=False)

    # figure
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3), gridspec_kw={"width_ratios": [1.45, 1]})
    xs = np.array(BLOCKS) * 10
    for pol, col in (("HH", "#1f4e79"), ("HV", "#8e44ad")):
        q = ns[(ns.pol == pol) & (ns.lag_days == 12) & ns.median_db.notna()].sort_values("scale_m")
        ax[0].plot(q.scale_m, q.median_db, "o-", color=col, label=f"NISAR {pol}, measured (stable 12-day pairs)")
        ax[0].fill_between(q.scale_m, q.q25_db, q.q75_db, color=col, alpha=0.15)
    xf = np.geomspace(10, 500, 60)
    ax[0].plot(xf, speckle_db((xf / 10) ** 2 * ENL_FIT), "k--", lw=1, label=f"speckle only, NISAR 40 MHz ({ENL_FIT:.1f} looks per 10 m pixel, fitted)")
    for (az, rg, lab, ls_) in ((48.0, 7.1, "speckle only, ARMO burst mode 8 beams, 40 MHz", "-."), (66.0, 14.3, "speckle only, ARMO burst mode 11 beams, 20 MHz", ":")):
        n_ = np.maximum((xf / az) * (xf / rg), 1.0); ok_ = xf >= az
        ax[0].plot(xf[ok_], speckle_db(n_[ok_]), color="#16a085", ls=ls_, lw=1.3, label=lab)
    for pol, mk in (("HH", "D"), ("HV", "v")):
        q = rs[rs.pol == pol]
        ax[0].plot(q.scale_m, q.resid_sd_db, mk, color="#c0392b" if pol == "HH" else "#e59866", ms=4, ls="none",
                   label=f"{pol} at the station: change not explained by soil moisture (all noise)")
    ax[0].axhline(0.38, color="#c0392b", lw=1, ls=":"); ax[0].text(11, 0.40, "team's 0.38 dB requirement", color="#c0392b", fontsize=8)
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("product resolution (m, pixels averaged)")
    ax[0].set_xticks([10, 20, 50, 100, 200, 500]); ax[0].set_xticklabels(["10", "20", "50", "100", "200", "500"])
    ax[0].set_yticks([0.1, 0.2, 0.5, 1, 2, 5]); ax[0].set_yticklabels(["0.1", "0.2", "0.5", "1", "2", "5"])
    ax[0].set_ylabel("backscatter change noise (dB)"); ax[0].legend(frameon=False, fontsize=6.8, loc="lower left")
    ax[0].set_ylim(0.05, 10)
    ax[0].set_title("NISAR change noise against product resolution (stable 12-day pairs)", fontsize=9)
    for pol, col in (("HH", "#1f4e79"), ("HV", "#8e44ad")):
        for scale, ls in ((100, "-"), (200, "--")):
            q = ns[(ns.pol == pol) & (ns.scale_m == scale) & ns.median_db.notna()].sort_values("lag_days")
            if len(q):
                ax[1].plot(q.lag_days, q.median_db, "o" + ls, color=col, label=f"{pol}, {scale} m")
    ax[1].axhline(0.38, color="#c0392b", lw=1, ls=":")
    ax[1].set_xlabel("days between the two passes"); ax[1].set_ylabel("change noise on stable pairs (dB)")
    ax[1].set_xticks([12, 24, 36]); ax[1].legend(frameon=False, fontsize=8)
    ax[1].set_title("Does the noise grow with the time between passes?", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_noise_vs_resolution.png")); plt.close(fig)
    print("figure written")
