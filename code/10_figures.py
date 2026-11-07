"""Step 10 - report figures and tables from the processed results (no downloads).
  fig_capture_and_wait.png   capture with 95 % CI, and how long each event waits for its first look
  fig_gap_vs_latitude.png    mean gap between looks at each station against latitude (the 1.4-day claim)
  fig_design_sweeps.png      satellites, usable swath and one-satellite options
  fig_fleet_availability.png the real look record: NISAR's 2026 hole, ALOS-2 absence, per-station gaps
  fig_altitude_scan.png      redrawn with only the main repeat-orbit traps labelled
  fig_capture_by_definition.png  the gain under four wetting-window rules
Tables: data/processed/table_*.csv (also printed as markdown)"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from armo_common import (PROC, FIG, ERA0, ERA1, load_events, stations_for, fleet_looks, constellation, merge, gap_table, score,
                         PHASES)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150,
                     "savefig.bbox": "tight"})
COL = {"Sentinel-1": "#7f8c8d", "NISAR": "#e67e22", "Sentinel-1 + NISAR": "#34495e", "ARMO D alone": "#2e86c1",
       "Fleet + ARMO D": "#1e8449"}
LAB = {"Sentinel-1": "Sentinel-1", "NISAR": "NISAR", "Sentinel-1 + NISAR": "Sentinel-1 + NISAR\n(today)",
       "ARMO D alone": "ARMO alone\n(2 satellites)", "Fleet + ARMO D": "Today + ARMO"}
TEAM = {"Sentinel-1": 30.2, "NISAR": 37.9, "Sentinel-1 + NISAR": 51.6, "ARMO D alone": 76.1, "Fleet + ARMO D": 81.5}


def md(df, floatfmt=1):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(f"{v:.{floatfmt}f}" if isinstance(v, float) else str(v) for v in r) + " |")
    return "\n".join(out)


if __name__ == "__main__":
    summ = pd.read_csv(os.path.join(PROC, "capture_summary.csv")).set_index("system")
    per = pd.read_parquet(os.path.join(PROC, "capture_events.parquet"))
    ev = load_events(); st = stations_for(ev); FL = fleet_looks()

    # ---------- capture and waiting time
    order = ["Sentinel-1", "NISAR", "Sentinel-1 + NISAR", "ARMO D alone", "Fleet + ARMO D"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0), gridspec_kw={"width_ratios": [1.05, 1]})
    x = np.arange(len(order))
    v = [100 * summ.loc[s, "capture"] for s in order]
    lo = [100 * (summ.loc[s, "capture"] - summ.loc[s, "ci_lo"]) for s in order]
    hi = [100 * (summ.loc[s, "ci_hi"] - summ.loc[s, "capture"]) for s in order]
    ax[0].bar(x, v, color=[COL[s] for s in order], yerr=[lo, hi], capsize=3, width=0.65)
    ax[0].scatter(x, [TEAM[s] for s in order], marker="_", s=400, color="#c0392b", zorder=5, label="team's report value")
    for i, val in enumerate(v):
        ax[0].text(i, val + hi[i] + 1.5, f"{val:.1f}%", ha="center", fontsize=8.5)
    ax[0].set_xticks(x); ax[0].set_xticklabels([LAB[s] for s in order], fontsize=8)
    ax[0].set_ylabel("wetting events with a look inside the window (%)"); ax[0].set_ylim(0, 100)
    ax[0].legend(frameon=False, fontsize=8, loc="upper left")
    ax[0].set_title(f"Capture, independent re-run ({len(ev)} events, {ev.wban.nunique()} stations, 17 Jun-31 Oct 2026)", fontsize=9)
    pooled = np.load(os.path.join(PROC, "pooled_waits.npz"))
    for s in ["Sentinel-1 + NISAR", "ARMO D alone", "Fleet + ARMO D"]:
        if s == "Sentinel-1 + NISAR":
            w_all = per[f"wait_{s}"].to_numpy(float)                       # one real record
        else:
            w_all = pooled[s.replace(" ", "_").replace("+", "plus")].ravel()  # 8 starting positions pooled
        w = np.sort(w_all[np.isfinite(w_all)])
        ax[1].plot(w / 24, np.arange(1, len(w) + 1) / len(w_all), color=COL[s], lw=1.8,
                   label=f"{LAB[s].replace(chr(10), ' ')}: median {np.nanmedian(w_all):.0f} h")
    ax[1].axvline(1, color="0.6", lw=0.8, ls=":"); ax[1].axvline(np.median(per.W_half_h) / 24, color="#c0392b", lw=0.8, ls="--")
    ax[1].text(np.median(per.W_half_h) / 24 + 0.05, 0.05, "median window\nlength", color="#c0392b", fontsize=7.5)
    ax[1].set_xlim(0, 7); ax[1].set_ylim(0, 1)
    ax[1].set_xlabel("days from the window opening to the first look"); ax[1].set_ylabel("share of events")
    ax[1].legend(frameon=False, fontsize=8, loc="lower right")
    ax[1].set_title("Waiting time: how long after the soil wets does the first look come?", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_capture_and_wait.png")); plt.close(fig)

    # ---------- mean gap against latitude (ARMO averaged over the 8 starting positions)
    rows = []
    for ph in PHASES:
        g = gap_table(constellation(st, n_sats=2, phase=ph, theta=(25, 45)), st.wban); g["phase"] = ph; rows.append(g)
    armo = pd.concat(rows).groupby("wban").agg(armo_gap=("mean_gap_d", "mean"), armo_max=("max_gap_d", "mean")).reset_index()
    comb_rows = []
    for ph in PHASES:
        g = gap_table(merge(constellation(st, n_sats=2, phase=ph, theta=(25, 45)), FL["Sentinel-1 + NISAR"]), st.wban)
        comb_rows.append(g)
    comb = pd.concat(comb_rows).groupby("wban").mean_gap_d.mean().rename("comb_gap").reset_index()
    gs = st.merge(armo, on="wban").merge(comb, on="wban")
    for name in ("Sentinel-1", "NISAR", "Sentinel-1 + NISAR"):
        g = gap_table(FL[name], st.wban)[["wban", "mean_gap_d", "max_gap_d"]].rename(
            columns={"mean_gap_d": f"{name}_gap", "max_gap_d": f"{name}_max"})
        gs = gs.merge(g, on="wban")
    gs.to_csv(os.path.join(PROC, "gap_by_station.csv"), index=False)
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    for name, colname, mk in (("Sentinel-1", "Sentinel-1_gap", "x"), ("NISAR", "NISAR_gap", "+"), ("Sentinel-1 + NISAR", "Sentinel-1 + NISAR_gap", "s"),
                              ("ARMO D alone", "armo_gap", "o"), ("Fleet + ARMO D", "comb_gap", "^")):
        ax.scatter(gs.lat, gs[colname], s=14, marker=mk, color=COL[name], label=LAB[name].replace("\n", " "), alpha=0.85)
    lat = gs.lat.to_numpy(); fit = np.polyfit(lat, gs.armo_gap, 1)
    xx = np.linspace(lat.min(), lat.max(), 50); ax.plot(xx, np.polyval(fit, xx), color=COL["ARMO D alone"], lw=1)
    ax.axhline(1.4, color="#c0392b", lw=0.9, ls="--"); ax.text(lat.max() + 0.3, 1.4, "1.4 d", color="#c0392b", va="center", fontsize=8)
    ax.set_yscale("log"); ax.set_yticks([0.5, 1, 1.4, 2, 3, 5, 10, 20]); ax.set_yticklabels(["0.5", "1", "1.4", "2", "3", "5", "10", "20"])
    ax.set_xlabel("station latitude (deg N)"); ax.set_ylabel("mean gap between looks, 17 Jun-31 Oct 2026 (days)")
    share = (gs.armo_gap < 1.4).mean()
    ax.set_title(f"ARMO alone: mean gap {gs.armo_gap.mean():.2f} d on average, below 1.4 d at {100 * share:.0f} % of stations "
                 f"(shorter further north)", fontsize=9)
    ax.legend(frameon=False, fontsize=7.5, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_gap_vs_latitude.png")); plt.close(fig)
    print(f"ARMO alone mean gap {gs.armo_gap.mean():.3f} d; < 1.4 d at {100 * share:.1f} % of stations; slope {fit[0]:.4f} d per deg;"
          f" south of 35N {gs[gs.lat < 35].armo_gap.mean():.2f} d, north of 42N {gs[gs.lat > 42].armo_gap.mean():.2f} d")

    # ---------- design sweeps
    sw = pd.read_csv(os.path.join(PROC, "design_sweeps.csv"))
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw={"width_ratios": [1, 1.15, 1.25]})
    ns = sw[sw.sweep == "n_sats"].sort_values("n_sats")
    ax[0].plot(ns.n_sats, 100 * ns.capture_alone, "o-", color=COL["ARMO D alone"], label="ARMO alone")
    ax[0].plot(ns.n_sats, 100 * ns.capture_with_fleet, "s-", color=COL["Fleet + ARMO D"], label="today + ARMO")
    ax[0].axhline(100 * sw.fleet_capture.iloc[0], color=COL["Sentinel-1 + NISAR"], ls=":", lw=1, label="today alone")
    for _, r in ns.iterrows():
        dy = 8 if r.n_sats == 1 else -27                     # keep the 1-satellite label clear of the "today" line
        ax[0].annotate(f"gap {r.mean_gap_d:.2f} d\nwait {r.wait_median_alone_h:.0f} h", (r.n_sats, 100 * r.capture_alone),
                       textcoords="offset points", xytext=(14 if r.n_sats == 1 else 0, dy), ha="left" if r.n_sats == 1 else "center",
                       fontsize=7, color=COL["ARMO D alone"])
    ax[0].set_xticks([1, 2, 3, 4]); ax[0].set_xlim(0.6, 4.4); ax[0].set_xlabel("satellites in one 620 km plane (362 km swath)")
    ax[0].set_ylabel("events captured (%)"); ax[0].set_ylim(30, 100)
    ax[0].legend(frameon=False, fontsize=7.5, loc="upper left"); ax[0].set_title("How many satellites?", fontsize=9)
    for n, mk, alpha in ((1, "s--", 0.55), (2, "o-", 1.0)):
        q = sw[(sw.sweep == "swath") & sw.label.str.startswith(f"{n} sat")].copy()
        q["sw"] = q.label.str.extract(r"swath (\d+)").astype(float); q = q.sort_values("sw")
        ax[1].plot(q.sw, 100 * q.capture_alone, mk, ms=4, alpha=alpha, color=COL["ARMO D alone"], label=f"{n} sat{'s' if n > 1 else ''}, alone")
        ax[1].plot(q.sw, 100 * q.capture_with_fleet, mk, ms=4, alpha=alpha, color=COL["Fleet + ARMO D"], label=f"{n} sat{'s' if n > 1 else ''}, today + ARMO")
    ax[1].axvspan(220, 300, color="#f5b041", alpha=0.15); ax[1].text(224, 92, "if only\n220-300 km\nis usable", fontsize=7.5, va="top")
    ax[1].axvline(362, color="0.5", lw=0.8, ls=":"); ax[1].text(366, 92, "362 km\naccess", fontsize=7.5, va="top")
    ax[1].set_xlabel("usable imaging swath (km)"); ax[1].set_ylim(30, 100); ax[1].legend(frameon=False, fontsize=7, loc="lower right")
    ax[1].set_title("What if the usable swath is narrower?", fontsize=9)
    names = {"2 satellite(s), one plane": "2 sats, 362 km (proposed)", "1 satellite(s), one plane": "1 sat, 362 km",
             "1 sat, alternate sides, 25-45 deg": "1 sat, alternate sides each orbit",
             "1 sat, 20-50 deg (574 km access)": "1 sat, 20-50 deg off-nadir",
             "1 sat at 800 km, 20-50 deg": "1 sat at 800 km, 20-50 deg",
             "1 sat, both sides every pass (bound)": "1 sat imaging both sides at once"}
    one = sw[sw.label.isin(names)].set_index("label").loc[list(names)].reset_index()
    yb = np.arange(len(one))[::-1]
    ax[2].barh(yb, 100 * one.capture_alone, color=[COL["ARMO D alone"] if "proposed" in names[l] else "#85c1e9" for l in one.label], height=0.6)
    for y_, (_, r) in zip(yb, one.iterrows()):
        ax[2].text(100 * r.capture_alone + 1, y_, f"{100 * r.capture_alone:.0f}%  ({r.swath_km:.0f} km, gap {r.mean_gap_d:.1f} d)",
                   va="center", fontsize=7.5)
    ax[2].set_yticks(yb); ax[2].set_yticklabels([names[l] for l in one.label], fontsize=8); ax[2].set_xlim(0, 115)
    ax[2].set_xlabel("events captured by ARMO alone (%)"); ax[2].set_title("Can one satellite do the job?", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_design_sweeps.png")); plt.close(fig)

    # ---------- the real fleet record
    s1 = pd.read_csv(os.path.join(PROC, "looks_s1.csv"), parse_dates=["time"]); ni = pd.read_csv(os.path.join(PROC, "looks_nisar.csv"), parse_dates=["time"])
    al = pd.read_csv(os.path.join(PROC, "looks_alos2.csv"), parse_dates=["time"])
    allst = pd.read_csv(os.path.join(PROC, "station_seasons.csv")); n_st = allst[allst.included].wban.nunique()
    months = pd.period_range("2025-04", "2026-10", freq="M")
    def per_station_month(d):
        c = d.groupby(d.time.dt.tz_localize(None).dt.to_period("M")).size().reindex(months, fill_value=0)
        return c / n_st
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 3.8), gridspec_kw={"width_ratios": [1.5, 1]})
    for d, name in ((s1, "Sentinel-1"), (ni, "NISAR")):
        ax[0].plot(months.to_timestamp(), per_station_month(d), "o-", ms=3, color=COL[name], label=name)
    ax[0].plot(months.to_timestamp(), per_station_month(al), "o-", ms=3, color="#8e44ad", label="ALOS-2 (NASA archive)")
    ax[0].axvspan(pd.Timestamp("2026-02-01"), pd.Timestamp("2026-06-17"), color="#e67e22", alpha=0.12)
    ax[0].text(pd.Timestamp("2026-02-05"), 8.6, "NISAR products\nmissing at most\nstations", fontsize=7.5, color="#a04000")
    ax[0].set_ylabel("looks per station per month"); ax[0].set_ylim(0, 10); ax[0].legend(frameon=False, fontsize=8, loc="upper left")
    ax[0].set_title(f"What field-scale radar actually looked at {n_st} USCRN stations", fontsize=9)
    era = lambda d: d[(d.time >= ERA0) & (d.time < ERA1)]
    for d, name in ((s1, "Sentinel-1"), (ni, "NISAR")):
        g = gap_table({w: np.sort(((x.time - pd.Timestamp("1970-01-01", tz="UTC")) / pd.Timedelta(hours=1)).to_numpy()) for w, x in era(d).groupby("wban")}, st.wban)
        ax[1].hist(g.mean_gap_d.clip(upper=30), bins=np.arange(0, 31, 1), color=COL[name], alpha=0.6, label=f"{name}: median {g.mean_gap_d.median():.1f} d")
    ax[1].hist(gs.armo_gap, bins=np.arange(0, 31, 1), color=COL["ARMO D alone"], alpha=0.8, label=f"ARMO alone: median {gs.armo_gap.median():.1f} d")
    ax[1].set_xlabel("mean gap between looks at a station, 17 Jun-31 Oct 2026 (days; 30 = 30 or more)"); ax[1].set_ylabel("stations")
    ax[1].legend(frameon=False, fontsize=7.5); ax[1].set_title("Spread between stations", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_fleet_availability.png")); plt.close(fig)

    # ---------- altitude scan, redrawn
    alt = pd.read_csv(os.path.join(PROC, "altitude_scan.csv"))
    fig, ax = plt.subplots(2, 1, figsize=(9, 5.4), sharex=True, gridspec_kw={"height_ratios": [1.4, 1]})
    ax[0].plot(alt.h_km, 100 * alt.capture_alone, "o-", ms=3, color=COL["ARMO D alone"], label="ARMO alone (2 satellites)")
    ax[0].fill_between(alt.h_km, 100 * alt.capture_alone_min_phase, 100 * alt.capture_alone_max_phase, color=COL["ARMO D alone"], alpha=0.12)
    ax[0].plot(alt.h_km, 100 * alt.capture_with_fleet, "s-", ms=3, color=COL["Fleet + ARMO D"], label="today + ARMO")
    ax[0].axvline(620, color="0.4", lw=0.8, ls=":"); ax[0].text(622, 46, "620 km\n(proposed)", fontsize=8)
    ax[0].set_ylabel("events captured (%)"); ax[0].legend(frameon=False, fontsize=8, loc="lower right")
    ax[0].set_title("Same 362 km swath at every altitude: capture depends on the orbit's repeat pattern", fontsize=9)
    ax[1].plot(alt.h_km, alt.worst_gap_d, "o-", ms=3, color="#c0392b")
    for _, r in alt[alt.worst_gap_d > 20].iterrows():
        rep = {"15/1": "repeats daily", "29/2": "repeats every 2 days", "59/4": "repeats every 4 days"}.get(r.near_fraction, r.near_fraction)
        ax[1].annotate(f"{r.h_km} km: {rep}", (r.h_km, r.worst_gap_d), textcoords="offset points", xytext=(6, -2), fontsize=7.5)
    ax[1].set_yscale("log"); ax[1].set_ylabel("worst gap at\nany station (days)"); ax[1].set_xlabel("orbit altitude (km)")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_altitude_scan.png")); plt.close(fig)

    # ---------- capture by wetting-window definition
    bd = pd.read_csv(os.path.join(PROC, "capture_by_definition.csv"))
    pv = bd.pivot_table(index="system", columns="definition", values="capture")
    defs = [("half", "50 % rule\n(report)"), ("half_p3", "50 %, 3 h\nbelow to end"), ("efold", "1/e of the\nrise"), ("abs", "anomaly\n>= 0.03")]
    fig, ax = plt.subplots(figsize=(7.5, 3.6)); xb = np.arange(len(defs)); wdt = 0.26
    for k, s in enumerate(["Sentinel-1 + NISAR", "ARMO D alone", "Fleet + ARMO D"]):
        vals = [100 * pv.loc[s, d] for d, _ in defs]
        ax.bar(xb + (k - 1) * wdt, vals, width=wdt, color=COL[s], label=LAB[s].replace("\n", " "))
        for i, val in enumerate(vals):
            ax.text(xb[i] + (k - 1) * wdt, val + 1, f"{val:.0f}", ha="center", fontsize=7)
    for i, (d, _) in enumerate(defs):
        ax.text(xb[i], 4, f"+{100 * (pv.loc['Fleet + ARMO D', d] - pv.loc['Sentinel-1 + NISAR', d]):.0f} pts", ha="center", fontsize=8,
                color="white", fontweight="bold")
    ax.set_xticks(xb); ax.set_xticklabels([l for _, l in defs], fontsize=8); ax.set_ylabel("events captured (%)"); ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=7.5, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    ax.set_title("ARMO's gain holds under every window rule tested (white: gain of today + ARMO over today)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_capture_by_definition.png")); plt.close(fig)

    # ---------- tables
    t1 = pd.DataFrame([dict(System=LAB[s].replace("\n", " "), **{"Capture %": 100 * summ.loc[s, "capture"],
                       "95% CI": f"{100 * summ.loc[s, 'ci_lo']:.1f}-{100 * summ.loc[s, 'ci_hi']:.1f}", "Team %": TEAM[s],
                       "Median wait h": summ.loc[s, "wait_median_h"], "90th pct wait h": summ.loc[s, "wait_p90_h"],
                       "Looked at within 24 h %": 100 * summ.loc[s, "wait_le_24h"], "Mean gap d": summ.loc[s, "mean_gap_d"],
                       "Stations < 1.4 d %": 100 * summ.loc[s, "share_stations_lt_1p4"], "Worst gap d": summ.loc[s, "worst_gap_d"]})
                       for s in order])
    t1.to_csv(os.path.join(PROC, "table_capture_wait.csv"), index=False); print(md(t1))
    print("figures written")
