"""Step 5 - does ARMO catch more wetting events, and sooner? Independent re-check of the report's Tables 3-4.

Systems: Sentinel-1, NISAR and the fleet (both) from real look times; ARMO designs A-E simulated with the orbit model in
04_orbit_sim.py and averaged over 8 starting positions, alone and merged with the real fleet.
  A  1 satellite, right-looking 25-45 deg off-nadir     B  1 satellite, left/right on alternate orbits, 25-45 deg
  C  1 satellite, left/right on alternate orbits, 20-50 deg (as in the team's code/06)
  D  2 satellites 180 deg apart, right-looking 25-45 deg (proposed)
  E  2 satellites 180 deg apart, alternate sides
All at 620 km, dawn-dusk sun-synchronous (ascending node 18:00).
Per system: capture (with station-bootstrap 95 % CI), waiting time from window opening to first look (mean, median, 90th
percentile, share within 24 h), mean gap between looks per station, worst gap. Repeated for the four window definitions.
Outputs: data/processed/capture_summary.csv, capture_by_definition.csv, capture_events.parquet, gaps_station.csv"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from armo_common import (PROC, DEFS, load_events, stations_for, fleet_looks, constellation, merge, score, gap_table,
                         gap_summary, boot_ci, phase_average, wait_stats, PHASES)

DESIGNS = {"A": dict(n_sats=1, theta=(25, 45), side="right"), "B": dict(n_sats=1, theta=(25, 45), side="alt"),
           "C": dict(n_sats=1, theta=(20, 50), side="alt"), "D": dict(n_sats=2, theta=(25, 45), side="right"),
           "E": dict(n_sats=2, theta=(25, 45), side="alt")}

if __name__ == "__main__":
    ev = load_events(); st = stations_for(ev); FLEET = fleet_looks()
    print(f"cohort: {len(ev)} events at {ev.wban.nunique()} stations")
    rows, by_def, gaps_all = [], [], []
    per_event = ev[["wban", "t_start", "amp", "rain_mm", "W_half_h"]].copy()
    pooled_waits = {}                     # phases x events waits for the simulated systems (for distributions)

    def record(name, cap, wait, looks_for_gaps, d="half"):
        lo, hi = boot_ci(ev, cap)
        r = dict(system=name, definition=d, events=len(ev), capture=cap.mean(), ci_lo=lo, ci_hi=hi, **wait_stats(wait))
        if looks_for_gaps is not None:
            gt = gap_table(looks_for_gaps, st.wban); gt["system"] = name
            if d == "half":
                gaps_all.append(gt.drop(columns="gaps"))
            r.update(gap_summary(gt))
        return r

    for d in DEFS:
        for name, L in FLEET.items():
            c, w = score(ev, L, d)
            (rows if d == "half" else by_def).append(record(name, c.astype(float), w, L if d == "half" else None, d))
            if d == "half":
                per_event[f"cap_{name}"], per_event[f"wait_{name}"] = c, w
        for dn, kw in DESIGNS.items():
            if d != "half" and dn != "D":
                continue
            kw2 = {k: v for k, v in kw.items() if k != "n_sats"}
            mk = lambda ph: constellation(st, n_sats=kw["n_sats"], phase=ph, **kw2)
            ca, wa = phase_average(ev, mk, d, per_phase=True)
            cf, wf = phase_average(ev, mk, d, base=FLEET["Sentinel-1 + NISAR"], per_phase=True)
            out = rows if d == "half" else by_def
            ra = record(f"ARMO {dn} alone", ca, wa, None, d); rf = record(f"Fleet + ARMO {dn}", cf, wf, None, d)
            if d == "half":                     # gap statistics averaged over the same 8 starting positions
                ga, gf = [], []
                for ph in PHASES:
                    L = mk(ph); ta = gap_table(L, st.wban); tf = gap_table(merge(L, FLEET["Sentinel-1 + NISAR"]), st.wban)
                    ga.append(gap_summary(ta)); gf.append(gap_summary(tf))
                    if ph == 0:
                        gaps_all.append(ta.drop(columns="gaps").assign(system=f"ARMO {dn} alone (phase 0)"))
                ra.update(pd.DataFrame(ga).mean().to_dict()); rf.update(pd.DataFrame(gf).mean().to_dict())
            out.append(ra); out.append(rf)
            if d == "half":
                per_event[f"cap_ARMO {dn} alone"], per_event[f"wait_ARMO {dn} alone"] = ca, np.nanmean(wa, 0)
                per_event[f"cap_Fleet + ARMO {dn}"], per_event[f"wait_Fleet + ARMO {dn}"] = cf, np.nanmean(wf, 0)
                pooled_waits[f"ARMO {dn} alone"], pooled_waits[f"Fleet + ARMO {dn}"] = wa, wf
        print("done definition", d, flush=True)

    summ = pd.DataFrame(rows); bd = pd.DataFrame(by_def)
    gain = per_event["cap_Fleet + ARMO D"] - per_event["cap_Sentinel-1 + NISAR"].astype(float)
    lo, hi = boot_ci(ev, gain.to_numpy())
    print(f"\npaired gain, fleet + ARMO D over fleet: {100 * gain.mean():.1f} points (95 % CI {100 * lo:.1f}-{100 * hi:.1f})")
    summ.to_csv(os.path.join(PROC, "capture_summary.csv"), index=False)
    pd.concat([summ, bd], ignore_index=True).to_csv(os.path.join(PROC, "capture_by_definition.csv"), index=False)
    per_event.to_parquet(os.path.join(PROC, "capture_events.parquet"))
    np.savez_compressed(os.path.join(PROC, "pooled_waits.npz"), **{k.replace(" ", "_").replace("+", "plus"): v for k, v in pooled_waits.items()})
    pd.concat(gaps_all, ignore_index=True).to_csv(os.path.join(PROC, "gaps_station.csv"), index=False)
    pd.set_option("display.width", 260); pd.set_option("display.max_columns", 30)
    cols = ["system", "capture", "ci_lo", "ci_hi", "wait_median_h", "wait_mean_h", "wait_p90_h", "wait_le_24h",
            "mean_gap_d", "mean_of_gaps_d", "share_stations_lt_1p4", "mean_max_gap_d", "worst_gap_d"]
    print(summ[cols].round(3).to_string(index=False))
    print("\nby definition:")
    allr = pd.concat([summ, bd], ignore_index=True)
    print(allr.pivot_table(index="system", columns="definition", values="capture").round(3).to_string())
