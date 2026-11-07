"""Step 13 - extra numbers for the report's method questions, from this independent re-run.
  1. paired difference, two satellites (design D) minus the best single satellite (design C), alone and with today's fleet
  2. Sentinel-1 + ARMO without NISAR (NISAR's three-year prime mission may end before ARMO flies)
  3. the gap distribution for ARMO alone: mean, median, 90th percentile, longest
  4. the 50 % choice itself: windows at 25 / 37 / 50 / 75 % of the rise -> median window, regular-sampler capture at 1, 2,
     3 and 6 days, the revisit needed for 80 % capture, and fleet / ARMO capture on the 2026 cohort
  5. capture counts for today's fleet (events with a Sentinel-1 look, a NISAR look, either, neither)
Output: data/processed/review_checks.csv (also printed)"""
import os, sys
import numpy as np, pandas as pd
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from armo_common import (PROC, PHASES, load_events, stations_for, fleet_looks, constellation, merge, score, gap_table, boot_ci,
                         phase_average)

rows = []
def put(item, value, note=""):
    rows.append(dict(item=item, value=value, note=note)); print(f"{item:70} {value}  {note}")

if __name__ == "__main__":
    ev = load_events(); st = stations_for(ev); FL = fleet_looks(); fleet = FL["Sentinel-1 + NISAR"]
    mkD = lambda ph: constellation(st, n_sats=2, phase=ph, theta=(25, 45), side="right")
    mkC = lambda ph: constellation(st, n_sats=1, phase=ph, theta=(20, 50), side="alt")   # team's design C
    # 1. paired D - C
    cD, _ = phase_average(ev, mkD); cC, _ = phase_average(ev, mkC)
    fD, _ = phase_average(ev, mkD, base=fleet); fC, _ = phase_average(ev, mkC, base=fleet)
    for lab, a, b in (("alone", cD, cC), ("with today's fleet", fD, fC)):
        d = a - b; lo, hi = boot_ci(ev, d)
        put(f"design D minus design C, {lab} (points)", f"{100 * d.mean():.1f} (95 % CI {100 * lo:.1f} to {100 * hi:.1f})",
            f"D {100 * a.mean():.1f} %, C {100 * b.mean():.1f} %")
    # 2. without NISAR
    s1 = FL["Sentinel-1"]
    sD, _ = phase_average(ev, mkD, base=s1); c_s1 = score(ev, s1)[0].astype(float)
    lo, hi = boot_ci(ev, sD)
    put("Sentinel-1 + ARMO (no NISAR) capture", f"{100 * sD.mean():.1f} % (95 % CI {100 * lo:.1f}-{100 * hi:.1f})",
        f"Sentinel-1 alone {100 * c_s1.mean():.1f} %")
    # 3. ARMO gap distribution (pooled over 8 starting positions)
    gaps = np.concatenate([np.concatenate([g for g in gap_table(mkD(ph), st.wban).gaps if len(g)]) for ph in PHASES])
    put("ARMO D alone gaps: mean / median / 90th pct / longest (days)",
        f"{gaps.mean():.2f} / {np.median(gaps):.2f} / {np.percentile(gaps, 90):.2f} / {gaps.max():.2f}")
    fg = np.concatenate([g for g in gap_table(fleet, st.wban).gaps if len(g)])
    put("today's fleet gaps: mean / median / 90th pct / longest (days)",
        f"{fg.mean():.2f} / {np.median(fg):.2f} / {np.percentile(fg, 90):.2f} / {fg.max():.2f}")
    # 4. sweep of the threshold fraction
    evall = pd.read_parquet(os.path.join(PROC, "events.parquet"))
    for d, lab in (("q25", "25 %"), ("efold", "37 % (1/e)"), ("half", "50 % (report)"), ("q75", "75 %"), ("abs", "absolute 0.03")):
        W = (evall[f"W_{d}_h"] / 24.0).dropna().to_numpy()
        cap = lambda R: np.minimum(1.0, W / R).mean()
        r80 = brentq(lambda R: cap(R) - 0.8, 0.05, 30)
        cf = score(ev, fleet, d)[0].astype(float).mean(); ca, _ = phase_average(ev, mkD, d); cb, _ = phase_average(ev, mkD, d, base=fleet)
        put(f"window rule {lab}: median window (h) | regular capture 1/2/3/6 d | revisit for 80 %",
            f"{np.median(W) * 24:.0f} h | {100 * cap(1):.1f} / {100 * cap(2):.1f} / {100 * cap(3):.1f} / {100 * cap(6):.1f} % | {r80:.2f} d",
            f"2026 cohort: today {100 * cf:.1f} %, ARMO alone {100 * ca.mean():.1f} %, today + ARMO {100 * cb.mean():.1f} %")
    # 5. counts
    a = score(ev, FL["Sentinel-1"])[0]; b = score(ev, FL["NISAR"])[0]
    put("today's fleet counts: S1 / NISAR / either / S1 only / NISAR only / both / none",
        f"{a.sum()} / {b.sum()} / {(a | b).sum()} / {(a & ~b).sum()} / {(b & ~a).sum()} / {(a & b).sum()} / {(~a & ~b).sum()}",
        f"of {len(ev)} events")
    gained = (phase_average(ev, mkD, base=fleet)[0] - (a | b)).sum()
    put("events gained by adding ARMO D (average over starting positions)", f"{gained:.0f}")
    pd.DataFrame(rows).to_csv(os.path.join(PROC, "review_checks.csv"), index=False)
