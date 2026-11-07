"""Step 13b - designs A, C and D when only a 220, 250 or 300 km strip is imaged on each pass (the strip centred on the
middle of each design's access). A: 1 satellite right-looking; C: 1 satellite left/right on alternate orbits; D: 2
satellites 180 deg apart, right-looking. 8 starting positions, 2026 cohort, half-amplitude windows.
Output: data/processed/usable_strip_designs.csv"""
import os, sys
import numpy as np, pandas as pd
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from armo_common import PROC, orb, load_events, stations_for, fleet_looks, constellation, phase_average, boot_ci

def band(h, access, width):
    gm = 0.5 * (orb.ground_range_km(h, access[0]) + orb.ground_range_km(h, access[1]))
    inv = lambda g: brentq(lambda th: orb.ground_range_km(h, th) - g, 0.01, 75.0)
    return inv(gm - width / 2), inv(gm + width / 2)

if __name__ == "__main__":
    ev = load_events(); st = stations_for(ev); fleet = fleet_looks()["Sentinel-1 + NISAR"]
    rows = []
    for name, n, access, side in (("A (1 sat, right)", 1, (25, 45), "right"), ("C (1 sat, alternate, 20-50 access)", 1, (20, 50), "alt"),
                                  ("D (2 sats, right)", 2, (25, 45), "right")):
        for width in (220, 250, 300):
            th = band(620.0, access, width)
            mk = lambda ph: constellation(st, n_sats=n, phase=ph, theta=th, side=side)
            ca, _ = phase_average(ev, mk); cf, _ = phase_average(ev, mk, base=fleet)
            lo, hi = boot_ci(ev, cf)
            rows.append(dict(design=name, usable_strip_km=width, off_nadir=f"{th[0]:.1f}-{th[1]:.1f}", capture_alone=100 * ca.mean(),
                             capture_with_fleet=100 * cf.mean(), ci_lo=100 * lo, ci_hi=100 * hi))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(PROC, "usable_strip_designs.csv"), index=False)
