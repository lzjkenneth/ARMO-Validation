"""Station-by-station mean look gap of design D (two satellites 180 deg apart, one plane, right-looking) for the
proposed 620 km orbit and the exact-repeat options of repeat_orbit_options.csv, 25-45 deg off-nadir unless stated.
Re-run orbit model and 2026 cohort (code/armo_common.py), 12 starting phases. Writes per-station gaps and a summary to
report-support/data/. Inputs: the committed files in data/processed/. Run from the repository root:
python report-support/gap_by_orbit.py (a few minutes)."""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "code"))
from armo_common import load_events, stations_for, fleet_looks, constellation, merge, score, gap_table, orb

ev = load_events(); st = stations_for(ev); FLEET = fleet_looks()["Sentinel-1 + NISAR"]
PH = range(0, 360, 30)
configs = [("620 km (proposed, no exact repeat)", 620.0, (25.0, 45.0)),
           ("623.5 km (74/5, 5-day repeat)", 623.5, (25.0, 45.0)),
           ("617.8 km (163/11, 11-day repeat)", 617.8, (25.0, 45.0)),
           ("613 km (89/6, 6-day repeat)", 613.0, (25.0, 45.0)),
           ("605.5 km (104/7, 7-day repeat)", 605.5, (25.0, 45.0)),
           ("599.9 km (119/8, 8-day repeat)", 599.9, (25.0, 45.0)),
           ("620 km, 24-46 deg", 620.0, (24.0, 46.0)),
           ("620 km, 23-46 deg", 620.0, (23.0, 46.0))]
rows, per_station = [], {}
for name, h, th in configs:
    per = []; ca = []; cf = []
    for ph in PH:
        L = constellation(st, n_sats=2, phase=ph, h_km=h, theta=th)
        per.append(gap_table(L, st.wban).set_index("wban")["mean_gap_d"])
        ca.append(score(ev, L)[0].mean()); cf.append(score(ev, merge(L, FLEET))[0].mean())
    P = pd.concat(per, axis=1)            # stations x phases
    g = P.mean(axis=1)                     # phase-averaged station mean gap
    worst_phase = P.max(axis=1)            # worst phase for each station
    la = st.set_index("wban").loc[g.index, "lat"]
    per_station[name] = g
    r = dict(config=name, h_km=h, off_nadir=f"{th[0]:.0f}-{th[1]:.0f}", access_km=round(orb.access_width_km(h, *th)),
             capture_alone=round(100 * np.mean(ca), 1), capture_with_fleet=round(100 * np.mean(cf), 1),
             conus_mean_gap_d=round(g.mean(), 3), p90_station_gap_d=round(g.quantile(0.9), 2), max_station_gap_d=round(g.max(), 2),
             stations_over_1p52=int((g > 1.52).sum()), stations_over_1p75=int((g > 1.75).sum()),
             worst_phase_max_gap_d=round(worst_phase.max(), 2), mean_gap_south_35N=round(g[la < 35].mean(), 2),
             mean_gap_north_35N=round(g[la >= 35].mean(), 2))
    rows.append(r); print(r, flush=True)
os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
pd.DataFrame(rows).to_csv(os.path.join(HERE, "data", "gap_by_orbit.csv"), index=False)
pd.DataFrame(per_station).assign(lat=st.set_index("wban")["lat"]).to_csv(os.path.join(HERE, "data", "gap_by_orbit_stations.csv"))
