"""Step 9 - design sweeps: what actually buys capture and a short waiting time?
Answers the questions from the presentation: why two satellites and not one; what if the usable imaging swath is smaller
than the 362 km access; does altitude or time of day matter; what does a wider swath or a second look side buy.

Every configuration is flown at 12 starting positions (30 deg steps) and the results averaged. Satellites in one plane are
equally spaced. "Usable swath" = a fixed band of ground range, centred on the middle of the 25-45 deg access at 620 km,
that the radar images on every pass (a realistic systematic mapping plan); 362 km reproduces the 25-45 deg access.
Scored on the same 2026 cohort (half-amplitude windows), alone and merged with the real Sentinel-1 + NISAR looks.
Output: data/processed/design_sweeps.csv"""
import os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from armo_common import PROC, orb, load_events, stations_for, fleet_looks, constellation, merge, score, gap_table, gap_summary

PH12 = list(range(0, 360, 30))
H0, ACC = 620.0, (25.0, 45.0)
G_MID = 0.5 * (orb.ground_range_km(H0, ACC[0]) + orb.ground_range_km(H0, ACC[1]))


def theta_for(h_km, swath_km):
    """Off-nadir band of ground width swath_km centred on the access midpoint (in ground range at 620 km, kept fixed in
    ground range for other altitudes)."""
    from scipy.optimize import brentq
    g1, g2 = G_MID - swath_km / 2, G_MID + swath_km / 2
    inv = lambda g: brentq(lambda th: orb.ground_range_km(h_km, th) - g, 0.01, 75.0)
    return (inv(max(g1, 1.0)), inv(g2))


def evaluate(ev, st, fleet, label, sweep, n_sats=2, planes=None, phases=PH12, **kw):
    caps_a, caps_f, waits_a, waits_f, gsum = [], [], [], [], []
    for ph in phases:
        L = constellation(st, n_sats=n_sats, phase=ph, planes=planes, **kw)
        ca, wa = score(ev, L); cf, wf = score(ev, merge(L, fleet))
        caps_a.append(ca); caps_f.append(cf); waits_a.append(wa); waits_f.append(wf)
        gsum.append(gap_summary(gap_table(L, st.wban)))
    ca, cf = np.mean(caps_a, 0), np.mean(caps_f, 0); wa, wf = np.array(waits_a).ravel(), np.array(waits_f).ravel()   # pooled
    g = pd.DataFrame(gsum).mean()
    th = kw.get("theta", ACC); h = kw.get("h_km", H0)
    return dict(sweep=sweep, label=label, n_sats=n_sats if planes is None else len(planes), h_km=h, ltan_h=kw.get("ltan_h", 18.0),
                theta1=th[0], theta2=th[1], swath_km=orb.access_width_km(h, *th), side=kw.get("side", "right"),
                capture_alone=ca.mean(), capture_with_fleet=cf.mean(), wait_median_alone_h=np.nanmedian(wa),
                wait_median_with_fleet_h=np.nanmedian(wf), wait_le_24h_alone=np.nanmean(wa <= 24),
                wait_le_24h_with_fleet=np.nanmean(wf <= 24), **g.to_dict())


if __name__ == "__main__":
    ev = load_events(); st = stations_for(ev); FLEET = fleet_looks()["Sentinel-1 + NISAR"]
    fleet_cap = score(ev, FLEET)[0].mean()
    rows = []; t0 = time.time()
    def run(*a, **k):
        r = evaluate(ev, st, FLEET, *a, **k); rows.append(r)
        print(f"{r['sweep']:10} {r['label']:42} alone {100 * r['capture_alone']:5.1f} | +fleet {100 * r['capture_with_fleet']:5.1f} | "
              f"mean gap {r['mean_gap_d']:.2f} d | median wait {r['wait_median_alone_h']:5.1f} h | {time.time() - t0:5.0f}s", flush=True)

    # 1. number of satellites in one plane (620 km, 25-45 deg right-looking)
    for n in (1, 2, 3, 4):
        run(f"{n} satellite(s), one plane", "n_sats", n_sats=n, theta=ACC)
    # 2. usable imaging swath, 1 and 2 satellites
    for n in (1, 2):
        for sw in ((150, 200, 250, 300, 362, 450, 574) if n == 1 else (150, 200, 250, 300, 362, 400, 420, 440, 450, 574)):
            run(f"{n} sat, usable swath {sw} km", "swath", n_sats=n, theta=theta_for(H0, sw))
    # 3. altitude with the same 25-45 deg off-nadir range (swath grows with altitude), 2 satellites
    for h in (450, 500, 560, 620, 700, 800):
        run(f"2 sats at {h} km, 25-45 deg", "altitude", n_sats=2, h_km=float(h), theta=ACC)
    # 4. time of day (local time of the ascending node), 2 satellites
    for lt, lab in ((18.0, "dawn-dusk, asc 18:00 (proposed)"), (6.0, "dawn-dusk, asc 06:00"), (22.5, "asc 22:30 (looks ~10:30/22:30)"),
                    (13.5, "asc 13:30 (looks ~01:30/13:30)")):
        run(f"2 sats, {lab}", "time_of_day", n_sats=2, ltan_h=lt, theta=ACC)
    # 5. look side, 2 satellites
    for side in ("left", "alt", "both"):
        run(f"2 sats, look side {side}", "look_side", n_sats=2, theta=ACC, side=side)
    # 6. two planes: second satellite in a different orbit plane / time of day
    for lt2, lab in ((22.5, "2nd sat in a 10:30/22:30 plane"), (13.5, "2nd sat in a 01:30/13:30 plane")):
        run(f"2 sats, {lab}", "two_planes", planes=[(0.0, dict(ltan_h=18.0)), (0.0, dict(ltan_h=lt2))], theta=ACC)
    # 7. one satellite alternatives the professor asked about
    run("1 sat, alternate sides, 25-45 deg", "one_sat", n_sats=1, theta=ACC, side="alt")
    run("1 sat, both sides every pass (bound)", "one_sat", n_sats=1, theta=ACC, side="both")
    run("1 sat, 20-50 deg (574 km access)", "one_sat", n_sats=1, theta=(20.0, 50.0))
    run("1 sat at 800 km, 20-50 deg", "one_sat", n_sats=1, h_km=800.0, theta=(20.0, 50.0))
    out = pd.DataFrame(rows); out["fleet_capture"] = fleet_cap
    out.to_csv(os.path.join(PROC, "design_sweeps.csv"), index=False)
    print("written", len(out), "rows; fleet alone captures", round(100 * fleet_cap, 1), "%")
