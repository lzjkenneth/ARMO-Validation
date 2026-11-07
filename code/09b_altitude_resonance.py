"""Step 9b - altitude traps: capture and mean gap against altitude in 5 km steps (2 satellites 180 deg apart, dawn-dusk,
right-looking), with the ground swath held at 362 km so only the orbit's repeat pattern changes.

Why: a sun-synchronous orbit whose revolutions per day are close to a simple fraction (15/1, 29/2, 44/3 ...) repeats its
ground track after a few days, so its swaths stack on the same strips and stations between them are rarely or never seen.
Output: data/processed/altitude_scan.csv, figures/fig_altitude_scan.png"""
import os, sys, math
from fractions import Fraction
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from armo_common import PROC, FIG, orb, load_events, stations_for, fleet_looks, constellation, merge, score, gap_table, gap_summary

PH = [0, 45, 90, 135, 180, 225, 270, 315]


def revs_per_day(h_km):
    a = orb.RE + h_km; n = math.sqrt(orb.MU / a ** 3); inc = math.radians(orb.sso_inclination(h_km))
    k = 0.75 * orb.J2 * (orb.RE / a) ** 2
    u_dot = n * (1 + k * (6 - 8 * math.sin(inc) ** 2))
    raan_dot = -1.5 * n * orb.J2 * (orb.RE / a) ** 2 * math.cos(inc)
    return u_dot / (orb.OMEGA_E - raan_dot)          # nodal revolutions per (nodal) day


def theta_for(h_km, swath_km, g_mid):
    from scipy.optimize import brentq
    inv = lambda g: brentq(lambda th: orb.ground_range_km(h_km, th) - g, 0.01, 75.0)
    return inv(g_mid - swath_km / 2), inv(g_mid + swath_km / 2)


if __name__ == "__main__":
    ev = load_events(); st = stations_for(ev); FLEET = fleet_looks()["Sentinel-1 + NISAR"]
    g_mid = 0.5 * (orb.ground_range_km(620, 25) + orb.ground_range_km(620, 45))
    rows = []
    for h in range(500, 805, 5):
        th = theta_for(float(h), 362.0, g_mid)
        caps, capf, gs, waits = [], [], [], []
        for ph in PH:
            L = constellation(st, n_sats=2, phase=ph, h_km=float(h), theta=th)
            c, w = score(ev, L); caps.append(c.mean()); waits.append(w)
            capf.append(score(ev, merge(L, FLEET))[0].mean()); gs.append(gap_summary(gap_table(L, st.wban)))
        q = revs_per_day(h); fr = Fraction(q).limit_denominator(10)
        g = pd.DataFrame(gs).mean()
        rows.append(dict(h_km=h, revs_per_day=q, near_fraction=f"{fr.numerator}/{fr.denominator}",
                         repeat_days_if_exact=fr.denominator, off_by=abs(q - float(fr)), theta1=th[0], theta2=th[1],
                         capture_alone=np.mean(caps), capture_alone_min_phase=np.min(caps), capture_alone_max_phase=np.max(caps),
                         capture_with_fleet=np.mean(capf), mean_gap_d=g.mean_gap_d, worst_gap_d=g.worst_gap_d,
                         share_stations_lt_1p4=g.share_stations_lt_1p4, median_wait_h=np.nanmedian(np.ravel(waits))))
        r = rows[-1]
        print(f"{h} km  {q:7.4f} rev/day (~{r['near_fraction']:>6})  alone {100 * r['capture_alone']:5.1f}  "
              f"+fleet {100 * r['capture_with_fleet']:5.1f}  mean gap {r['mean_gap_d']:.2f} d  worst {r['worst_gap_d']:.1f} d", flush=True)
    out = pd.DataFrame(rows); out.to_csv(os.path.join(PROC, "altitude_scan.csv"), index=False)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})
    fig, ax = plt.subplots(2, 1, figsize=(9, 5.6), sharex=True)
    ax[0].plot(out.h_km, 100 * out.capture_alone, "o-", ms=3, color="#1f4e79", label="ARMO alone (2 satellites)")
    ax[0].fill_between(out.h_km, 100 * out.capture_alone_min_phase, 100 * out.capture_alone_max_phase, color="#1f4e79", alpha=0.12,
                       label="range over 8 starting positions")
    ax[0].plot(out.h_km, 100 * out.capture_with_fleet, "s-", ms=3, color="#27ae60", label="Sentinel-1 + NISAR + ARMO")
    ax[0].axvline(620, color="0.4", lw=0.8, ls=":"); ax[0].text(622, 45, "620 km\n(proposed)", fontsize=8)
    ax[0].set_ylabel("events captured (%)"); ax[0].legend(frameon=False, fontsize=8, loc="lower right")
    ax[0].set_title("Same 362 km swath at every altitude: capture depends on the orbit's repeat pattern", fontsize=9)
    ax[1].plot(out.h_km, out.worst_gap_d, "o-", ms=3, color="#c0392b", label="worst gap at any station")
    ax[1].plot(out.h_km, out.mean_gap_d, "o-", ms=3, color="#1f4e79", label="mean gap between looks")
    for _, r in out[out.off_by < 0.004].iterrows():
        ax[1].annotate(r.near_fraction, (r.h_km, r.worst_gap_d), textcoords="offset points", xytext=(0, 4), fontsize=7, ha="center")
    ax[1].set_yscale("log"); ax[1].set_ylabel("days"); ax[1].set_xlabel("orbit altitude (km)"); ax[1].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_altitude_scan.png")); plt.close(fig)
    print("figure written")
