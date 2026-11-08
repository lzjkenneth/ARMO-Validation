"""Validation figure: station mean look gap against latitude for design D (two satellites), with the 1.5-day Level-1 line.
Inputs: report-support/data/gap_by_orbit_stations.csv (from gap_by_orbit.py) and data/processed/gap_by_station.csv
(fleet + ARMO, from code/05_capture.py). Output: report-support/figures/val_gap_by_latitude.png.
Run from the repository root: python report-support/fig_gap_latitude.py"""
import os
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
S = pd.read_csv(os.path.join(HERE, "data", "gap_by_orbit_stations.csv"), index_col=0)
C = pd.read_csv(os.path.join(ROOT, "data", "processed", "gap_by_station.csv")).set_index("wban")
fig, ax = plt.subplots(figsize=(7.2, 3.9), dpi=200)
series = [("620 km (proposed, no exact repeat)", "620 km, 25–45° (proposed)", "#1f4e79", "o"),
          ("613 km (89/6, 6-day repeat)", "613 km exact repeat, 25–45°", "#c0504d", "s"),
          ("620 km, 24-46 deg", "620 km, 24–46° (wider access)", "#2e8b57", "^")]
for col, lab, c, m in series:
    ax.scatter(S["lat"], S[col], s=13, marker=m, color=c, alpha=0.8, label=f"{lab}: CONUS mean {S[col].mean():.2f} d, max {S[col].max():.2f} d")
ax.scatter(C.loc[S.index, "lat"], C.loc[S.index, "comb_gap"], s=10, marker="x", color="0.45",
           label=f"Existing fleet + ARMO (620 km): max {C.loc[S.index, 'comb_gap'].max():.2f} d")
ax.axhline(1.5, color="k", lw=1, ls="--"); ax.text(48.6, 1.52, "Level-1: 1.5 days", ha="right", va="bottom", fontsize=7.5)
ax.set_xlabel("Station latitude (°N), 107 USCRN stations"); ax.set_ylabel("Mean gap between ARMO looks (days)")
ax.set_xlim(25, 49); ax.set_ylim(0.6, 2.45); ax.grid(alpha=0.25, lw=0.5)
ax.legend(fontsize=6.6, loc="upper right", frameon=True, framealpha=0.95)
ax.set_title("Revisit gap by station: two ARMO satellites, 2026 season, 12 orbital phases averaged", fontsize=8.5)
out = os.path.join(HERE, "figures", "val_gap_by_latitude.png"); os.makedirs(os.path.dirname(out), exist_ok=True)
fig.tight_layout(); fig.savefig(out); print("wrote", out)
