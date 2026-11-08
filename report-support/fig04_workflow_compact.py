"""Report Figure 4 (Section 3.2): the analysis workflow, with the three phases, the objectives they serve and the
validation steps, laid out to print legibly at 15-16 cm width."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

W, H = 16.0 / 2.54, 5.6 / 2.54
fig = plt.figure(figsize=(W, H), dpi=400)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 160); ax.set_ylim(0, 56); ax.axis("off")
plt.rcParams["font.family"] = "DejaVu Sans"
cols = [(1.5, "Phase 1  Data assembly", "#eaf2fb"), (55.0, "Phase 2  Analysis", "#fdf3e6"), (108.5, "Phase 3  Results", "#ebf6ee")]
cw = 50.0
for x, title, bg in cols:
    ax.add_patch(FancyBboxPatch((x, 9.2), cw, 45.8, boxstyle="round,pad=0,rounding_size=1.5", fc=bg, ec="none"))
    ax.text(x + 2, 52.2, title, fontsize=7.4, fontweight="bold", va="center")
rows_y = [41.6, 31.2, 20.8, 10.4]
rh = [8.2, 8.2, 8.2, 8.2]
def box(x, y, h, t1, t2, ec="#8a8a8a"):
    ax.add_patch(FancyBboxPatch((x + 1.5, y), cw - 3, h, boxstyle="round,pad=0,rounding_size=1.0", fc="white", ec=ec, lw=0.7))
    ax.text(x + 3.2, y + h - 2.3, t1, fontsize=6.9, fontweight="bold", va="center")
    ax.text(x + 3.2, y + 2.4, t2, fontsize=6.5, va="center", color="#333333")
p1 = [("USCRN hourly ground truth", "121 stations, Apr–Oct 2025–26, QC"),
      ("Real look times (observed)", "Sentinel-1 STAC, NISAR CMR passes"),
      ("Earth Engine, ancillary data", "SMAP, gridMET, WorldCover, CDL")]
p2 = [("Wetting events (O1)", "≥ 5 mm storm, ≥ 0.03 m³/m³ rise"),
      ("Real-fleet capture (O2)", "≥ 1 look in window; 979 events"),
      ("Orbit simulator (O3)", "J2, 620 km, 25–45°; designs A–E"),
      ("Revisit, maps, sizing (O3–O4)", "1–12-day sampler; 0.1° maps")]
p3 = [("O1  Windows are short", "3,340 events; median 1.8 d"),
      ("O2  The fleet misses half", "S1 30.2, NISAR 37.9, either 51.6 %"),
      ("O3  Requirement and design", "80 % needs ≈ 1.5 d; D: 81.5 %"),
      ("O4  Where, and at what cost", "need maps; Table 4; US$ 115 M")]
for i, (a, b) in enumerate(p1):
    box(cols[0][0], rows_y[i], rh[i], a, b)
for i, (a, b) in enumerate(p2):
    box(cols[1][0], rows_y[i], rh[i], a, b)
for i, (a, b) in enumerate(p3):
    box(cols[2][0], rows_y[i], rh[i], a, b, ec="#2e8b57")
def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=5, lw=0.6, color="#555555"))
xa, xb = cols[0][0] + cw - 1.5, cols[1][0] + 1.5
arrow(xa, rows_y[0] + 4.3, xb, rows_y[0] + 4.3)                 # USCRN -> O1
arrow(xa, rows_y[1] + 4.3, xb, rows_y[1] + 4.3)                 # looks -> O2
arrow(xa, rows_y[1] + 3.0, xb, rows_y[2] + 5.0)                 # looks -> orbit sim
arrow(xa, rows_y[2] + 4.3, xb, rows_y[3] + 3.6)                 # EE -> revisit/maps
for i in range(4):                                             # phase 2 -> phase 3
    arrow(cols[1][0] + cw - 1.5, rows_y[i] + rh[i] / 2, cols[2][0] + 1.5, rows_y[i] + rh[i] / 2)
ax.add_patch(FancyBboxPatch((1.5, 1.2), 157, 6.6, boxstyle="round,pad=0,rounding_size=1.2", fc="#f3f2ee", ec="none"))
ax.text(4, 4.5, "Validation", fontsize=7.0, fontweight="bold", va="center")
ax.text(22, 4.5, "bootstrap (2,000 draws) · paired gains · held-out map stations · 3 × 3 thresholds · 10 cm probes · re-run of O1–O3",
        fontsize=6.2, va="center", color="#333333")
import os
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(out, exist_ok=True)
fig.savefig(os.path.join(out, "fig04_workflow_compact.png"), dpi=400)
print("ok")
