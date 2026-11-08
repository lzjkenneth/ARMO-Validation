"""Print-size redraws of Figures 6, 9 and 10 of the ARMO report.

The team's PNGs were drawn about 19.5 cm wide and shrunk to 10-12.5 cm in the report, so their smallest labels printed at
3-4 pt. These redraws use the team's own output tables (figures/generated/*.csv on origin/main of the team repository, read
with `git show` from a clone named in TEAM_REPO, or from the copies in report-support/data/team_generated/) and are drawn at the printed width, dpi 400, with every text element
at least 6.5 pt and axis labels 7-8 pt. Data, panels, panel letters, error bars, ordering and colours are the team's
(figstyle.py); only the layout is new. The small-print footnotes under the plots are dropped (the captions carry them).

Outputs (report-support/figures/): fig06_print.png, fig09_print.png, fig10_print.png
Run:  python report-support/fig_print_redraws.py [6] [9] [10]

Provenance of the one number that is not in a CSV: the 95 % interval of the observed-fleet bar in Figure 9a (47.3-55.7 %)
is not written to any figures/generated table (fig_capture_gain_paired_ci.py computes it on the fly from event-level data
that is not in the repository). It is read from the team's committed 300-dpi figS12_capture_gain_paired_ci.png, calibrated
against the known interval caps of the other two bars (see measure_fleet_ci). Cross-check: a station-cluster bootstrap on
the independent re-run's event table gives 47.2-55.3 %.
"""
import io
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.patheffects as pe
from matplotlib.legend import Legend
from matplotlib.text import Text
from matplotlib.ticker import NullLocator

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
CACHE = os.path.join(HERE, "data", "team_generated")
REPO = os.environ.get("TEAM_REPO", "")   # optional local clone of https://github.com/zheniannn/ARMO
REF = "origin/main"
DPI = 400
CM = 1 / 2.54

# ---- team colours (figures/figstyle.py on origin/main) -------------------------------------------------------------
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
# colours the team scripts set inline
SME2, SMAPREC = "#a99be8", "#b0476f"                       # fig_quality_attrition.py
GOLD_EDGE, PALE_GOLD, SAND = "#b07a00", "#f3d38a", "#e8cf98"   # fig_capture_gain_paired_ci.py
CYCLED = "#e8b04a"                                         # fig_swath_capture_sensitivity.py

# ---- type sizes (pt, at printed size) --------------------------------------------------------------------------------
FS_TICK = 6.5      # tick labels, legends, annotations, bar labels (the floor)
FS_AXLAB = 7.5     # axis labels (7-8 pt)
FS_TITLE = 7.5     # panel titles
FS_LETTER = 9.5    # bold panel letters
LW_AX = 0.5
HALO = [pe.withStroke(linewidth=1.8, foreground="white")]    # keeps a reference line from running through a label
TOP = 2.5          # clearance of the first text row to the top border (pt)
BOTTOM = 2.5       # clearance of the last text row to the bottom border (pt)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": FS_TICK, "axes.labelsize": FS_AXLAB, "axes.titlesize": FS_TITLE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.linewidth": LW_AX, "xtick.major.width": LW_AX, "ytick.major.width": LW_AX,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2, "xtick.major.pad": 1.8, "ytick.major.pad": 1.8,
    "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
    "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.45,
    "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": FS_TICK, "lines.linewidth": 1.0,
    "savefig.dpi": DPI, "figure.dpi": 100, "axes.unicode_minus": True,
})


# =====================================================================================================================
# data access
# =====================================================================================================================
def team_csv(name):
    """figures/generated/<name> of the team repository (origin/main): cached copy, else `git show` (read only)."""
    path = os.path.join(CACHE, name)
    if not os.path.exists(path):
        raw = subprocess.run(["git", "-C", REPO, "show", f"{REF}:figures/generated/{name}"],
                             capture_output=True, check=True).stdout
        os.makedirs(CACHE, exist_ok=True)
        with open(path, "wb") as fh:
            fh.write(raw)
    return pd.read_csv(path)


# =====================================================================================================================
# layout helpers: everything is placed in points from the top-left corner of the printed figure
# =====================================================================================================================
class Canvas:
    def __init__(self, w_cm, h_pt):
        self.w = w_cm * CM * 72.0
        self.h = float(h_pt)
        self.fig = plt.figure(figsize=(w_cm * CM, self.h / 72.0), dpi=DPI)
        self.fig.patch.set_facecolor("white")

    def ax(self, left, top, width, height):
        return self.fig.add_axes([left / self.w, 1 - (top + height) / self.h, width / self.w, height / self.h])

    def text(self, x, y, s, **kw):
        return self.fig.text(x / self.w, 1 - y / self.h, s, **kw)

    def letter(self, x, y_top, s):
        """bold panel letter; (x, y_top) = left edge and top of the letter, in pt"""
        return self.text(x, y_top, s, fontsize=FS_LETTER, fontweight="bold", ha="left", va="top")


def hide_grid_spines(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def right_text(ax, y, s, dx=4.0, **kw):
    """text in a column just right of the axes, vertically at data y"""
    kw.setdefault("fontsize", FS_TICK)
    return ax.annotate(s, xy=(1.0, y), xycoords=("axes fraction", "data"), xytext=(dx, 0), textcoords="offset points",
                       ha="left", va="center", annotation_clip=False, **kw)


# =====================================================================================================================
# audit: font sizes, overlaps between text elements, text outside the canvas
# =====================================================================================================================
def _texts(fig):
    out = []
    seen = set()

    def add(t, kind):
        if id(t) in seen or not isinstance(t, Text) or not t.get_visible() or not t.get_text().strip():
            return
        seen.add(id(t))
        out.append((t, kind))

    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for tick in axis._update_ticks():
                add(tick.label1, "tick")
                add(tick.label2, "tick")
            add(axis.label, "axislabel")
        add(ax.title, "title")
        add(ax._left_title, "title")
        add(ax._right_title, "title")
        for t in ax.texts:
            add(t, "text")
    for t in fig.texts:
        add(t, "text")
    for lg in fig.findobj(Legend):
        for t in lg.get_texts():
            add(t, "legend")
    return out


def audit(fig, name):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    items = _texts(fig)
    wpx, hpx = fig.bbox.width, fig.bbox.height
    boxes = [(t, kind, t.get_window_extent(r)) for t, kind in items]
    sizes = [(t.get_fontsize(), kind, t.get_text().replace("\n", " / ")) for t, kind, _ in boxes]
    smin = min(s for s, _, _ in sizes)
    axl = sorted({round(t.get_fontsize(), 2) for t, kind, _ in boxes if kind == "axislabel"})
    problems = []
    pad = 0.35 * DPI / 72.0 * 0  # no tolerance on canvas bounds
    for t, kind, bb in boxes:
        if bb.x0 < -0.5 or bb.y0 < -0.5 or bb.x1 > wpx + 0.5 or bb.y1 > hpx + 0.5:
            problems.append(f"OUTSIDE canvas: {kind} '{t.get_text()[:40]}' bbox=({bb.x0:.0f},{bb.y0:.0f},{bb.x1:.0f},{bb.y1:.0f}) canvas=({wpx:.0f},{hpx:.0f})")
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][2], boxes[j][2]
            ox = min(a.x1, b.x1) - max(a.x0, b.x0)
            oy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ox > 0.8 and oy > 0.8:
                problems.append(f"OVERLAP {ox * 72 / DPI:.1f}x{oy * 72 / DPI:.1f} pt: [{boxes[i][1]}] '{boxes[i][0].get_text()[:32]}' vs "
                                f"[{boxes[j][1]}] '{boxes[j][0].get_text()[:32]}'")
    w_cm, h_cm = fig.get_figwidth() / CM, fig.get_figheight() / CM
    k = 72.0 / DPI
    sides = {"left": min(bb.x0 for _, _, bb in boxes) * k, "right": min(wpx - bb.x1 for _, _, bb in boxes) * k,
             "top": min(hpx - bb.y1 for _, _, bb in boxes) * k, "bottom": min(bb.y0 for _, _, bb in boxes) * k}
    print(f"[{name}] canvas {w_cm:.2f} x {h_cm:.2f} cm; {len(boxes)} text elements; smallest font {smin:.2f} pt; axis-label sizes {axl}")
    print("      clearance of the text to the border (pt): " + ", ".join(f"{a} {v:.1f}" for a, v in sides.items()))
    for a, v in sides.items():
        if v < 2.0:
            problems.append(f"text within {v:.1f} pt of the {a} border")
    for s, kind, txt in sorted(sizes)[:3]:
        print(f"      smallest: {s:.2f} pt [{kind}] '{txt[:50]}'")
    if problems:
        for p in problems:
            print("   !!", p)
    else:
        print("   no text overlaps, nothing outside the canvas")
    return smin, problems


def save(cv, fname):
    """write the figure at its printed size (no bbox_inches="tight": the canvas IS the printed size) as an RGB PNG, dpi 400"""
    from PIL import Image
    path = os.path.join(FIGDIR, fname)
    buf = io.BytesIO()
    cv.fig.savefig(buf, format="png", dpi=DPI, facecolor="white")
    buf.seek(0)
    with Image.open(buf) as im:
        im = im.convert("RGB")
        im.save(path, dpi=(DPI, DPI), optimize=True)
        wpx, hpx = im.size
    print(f"      wrote {path}: {wpx} x {hpx} px = {wpx / DPI * 2.54:.2f} x {hpx / DPI * 2.54:.2f} cm at {DPI} dpi")
    return wpx, hpx


def content_bottom_pt(fig):
    """lowest text element (pt from the top of the canvas) - used to trim the canvas height to the content"""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    hpx = fig.bbox.height
    low = 0.0
    for t, _ in _texts(fig):
        bb = t.get_window_extent(r)
        low = max(low, (hpx - bb.y0) * 72.0 / DPI)
    return low


def fit_height(build, margin=BOTTOM, **kw):
    """build(H_PT=...) twice: once tall, then with the canvas trimmed to the lowest text plus `margin` pt"""
    tall = build(H_PT=400.0, **kw)
    h = content_bottom_pt(tall.fig) + margin
    plt.close(tall.fig)
    return build(H_PT=h, **kw)


# =====================================================================================================================
# Figure 6: from a real look to a detected change (fig_quality_attrition.py, quality_attrition_stats.csv)
# =====================================================================================================================
SENS = ["S1 VV", "NISAR HH", "NISAR SME2", "SMAP all quality", "SMAP recommended"]
SENS_COL = {"S1 VV": ORANGE, "NISAR HH": VIOLET, "NISAR SME2": SME2, "SMAP all quality": MAGENTA, "SMAP recommended": SMAPREC}
STAGES = ["Real look in window", "Valid same- geometry pair", "Valid retrieval in pair", "Change > 1σ of no-change noise"]
STAGE_LABELS = ["Real look\nin window", "Valid same-\ngeometry pair", "Valid retrieval\nin pair", "Change > 1σ of\nno-change noise"]


def fig6(W_CM=11.8, H_PT=214.0, title_fs=7.0):
    q = team_csv("quality_attrition_stats.csv")
    qa = q[(q.panel == "a") & q.stage.isin(STAGES)]
    union = q[q.stage == "detected change (union)"].iloc[0]
    N = int(union.denominator)
    bq = q[q.panel == "b"].set_index("sensor")
    cq = q[q.panel == "c"].set_index("sensor")

    cv = Canvas(W_CM, H_PT)
    W = cv.w
    # ---- vertical layout (pt from the top)
    T1 = TOP                      # title row top
    LEG = T1 + 10.5               # legend row top
    a_top, a_h = LEG + 10.5, 69.5 # (a) axes
    bc_title = a_top + a_h + 4.0 + 2 * 7.5 + 4.0      # (b)/(c) title row top (after the 2-line tick labels of (a))
    row = 11.3
    bc_top = bc_title + 11.0

    # ------------------------------------------------------------------ (a)
    a_l = 41.0
    a = cv.ax(a_l, a_top, W - a_l - 4.0, a_h)
    x = np.arange(4)
    bw = 0.16
    for j, s in enumerate(SENS):
        r = qa[qa.sensor == s].set_index("stage").loc[STAGES]
        n = r.n.to_numpy()
        assert (r.denominator == N).all()
        xx = x + (j - 2) * bw
        a.bar(xx, n / N * 100, bw * 0.92, color=SENS_COL[s], linewidth=0, zorder=3)
        for xi, k in zip(xx, n):
            a.annotate(str(int(k)), (xi, k / N * 100), xytext=(0, 1.4), textcoords="offset points", ha="center", va="bottom",
                       rotation=90, fontsize=FS_TICK, color=INK, annotation_clip=False)
    a.set_xlim(-0.55, 3.55)
    a.set_ylim(0, 100)
    a.set_yticks(range(0, 101, 20))
    a.set_xticks(x)
    a.set_xticklabels(STAGE_LABELS, fontsize=FS_TICK, linespacing=1.15)
    a.set_ylabel(f"% of the {N}\nmatched events", fontsize=FS_AXLAB, labelpad=2, linespacing=1.1)
    a.grid(axis="y", color=GRID, lw=0.45, zorder=0)
    a.text(0.995, 0.97, "numbers above\nbars = events", transform=a.transAxes, ha="right", va="top", fontsize=FS_TICK,
           color=INK2, linespacing=1.15)
    cv.letter(3.0, T1, "a")
    cv.text(14.0, T1 + 0.6, f"Same {N} events: a radar-detected change in only {int(union.n)} "
                            f"({union.value * 100:.1f} %, any field-scale radar)", fontsize=title_fs, ha="left", va="top")
    handles = [Patch(facecolor=SENS_COL[s], edgecolor="none", label=s) for s in SENS]
    cv.fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(a_l / W, 1 - LEG / cv.h), ncol=5, handlelength=0.9,
                  handleheight=0.9, handletextpad=0.4, columnspacing=1.0, borderpad=0, borderaxespad=0, fontsize=FS_TICK)

    # ------------------------------------------------------------------ (b) and (c)
    b_l, b_w = 77.0, 42.0
    c_l, c_w = 218.0, 47.0
    y = np.arange(5)[::-1]
    b = cv.ax(b_l, bc_top, b_w, 5 * row)
    c = cv.ax(c_l, bc_top, c_w, 5 * row)
    for ax, xmax in ((b, 76.0), (c, 13.5)):
        ax.set_ylim(-0.5, 4.5)
        ax.set_xlim(0, xmax)
        ax.set_yticks(y)
        ax.grid(axis="x", color=GRID, lw=0.45, zorder=0)
        ax.spines["left"].set_visible(True)
    b.set_yticklabels(SENS, fontsize=FS_TICK)
    c.set_yticklabels([])
    c.tick_params(axis="y", length=0)
    b.set_xticks([0, 25, 50, 75])
    c.set_xticks([0, 5, 10])
    for yi, s in zip(y, SENS):
        r = bq.loc[s]
        p, lo, hi = r.value, r.ci_lo, r.ci_hi
        b.barh(yi, p * 100, height=0.6, color=SENS_COL[s], linewidth=0, zorder=3)
        b.errorbar(p * 100, yi, xerr=[[(p - lo) * 100], [(hi - p) * 100]], fmt="none", ecolor=INK, capsize=1.6, lw=0.7,
                   capthick=0.7, zorder=4)
        period = r.cohort.split(";")[0]
        right_text(b, yi, f"{p * 100:.0f} % (n={int(r.denominator)}, {period})")
        r = cq.loc[s]
        p, lo, hi = r.value, r.ci_lo, r.ci_hi
        c.barh(yi, p * 100, height=0.6, color=SENS_COL[s], alpha=0.55, linewidth=0, zorder=3)
        c.errorbar(p * 100, yi, xerr=[[(p - lo) * 100], [(hi - p) * 100]], fmt="none", ecolor=INK, capsize=1.6, lw=0.7,
                   capthick=0.7, zorder=4)
        right_text(c, yi, f"{p * 100:.1f} % (n={int(r.denominator):,})")
    b.axvline(70, color=BLUE, ls=":", lw=0.9, zorder=5)
    c.axvline(10, color=BLUE, ls=":", lw=0.9, zorder=5)
    b.set_xlabel("Detected when sampled (%)\n(dotted: ARMO 70 % flight target)", fontsize=FS_AXLAB, labelpad=1.5, linespacing=1.1)
    c.set_xlabel("False alarms on unchanged pairs (%)\n(dotted: 10 % target)", fontsize=FS_AXLAB, labelpad=1.5, linespacing=1.1)
    cv.letter(3.0, bc_title, "b")
    cv.text(14.0, bc_title + 0.6, "Detection, published denominators", fontsize=title_fs, ha="left", va="top")
    cv.letter(c_l - 8.0, bc_title, "c")
    cv.text(c_l + 3.0, bc_title + 0.6, "False alarms, own denominator", fontsize=title_fs, ha="left", va="top")
    audit(cv.fig, "fig06")
    return cv


# =====================================================================================================================
# Figure 10: usable-strip sensitivity (fig_swath_capture_sensitivity.py, swath_sensitivity_stats.csv)
# =====================================================================================================================
def fig10(W_CM=13.1, H_PT=168.0, title_fs=7.0):
    R = team_csv("swath_sensitivity_stats.csv")
    P = team_csv("capture_gain_paired_stats.csv")
    fleet = float(P.fleet_capture.iloc[0]) * 100            # observed fleet, 51.6 %
    n = len(R)
    assert n == 11 and R["mode"].iloc[0] == "baseline"
    labels = ["362 km (baseline)" if m == "baseline" else sc for sc, m in zip(R.scenario, R["mode"])]

    cv = Canvas(W_CM, H_PT)
    W = cv.w
    p = 9.4                                    # row pitch (pt)
    BAND = 2                                   # two rows above the bars for the notes / legend
    a_l, a_w = 66.0, 96.0
    b_l, b_w = a_l + a_w + 12.0, 88.0
    c_l, c_w = b_l + 106.0, 88.0
    top_axes = TOP + 11.5
    h_axes = (n + BAND) * p
    y = np.arange(n)[::-1].astype(float)       # 10 ... 0
    ylim = (-0.5, n - 0.5 + BAND)
    yb1, yb2 = n - 0.5 + 0.5, n - 0.5 + 1.5    # band row centres (11, 12)

    def style(ax, xlim, ticks):
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xticks(ticks)
        ax.set_yticks(y)
        ax.vlines(ticks, -0.5, n - 0.5, color=GRID, lw=0.45, zorder=0)
        ax.spines["left"].set_bounds(-0.5, n - 0.5)
        ax.tick_params(axis="y", length=0)

    def band_y(yy):
        return (yy - ylim[0]) / (ylim[1] - ylim[0])

    a = cv.ax(a_l, top_axes, a_w, h_axes)
    b = cv.ax(b_l, top_axes, b_w, h_axes)
    c = cv.ax(c_l, top_axes, c_w, h_axes)
    style(a, (40, 100), [40, 60, 80, 100])
    style(b, (-27.5, 2.5), [-20, -10, 0])
    style(c, (1.0, 100.0), [1, 2, 5, 10, 20, 40])
    a.set_yticklabels(labels, fontsize=FS_TICK)
    b.set_yticklabels([])
    c.set_yticklabels([])
    c.set_xscale("log")
    c.set_xlim(1.0, 100.0)
    c.set_xticks([1, 2, 5, 10, 20, 40])
    c.set_xticklabels(["1", "2", "5", "10", "20", "40"])
    c.xaxis.set_minor_locator(NullLocator())
    c.vlines([1, 2, 5, 10, 20, 40], -0.5, n - 0.5, color=GRID, lw=0.45, zorder=0)

    for yi, r in zip(y, R.itertuples()):
        # (a) capture with the range over the 8 orbital phases
        col = YELLOW if r.mode == "baseline" else (CYCLED if r.mode == "cycled" else PALE_GOLD)
        a.barh(yi, r.combined * 100, height=0.6, color=col, edgecolor=GOLD_EDGE, linewidth=0.5, zorder=3)
        a.plot([r.combined_phase_min * 100, r.combined_phase_max * 100], [yi, yi], color=INK, lw=0.9, solid_capstyle="butt", zorder=5)
        a.annotate(f"{r.combined * 100:.1f} %", xy=(r.combined_phase_max * 100, yi), xytext=(2.4, 0), textcoords="offset points",
                   ha="left", va="center", fontsize=FS_TICK, annotation_clip=False, zorder=6, path_effects=HALO)
        # (b) paired change vs the 362 km baseline with its 95 % station-bootstrap CI
        base = r.mode == "baseline"
        b.errorbar(r.vs_baseline_pp, yi, xerr=[[r.vs_baseline_pp - r.vs_baseline_ci_lo], [r.vs_baseline_ci_hi - r.vs_baseline_pp]],
                   fmt="o", color=INK2 if base else RED, ms=3.2, capsize=1.8, lw=0.8, capthick=0.8, zorder=4)
        # (c) ARMO-only gaps: median and 90th percentile, longest gap at the median phase
        c.plot([r.gap_p50_d, r.gap_p90_d], [yi, yi], color=INK2, lw=1.0, zorder=3, solid_capstyle="butt")
        c.scatter([r.gap_p50_d], [yi], s=3.6 ** 2, color=YELLOW, edgecolor=INK, linewidth=0.5, zorder=4)
        c.scatter([r.gap_p90_d], [yi], s=3.6 ** 2, color=ORANGE, edgecolor=INK, linewidth=0.5, zorder=4)
        c.scatter([r.gap_max_phase_median_d], [yi], s=3.8 ** 2, marker="x", color=RED, linewidths=0.9, zorder=4)

    # (a) reference lines and their labels (band row 1), note on the black line (band row 2)
    top_line = n - 0.5 + 1.0
    a.plot([fleet, fleet], [-0.5, top_line], color=BLUE, lw=0.9, zorder=4, solid_capstyle="butt")
    a.plot([80, 80], [-0.5, top_line], color=MUTED, lw=0.8, ls=(0, (1.2, 1.6)), zorder=4)
    a.annotate(f"real fleet {fleet:.1f} %", xy=(fleet, yb1), xytext=(-2.4, 0), textcoords="offset points", ha="right", va="center",
               color=BLUE, fontsize=FS_TICK)
    a.annotate("80 %", xy=(80, yb1), xytext=(2.4, 0), textcoords="offset points", ha="left", va="center", color=INK2, fontsize=FS_TICK)
    a.annotate("line = range over 8 orbital phases", xy=(1.0, yb2), xycoords=("axes fraction", "data"), xytext=(0, 0),
               textcoords="offset points", ha="right", va="center", color=INK2, fontsize=FS_TICK)
    b.plot([0, 0], [-0.5, n - 0.5], color=MUTED, lw=0.7, zorder=2)
    b.annotate("95 % station-\nbootstrap CI", xy=(0.5, (yb1 + yb2) / 2), xycoords=("axes fraction", "data"), ha="center", va="center",
               color=INK2, fontsize=FS_TICK, linespacing=1.1)
    h1 = [Line2D([], [], marker="o", ls="", mfc=YELLOW, mec=INK, mew=0.5, ms=3.6, label="median"),
          Line2D([], [], marker="o", ls="", mfc=ORANGE, mec=INK, mew=0.5, ms=3.6, label="90th pct")]
    h2 = [Line2D([], [], marker="x", ls="", mec=RED, mew=0.9, ms=3.8, label="longest (median phase)")]
    kw = dict(loc="center left", borderpad=0, borderaxespad=0, handlelength=0.9, handletextpad=0.3, labelspacing=0, frameon=False,
              fontsize=FS_TICK, columnspacing=0.9)
    leg1 = c.legend(handles=h1, ncol=2, bbox_to_anchor=(0.0, band_y(yb2)), **kw)
    c.add_artist(leg1)
    c.legend(handles=h2, ncol=1, bbox_to_anchor=(0.0, band_y(yb1)), **kw)

    a.set_xlabel("Fleet + ARMO D:\nevents with a look (%)", fontsize=FS_AXLAB, labelpad=1.5, linespacing=1.1)
    b.set_xlabel("Paired change vs\nbaseline (pp)", fontsize=FS_AXLAB, labelpad=1.5, linespacing=1.1)
    c.set_xlabel("ARMO-only gap\nbetween looks (d)", fontsize=FS_AXLAB, labelpad=1.5, linespacing=1.1)
    ty = TOP + 0.6
    cv.letter(3.0, TOP, "a")
    cv.text(14.0, ty, "Look opportunity vs usable strip", fontsize=title_fs, ha="left", va="top")
    cv.letter(b_l - 11.0, TOP, "b")
    cv.text(b_l, ty, "Loss from a narrower strip", fontsize=title_fs, ha="left", va="top")
    cv.letter(c_l - 11.0, TOP, "c")
    cv.text(c_l, ty, "ARMO-only revisit gaps", fontsize=title_fs, ha="left", va="top")
    audit(cv.fig, "fig10")
    return cv


# =====================================================================================================================
# Figure 9: capture gain from ARMO design D (fig_capture_gain_paired_ci.py, capture_gain_paired_stats.csv,
#           capture_gain_by_duration_D.csv)
# =====================================================================================================================
# Two things in this figure are not written to any figures/generated table and are read from the team's committed 300-dpi
# figS12_capture_gain_paired_ci.png (verify with: fig_print_redraws.py --verify-png):
#  1. the 95 % station-bootstrap interval of the observed-fleet bar in panel (a), in %: 47.33-55.74 (the bar itself, 51.6 %, and the
#     other two intervals are in capture_gain_paired_stats.csv; the pixel reading is calibrated on those known interval caps);
#  2. the eight per-phase gains behind the dots of panel (b). The CSV keeps only their minimum and maximum (gain_phase_min/max);
#     the dot positions below were decoded from the PNG (fit of equal discs to the outline of each dot cluster, extremes set to
#     the CSV minimum and maximum). Dots that lie under another dot in the team's PNG cannot be told apart and are not repeated
#     (A 8, B 6, C 5, D 3, E 6 visible dots); the hidden ones needed to reproduce the CSV mean gain fall inside the covered range.
FLEET_CI_PCT = (47.33, 55.74)
PHASE_GAINS_PP = {
    "A": [18.08, 18.28, 18.93, 19.09, 19.59, 20.23, 20.51, 21.14],
    "B": [17.98, 18.39, 18.48, 18.66, 19.13, 19.31],
    "C": [25.13, 25.66, 25.83, 26.19, 26.35],
    "D": [29.01, 29.74, 30.44],
    "E": [25.23, 25.51, 25.79, 25.85, 26.38, 26.56],
}
PHASE_DOT_PT = 2.1       # diameter of a phase dot (pt); the team's dots are 0.61 pp wide, about 1.8 pt at this scale
DESIGN_NAMES = {"A": "A 1 sat", "B": "B 1 sat roll", "C": "C 1 sat wide", "D": "D 2 sats", "E": "E 2 sats roll"}


def fig9(W_CM=11.0, H_PT=201.0, title_fs=7.0):
    P = team_csv("capture_gain_paired_stats.csv").set_index("design")
    DR = team_csv("capture_gain_by_duration_D.csv")
    D = P.loc["D"]
    vals = np.array([D.fleet_capture, D.combined_capture, D.armo_alone]) * 100
    ci_lo = np.array([FLEET_CI_PCT[0], D.combined_ci_lo * 100, D.armo_alone_ci_lo * 100])
    ci_hi = np.array([FLEET_CI_PCT[1], D.combined_ci_hi * 100, D.armo_alone_ci_hi * 100])

    cv = Canvas(W_CM, H_PT)
    W = cv.w
    # ------------------------------------------------------------------ geometry (pt from the top-left)
    a_l, a_w, a_top, a_h = 31.0, 124.0, TOP + 11.5, 148.0
    XR = a_l + a_w + 11.0                      # left edge of the right-hand column
    R_right = W - 6.0
    b_l = XR + 44.0
    c_l = XR + 36.5
    b_top, b_h = TOP + 19.5, 45.0
    c_title = b_top + b_h + 2.2 + 1.8 + 6.6 + 1.5 + 8.5 + 3.5
    c_leg = c_title + 16.1 + 1.8
    c_top = c_leg + 7.3 + 2.5
    c_h = 50.0

    # ------------------------------------------------------------------ (a)
    a = cv.ax(a_l, a_top, a_w, a_h)
    x = np.array([0.0, 1.0, 2.15])
    cols = [BLUE, YELLOW, PALE_GOLD]
    edges = [BLUE, GOLD_EDGE, GOLD_EDGE]
    a.bar(x, vals, width=0.72, color=cols, edgecolor=edges, linewidth=0.8, zorder=3)
    a.errorbar(x, vals, yerr=[vals - ci_lo, ci_hi - vals], fmt="none", ecolor=INK, capsize=2.2, lw=0.8, capthick=0.8, zorder=4)
    for xi, v in zip(x, vals):
        a.text(xi, 3.5, f"{v:.1f} %", ha="center", va="bottom", fontsize=7.0, fontweight="bold", color=INK, zorder=5)
    a.annotate("", xy=(1, vals[1] + 6), xytext=(0, vals[0] + 6),
               arrowprops=dict(arrowstyle="->", color=INK, lw=1.0, connectionstyle="arc3,rad=-0.3", shrinkA=0, shrinkB=0,
                               mutation_scale=7), zorder=5)
    a.text(0.97, 91.5, f"paired +{D.gain_pp:.1f} pp\n95 % CI {D.gain_ci_station_lo:.1f}–{D.gain_ci_station_hi:.1f} (stations)\n"
                       f"{D.gain_ci_station_and_phase_lo:.1f}–{D.gain_ci_station_and_phase_hi:.1f} (stations + phase)",
           ha="center", va="bottom", fontsize=FS_TICK, linespacing=1.15)
    a.set_xlim(-0.48, 2.62)
    a.set_ylim(0, 118)
    a.set_yticks(range(0, 101, 20))
    a.set_xticks(x)
    a.set_xticklabels(["Real fleet\n(S1C/D +\nNISAR)\nOBSERVED", "Fleet +\nARMO D\nSIMULATED\nopportunity",
                       "ARMO D\nalone\nSIMULATED\nopportunity"], fontsize=FS_TICK, linespacing=1.15)
    a.grid(axis="y", color=GRID, lw=0.45, zorder=0)
    a.spines["left"].set_bounds(0, 118)
    a.set_ylabel("Events with ≥ 1 look in the window (%)", fontsize=FS_AXLAB, labelpad=2)
    cv.letter(3.0, TOP, "a")
    cv.text(14.0, TOP + 0.6, "Look opportunity, same 979 events", fontsize=title_fs, ha="left", va="top")

    # ------------------------------------------------------------------ (b)
    b = cv.ax(b_l, b_top, R_right - b_l, b_h)
    XLIM_B = (10.0, 42.2)
    yb = np.arange(5)[::-1].astype(float)
    for yi, d in zip(yb, "ABCDE"):
        r = P.loc[d]
        gain = r.gain_pp
        ph = np.array(PHASE_GAINS_PP[d])
        assert abs(ph.min() - r.gain_phase_min) < 0.006 and abs(ph.max() - r.gain_phase_max) < 0.006, d
        b.plot([r.gain_ci_station_and_phase_lo, r.gain_ci_station_and_phase_hi], [yi, yi], color=MUTED, lw=3.2, alpha=0.6,
               solid_capstyle="butt", zorder=2)
        b.plot([r.gain_ci_station_lo, r.gain_ci_station_hi], [yi, yi], color=INK, lw=0.9, zorder=3, solid_capstyle="butt")
        b.scatter(ph, np.full(ph.size, yi + 0.45), s=PHASE_DOT_PT ** 2, color=INK2, linewidths=0, zorder=3)
        b.scatter([gain], [yi], s=4.4 ** 2, color=YELLOW if d == "D" else SAND, edgecolor=INK, linewidth=0.6, zorder=4)
        b.annotate(f"+{gain:.1f}", xy=(1.0, yi), xycoords=("axes fraction", "data"), xytext=(0, 0), textcoords="offset points",
                   ha="right", va="center", fontsize=FS_TICK, annotation_clip=False, zorder=6, path_effects=HALO)
    b.set_xlim(*XLIM_B)
    b.set_ylim(-0.55, 4.7)
    b.set_yticks(yb)
    b.set_yticklabels([DESIGN_NAMES[d] for d in "ABCDE"], fontsize=FS_TICK)
    b.set_xticks(range(10, 41, 5))
    b.tick_params(axis="y", length=0)
    b.grid(axis="x", color=GRID, lw=0.45, zorder=0)
    b.set_xlabel("Gain over real fleet (pp)", fontsize=FS_AXLAB, labelpad=1.5)
    cv.letter(XR, TOP, "b")
    cv.text(XR + 11.0, TOP + 0.6, "All designs: CI (line),\n+ phase (grey), 8 phases (dots)", fontsize=title_fs, ha="left", va="top", linespacing=1.15)

    # ------------------------------------------------------------------ (c)
    c = cv.ax(c_l, c_top, R_right - c_l, c_h)
    xx = np.arange(len(DR), dtype=float)
    c.bar(xx - 0.18, DR.fleet * 100, 0.34, color=BLUE, linewidth=0, zorder=3)
    c.bar(xx + 0.18, DR.combined * 100, 0.34, color=YELLOW, linewidth=0, zorder=3)
    for xi, r in zip(xx, DR.itertuples()):
        c.annotate(f"{int(r.n)}", (xi + 0.18, 0), xytext=(0, 2.0), textcoords="offset points", ha="center", va="bottom", rotation=90,
                   fontsize=FS_TICK, color=INK, zorder=5)
    c.set_xlim(-0.55, len(DR) - 0.45)
    c.set_ylim(0, 104)
    c.set_yticks([0, 25, 50, 75, 100])
    c.set_xticks(xx)
    c.set_xticklabels(list(DR["bin"]), fontsize=FS_TICK)
    c.grid(axis="y", color=GRID, lw=0.45, zorder=0)
    c.set_ylabel("Events with\na look (%)", fontsize=FS_AXLAB, labelpad=2, linespacing=1.1)
    c.set_xlabel("Wetting-window duration\n(numbers in bars = events)", fontsize=FS_AXLAB, labelpad=1.5, linespacing=1.1)
    cv.letter(XR, c_title - 0.6, "c")
    cv.text(XR + 11.0, c_title, "Look opportunity by window\nduration (design D)", fontsize=title_fs, ha="left", va="top", linespacing=1.15)
    hs = [Patch(facecolor=BLUE, edgecolor="none", label="Real fleet"), Patch(facecolor=YELLOW, edgecolor="none", label="+ ARMO D (simulated)")]
    cv.fig.legend(handles=hs, loc="upper left", bbox_to_anchor=((XR + 11.0) / W, 1 - c_leg / cv.h), ncol=2, handlelength=0.9,
                  handleheight=0.9, handletextpad=0.4, columnspacing=1.2, borderpad=0, borderaxespad=0, fontsize=FS_TICK)
    audit(cv.fig, "fig09")
    return cv


if __name__ == "__main__":
    which = [a for a in sys.argv[1:]] or ["6", "9", "10"]
    os.makedirs(FIGDIR, exist_ok=True)
    if "6" in which:
        save(fit_height(fig6), "fig06_print.png")
    if "9" in which:
        save(fit_height(fig9), "fig09_print.png")
    if "10" in which:
        save(fit_height(fig10), "fig10_print.png")
