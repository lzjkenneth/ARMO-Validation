"""Figure 8 of the rewrite: design D on one day over CONUS (map with scale bar and north arrow) and the Stillwater timeline.

Uses the independent orbit simulator of this repository (code/04_orbit_sim.py: circular J2 orbit, 620 km,
ascending node 18:00, right-looking 25-45 deg off-nadir) and the re-run's real look times and USCRN hourly data.
Satellite A flies at orbital phase 0 deg and satellite B at 180 deg.
Inputs: data/processed/ (events.parquet, looks_s1.csv, looks_nisar.csv; uscrn_hourly.parquet is written by
code/01_ground_truth.py and is not committed) and report-support/data/ne_50m_admin_1_states_provinces.geojson
(Natural Earth, public domain; download it first, see report-support/README.md).
Output: report-support/figures/fig08_design_d.png
Run from the repository root: python report-support/fig08_design_d.py
"""
import os, sys, json, math, importlib.util
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MPoly, FancyArrow
from pyproj import Transformer

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                                   # repository root
spec = importlib.util.spec_from_file_location("orbit", os.path.join(ROOT, "code", "04_orbit_sim.py"))
orb = importlib.util.module_from_spec(spec); spec.loader.exec_module(orb)
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(HERE, "figures", "fig08_design_d.png")
STATES = os.path.join(HERE, "data", "ne_50m_admin_1_states_provinces.geojson")

H_KM, LTAN, TH = 620.0, 18.0, (25.0, 45.0)
DAY0, DAY1 = "2026-08-22 00:00", "2026-08-23 00:00"
T0, T1 = "2026-08-10 00:00", "2026-08-30 00:00"
plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False, "font.family": "DejaVu Sans"})
to5070 = Transformer.from_crs("EPSG:4326", "EPSG:5070", always_xy=True)


def track_and_swath(phase_deg, t0, t1, step=20.0):
    a = orb.RE + H_KM; n = math.sqrt(orb.MU / a ** 3); inc = math.radians(orb.sso_inclination(H_KM))
    k = 0.75 * orb.J2 * (orb.RE / a) ** 2
    raan_dot = -1.5 * n * orb.J2 * (orb.RE / a) ** 2 * math.cos(inc); u_dot = n * (1 + k * (6 - 8 * math.sin(inc) ** 2))
    Ts = pd.Timestamp(t0, tz="UTC"); total = (pd.Timestamp(t1, tz="UTC") - Ts).total_seconds()
    # the look simulator starts its phase at T0 of the look run; keep the same epoch so phases match panel (b)
    E = pd.Timestamp(T0, tz="UTC")
    d0 = (E - orb.J2000).total_seconds() / 86400.0
    sun_ra = math.radians((280.460 + 0.9856474 * d0) % 360); raan0 = sun_ra + math.radians((LTAN - 12.0) * 15.0)
    gmst0 = math.radians((280.46061837 + 360.98564736629 * d0) % 360)
    off = (Ts - E).total_seconds()
    tt = off + np.arange(0.0, total + step, step)
    R, V = orb._state(tt, a, inc, raan0, raan_dot, math.radians(phase_deg), u_dot, gmst0)
    r = R / np.linalg.norm(R, axis=1)[:, None]
    c = np.cross(V, r); c /= np.linalg.norm(c, axis=1)[:, None]          # right of track (same test as the simulator)
    g1, g2 = [orb.ground_range_km(H_KM, x) / orb.RE_MEAN for x in TH]
    ll = lambda p: (np.degrees(np.arctan2(p[:, 1], p[:, 0])), np.degrees(np.arcsin(np.clip(p[:, 2], -1, 1))))
    sub = ll(r)
    e1 = ll(np.cos(g1) * r + np.sin(g1) * c); e2 = ll(np.cos(g2) * r + np.sin(g2) * c)
    return sub, e1, e2


def segments(lon, lat, box=(-135, -55, 15, 58)):
    inside = (lon > box[0]) & (lon < box[1]) & (lat > box[2]) & (lat < box[3])
    idx = np.nonzero(inside)[0]
    if len(idx) == 0:
        return []
    splits = np.nonzero(np.diff(idx) > 1)[0] + 1
    return [s for s in np.split(idx, splits) if len(s) > 3]


def main():
    fig = plt.figure(figsize=(16.0 / 2.54, 6.2 / 2.54), dpi=300)
    ax = fig.add_axes([0.0, 0.0, 0.47, 1.0]); ax.set_aspect("equal"); ax.axis("off")
    # states
    if not os.path.exists(STATES):
        sys.exit(f"missing {STATES}: download the Natural Earth file first (see report-support/README.md)")
    gj = json.load(open(STATES))
    for f in gj["features"]:
        p = f["properties"]
        if p.get("iso_a2") != "US" or p.get("postal") in ("AK", "HI"):
            continue
        g = f["geometry"]; polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        for poly in polys:
            x, y = to5070.transform(*np.array(poly[0]).T)
            ax.add_patch(MPoly(np.c_[x, y], closed=True, fc="#f2f0eb", ec="#9a9a9a", lw=0.3, zorder=1))
    cols = {0.0: "#2f6db3", 180.0: "#d9622b"}
    for ph, col in cols.items():
        (slon, slat), (l1, b1), (l2, b2) = track_and_swath(ph, DAY0, DAY1)
        for s in segments(slon, slat):
            x1, y1 = to5070.transform(l1[s], b1[s]); x2, y2 = to5070.transform(l2[s], b2[s])
            ax.add_patch(MPoly(np.r_[np.c_[x1, y1], np.c_[x2, y2][::-1]], closed=True, fc=col, alpha=0.22, lw=0, zorder=2))
            xs, ys = to5070.transform(slon[s], slat[s]); ax.plot(xs, ys, ls=(0, (3, 2)), lw=0.5, color=col, zorder=3)
    # Stillwater
    h = pd.read_parquet(os.path.join(PROC, "uscrn_hourly.parquet"), columns=["wban", "lon", "lat", "name"])
    st = h[h.name.str.contains("Stillwater", case=False, na=False)].drop_duplicates("wban").iloc[0]
    sx, sy = to5070.transform(st.lon, st.lat)
    ax.plot(sx, sy, marker="*", ms=6, color="k", zorder=5); ax.text(sx + 6e4, sy - 1.2e5, "Stillwater, OK", fontsize=6, zorder=5)
    ax.set_xlim(-2.40e6, 2.30e6); ax.set_ylim(0.25e6, 3.20e6)
    # scale bar 500 km and north arrow
    x0, y0 = -2.25e6, 0.42e6
    ax.add_patch(plt.Rectangle((x0, y0), 5e5, 3.5e4, fc="k", zorder=6)); ax.add_patch(plt.Rectangle((x0 + 2.5e5, y0), 2.5e5, 3.5e4, fc="w", ec="k", lw=0.4, zorder=6))
    ax.text(x0, y0 + 6e4, "0", fontsize=5.5, ha="center"); ax.text(x0 + 5e5, y0 + 6e4, "500 km", fontsize=5.5, ha="center")
    ax.annotate("N", xy=(-2.05e6, 1.10e6), xytext=(-2.05e6, 0.78e6), ha="center", va="center", fontsize=7, fontweight="bold",
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color="k"), zorder=6,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
    ax.plot([], [], color=cols[0.0], lw=4, alpha=0.4, label="Satellite A access swath")
    ax.plot([], [], color=cols[180.0], lw=4, alpha=0.4, label="Satellite B (180° behind A)")
    ax.plot([], [], color="0.3", lw=0.6, ls=(0, (3, 2)), label="Ground track")
    ax.legend(loc="lower right", fontsize=5.5, frameon=False, bbox_to_anchor=(1.0, 0.0))
    ax.text(-2.38e6, 3.12e6, "a", fontsize=8, fontweight="bold")
    ax.text(-2.15e6, 3.12e6, "22 Aug 2026: 620 km dawn–dusk orbit, 362 km access swaths", fontsize=6)
    # timeline at Stillwater
    w = int(st.wban)
    hs = pd.read_parquet(os.path.join(PROC, "uscrn_hourly.parquet"), filters=[("wban", "==", w)], columns=["time", "sm5", "rain"])
    hs["time"] = pd.to_datetime(hs.time, utc=True)
    hs = hs[(hs.time >= pd.Timestamp(T0, tz="UTC")) & (hs.time < pd.Timestamp(T1, tz="UTC"))]
    ev = pd.read_parquet(os.path.join(PROC, "events.parquet")); ev = ev[ev.wban == w]
    ev = ev[(ev.we_half > pd.Timestamp(T0, tz="UTC")) & (ev.ws_half < pd.Timestamp(T1, tz="UTC"))]
    s1 = pd.read_csv(os.path.join(PROC, "looks_s1.csv"), parse_dates=["time"]); s1 = s1[s1.wban == w]
    ni = pd.read_csv(os.path.join(PROC, "looks_nisar.csv"), parse_dates=["time"]); ni = ni[ni.wban == w]
    stn = pd.DataFrame(dict(lat=[st.lat], lon=[st.lon]))
    armo = {}
    for ph, lab in ((0.0, "ARMO A"), (180.0, "ARMO B")):
        hrs = orb.looks(stn, T0, T1, h_km=H_KM, ltan_h=LTAN, phase_deg=ph, theta=TH, side="right", step_s=30.0)[0]
        armo[lab] = pd.to_datetime(hrs * 3600.0, unit="s", utc=True)
    a1 = fig.add_axes([0.585, 0.47, 0.40, 0.43]); a2 = fig.add_axes([0.585, 0.14, 0.40, 0.27], sharex=a1)
    for _, e in ev.iterrows():
        for a_ in (a1, a2):
            a_.axvspan(e.ws_half, e.we_half, color="#9ecae1", alpha=0.55, lw=0, zorder=0)
    a1.plot(hs.time, hs.sm5, color="k", lw=0.8)
    a1.set_ylabel("5 cm soil moisture\n(m³/m³)", fontsize=5.5); a1.tick_params(labelbottom=False, labelsize=6)
    a1.text(0.0, 1.06, "b", transform=a1.transAxes, fontsize=8, fontweight="bold")
    a1.text(0.06, 1.06, "Stillwater 2 W, 10–29 Aug 2026 (shaded: wetting windows)", transform=a1.transAxes, fontsize=6)
    rows = [("Sentinel-1", s1.time, "#6a3d9a"), ("NISAR", ni.time, "#238b45"), ("ARMO A", armo["ARMO A"], cols[0.0]),
            ("ARMO B", armo["ARMO B"], cols[180.0])]
    report = {}
    for k, (lab, tt, col) in enumerate(rows):
        tt = pd.Series(pd.to_datetime(tt, utc=True))
        tt = tt[(tt >= pd.Timestamp(T0, tz="UTC")) & (tt < pd.Timestamp(T1, tz="UTC"))]
        a2.scatter(tt, np.full(len(tt), len(rows) - k), marker="|", s=40, lw=1.2, color=col)
        report[lab] = [str(x) for x in tt]
    a2.set_yticks(range(len(rows), 0, -1)); a2.set_yticklabels([r[0] for r in rows], fontsize=5.5)
    a2.set_ylim(0.4, len(rows) + 0.6); a2.spines["left"].set_visible(False); a2.tick_params(axis="y", length=0)
    import matplotlib.dates as mdates
    a2.xaxis.set_major_locator(mdates.DayLocator(bymonthday=[10, 14, 18, 22, 26])); a2.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    a2.tick_params(axis="x", labelsize=6); a1.set_xlim(pd.Timestamp(T0, tz="UTC"), pd.Timestamp(T1, tz="UTC"))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=300); plt.close(fig)
    # which looks fall inside each window
    out = []
    for _, e in ev.iterrows():
        inside = {lab: [t for t in tl if e.ws_half <= pd.Timestamp(t) <= e.we_half] for lab, tl in report.items()}
        out.append(dict(window=f"{e.ws_half} to {e.we_half} ({e.W_half_h:.0f} h)", **{k: len(v) for k, v in inside.items()}))
    print("station", w, st["name"])
    print(pd.DataFrame(out).to_string(index=False))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
