"""Step 12 - does a shorter gap between the "before" and "after" images lower the change noise? (real Sentinel-1 data)

ARMO's argument is not only more looks but shorter same-geometry intervals: the scene (crops, tillage, dew) changes less in a
few days than in 12. NISAR only gives 12-day repeats, but in 2026 Sentinel-1C and 1D fly 180 deg apart on the same tracks,
so many US stations have same-track pairs 6 days apart as well as 12 days. C-band is more sensitive to vegetation than
L-band, so this is a cautious stand-in for ARMO, not a measurement of it.

Data: Sentinel-1 RTC gamma0 (10 m, VV and VH) from Microsoft Planetary Computer (anonymous token; CC BY 4.0), a 5.12 km
window around each station, Apr-Oct 2026. Pairs = same relative orbit, 6, 12, 18 or 24 days apart. Stable = USCRN 5 cm
soil-moisture change below 0.01 m3/m3 and no rain in the 3 h before either pass. Same noise measures as step 8.
Outputs: data/raw/s1_windows/ (not in Git), data/processed/s1_noise_pairs.csv, s1_noise_by_lag.csv,
figures/fig_s1_noise_vs_interval.png"""
import os, sys, json, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROC = os.path.join(ROOT, "data", "processed")
RAW = os.path.join(ROOT, "data", "raw", "s1_windows"); FIG = os.path.join(ROOT, "figures")
STAC = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
TOKEN_URL = "https://planetarycomputer.microsoft.com/api/sas/v1/token/sentinel-1-rtc"
HALF = 256
os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif,.tiff")
_tok = {"t": None, "at": 0.0}


def token():
    if _tok["t"] is None or time.time() - _tok["at"] > 600:          # tokens last about an hour; renew every 10 min
        r = json.loads(urllib.request.urlopen(urllib.request.Request(TOKEN_URL, headers={"User-Agent": "ES4304"}), timeout=60).read())
        _tok.update(t=r["token"], at=time.time())
    return _tok["t"]


def items(lon, lat):
    body = {"collections": ["sentinel-1-rtc"], "intersects": {"type": "Point", "coordinates": [lon, lat]},
            "datetime": "2026-04-01T00:00:00Z/2026-11-01T00:00:00Z", "limit": 200}
    out = []
    while True:
        req = urllib.request.Request(STAC, data=json.dumps(body).encode(), headers={"Content-Type": "application/json", "User-Agent": "ES4304"})
        j = json.loads(urllib.request.urlopen(req, timeout=120).read())
        for f in j["features"]:
            a = f["assets"]
            if "vv" not in a:
                continue
            out.append(dict(id=f["id"], t=f["properties"]["datetime"], rel_orbit=f["properties"].get("sat:relative_orbit"),
                            vv=a["vv"]["href"], vh=a.get("vh", {}).get("href")))
        nxt = next((l for l in j.get("links", []) if l.get("rel") == "next"), None)
        if not nxt:
            return out
        body = nxt.get("body", body)


def read_window(job):
    import rasterio
    from rasterio.windows import Window
    from pyproj import Transformer
    wban, lat, lon, it = job
    dst = os.path.join(RAW, str(wban), it["id"] + ".npz")
    if os.path.exists(dst):
        return "cached"
    try:
        out = {}
        for pol in ("vv", "vh"):
            if not it.get(pol):
                continue
            with rasterio.open(it[pol] + "?" + token()) as src:
                X, Y = Transformer.from_crs(4326, src.crs, always_xy=True).transform(lon, lat)
                r, c = src.index(X, Y)
                if not (HALF <= r < src.height - HALF and HALF <= c < src.width - HALF):
                    return "near edge"
                win = Window(c - HALF, r - HALF, 2 * HALF, 2 * HALF)
                a = src.read(1, window=win).astype(np.float32)
                if src.nodata is not None:
                    a[a == src.nodata] = np.nan
                out[pol] = a
                if pol == "vv":
                    tr = src.window_transform(win); out.update(x0=tr.c, y0=tr.f, dx=tr.a, dy=tr.e, crs=str(src.crs))
        if "vv" not in out or not np.isfinite(out["vv"][HALF - 5:HALF + 5, HALF - 5:HALF + 5]).any():
            return "no data"
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        np.savez_compressed(dst, vv=out["vv"], vh=out.get("vh", np.full_like(out["vv"], np.nan)), x0=out["x0"], y0=out["y0"],
                            dx=out["dx"], dy=out["dy"], crs=out["crs"], t=it["t"][:19], rel_orbit=it["rel_orbit"] or -1)
        return "ok"
    except Exception as e:
        return "error " + str(e)[:80]


if __name__ == "__main__":
    s1 = pd.read_csv(os.path.join(PROC, "looks_s1.csv"), parse_dates=["time"])
    s1 = s1[(s1.time >= "2026-04-01") & (s1.time < "2026-11-01")].sort_values(["wban", "rel_orbit", "time"])
    s1["dt"] = s1.groupby(["wban", "rel_orbit"]).time.diff().dt.total_seconds().div(86400).round()
    cnt = s1.groupby("wban").dt.agg(n6=lambda x: (x == 6).sum(), n12=lambda x: (x == 12).sum())
    pick = cnt[(cnt.n6 >= 10) & (cnt.n12 >= 5)].sort_values("n6", ascending=False).index[:int(os.environ.get("N_STATIONS", "20"))]
    st = pd.read_csv(os.path.join(PROC, "station_seasons.csv")).drop_duplicates("wban").set_index("wban")
    jobs = []
    for w in pick:
        for it in items(st.loc[w, "lon"], st.loc[w, "lat"]):
            jobs.append((int(w), float(st.loc[w, "lat"]), float(st.loc[w, "lon"]), it))
    print(len(pick), "stations,", len(jobs), "Sentinel-1 RTC windows", flush=True)
    res = []
    with ThreadPoolExecutor(8) as ex:
        for k, r in enumerate(ex.map(read_window, jobs), 1):
            res.append(r)
            if k % 100 == 0:
                print(k, pd.Series([x[:8] for x in res]).value_counts().to_dict(), flush=True)
    print("done", pd.Series([x[:8] for x in res]).value_counts().to_dict(), flush=True)
