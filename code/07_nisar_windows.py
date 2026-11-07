"""Step 7 - real L-band data at the stations: read 5.12 x 5.12 km windows of NISAR GCOV (RTC gamma0, 10 m, HH and HV)
centred on USCRN stations, for every provisional pass 17 Jun-31 Oct 2026. Used by 08_noise_vs_resolution.py to measure how
backscatter-change noise falls as pixels are averaged to coarser product resolution, i.e. "how fine can ARMO's soil-moisture
product be?". Reads only the HDF5 chunks it needs from ASF (Earthdata login in a netrc file readable only by you; never printed).
Output: data/raw/nisar_windows/<wban>/<granule>.npz (raw imagery stays out of Git), data/processed/nisar_windows_index.csv"""
import os, sys, json, re, time, urllib.request, urllib.parse
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "data", "raw", "nisar_windows")
COLL = "NISAR_L2_GCOV_PROVISIONAL_V1"; HALF = 256
N_STATIONS = int(os.environ.get("N_STATIONS", "16"))


def granules(lon, lat):
    q = urllib.parse.urlencode({"short_name": COLL, "point": f"{lon},{lat}", "page_size": 500,
                                "temporal": "2026-06-17T00:00:00Z,2026-11-01T00:00:00Z"})
    r = urllib.request.urlopen(urllib.request.Request("https://cmr.earthdata.nasa.gov/search/granules.json?" + q,
                                                      headers={"User-Agent": "ES4304-student-validation"}), timeout=120)
    return [e["title"] for e in json.loads(r.read())["feed"]["entry"]]


def read_one(job):
    """Read the one 512 x 512 storage tile (10 m pixels, 5.12 km) that contains the station, HH and HV only: one HDF5
    chunk per layer keeps each read to a few range requests. The tile's UTM corner is saved so pairs can be co-registered."""
    sys.path.insert(0, os.path.join(ROOT, "code"))
    import h5py, pyproj
    from nisar_io import HTTPRange, signed_url
    wban, lat, lon, g = job
    dst = os.path.join(OUT, str(wban), g + ".npz")
    if os.path.exists(dst):
        return dict(wban=wban, granule=g, status="cached")
    try:
        fo = HTTPRange(signed_url(COLL, g)); f = h5py.File(fo, "r")
        grp = f["science/LSAR/GCOV/grids/frequencyA"]
        epsg = int(grp["projection"].attrs["epsg_code"])
        x01 = grp["xCoordinates"][0:2]; y01 = grp["yCoordinates"][0:2]
        nx, ny = grp["xCoordinates"].shape[0], grp["yCoordinates"].shape[0]
        dx, dy = x01[1] - x01[0], y01[1] - y01[0]
        X, Y = pyproj.Transformer.from_crs(4326, epsg, always_xy=True).transform(lon, lat)
        c = int(round((X - x01[0]) / dx)); r = int(round((Y - y01[0]) / dy))
        if not (0 <= r < ny and 0 <= c < nx):
            f.close(); return dict(wban=wban, granule=g, status="outside frame")
        ch = grp["HHHH"].chunks[0]
        r0, c0 = (r // ch) * ch, (c // ch) * ch
        sl = (slice(r0, min(r0 + ch, ny)), slice(c0, min(c0 + ch, nx)))
        hh = grp["HHHH"][sl].astype(np.float32)
        rr, cc = r - r0, c - c0
        if hh.shape != (ch, ch) or not np.isfinite(hh[max(rr - 5, 0):rr + 5, max(cc - 5, 0):cc + 5]).any():
            f.close(); return dict(wban=wban, granule=g, status="no data at station or partial tile")
        hv = grp["HVHV"][sl].astype(np.float32) if "HVHV" in grp else np.full_like(hh, np.nan)
        ident = f["science/LSAR/identification"]
        meta = dict(track=int(ident["trackNumber"][()]), frame=int(ident["frameNumber"][()]),
                    direction=ident["orbitPassDirection"][()].decode(), look=ident["lookDirection"][()].decode(),
                    t=ident["zeroDopplerStartTime"][()].decode()[:19], epsg=epsg, x0=float(x01[0] + c0 * dx),
                    y0=float(y01[0] + r0 * dy), dx=float(dx), dy=float(dy), sr=rr, sc=cc)
        f.close()
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        np.savez_compressed(dst, hh=hh, hv=hv, **{k: np.array(v) for k, v in meta.items()})
        return dict(wban=wban, granule=g, status="ok", mb=round(fo.nread / 1e6, 1), **meta)
    except Exception as e:
        return dict(wban=wban, granule=g, status="error: " + str(e)[:100])


if __name__ == "__main__":
    from multiprocessing import Pool
    looks = pd.read_csv(os.path.join(PROC, "looks_nisar.csv"), parse_dates=["time"])
    looks = looks[(looks.time >= "2026-06-17") & (looks.time < "2026-11-01")]
    st = pd.read_csv(os.path.join(PROC, "station_seasons.csv")).drop_duplicates("wban").set_index("wban")
    n = looks.groupby("wban").size().sort_values(ascending=False)
    # one station per state first (spread), most NISAR looks first; then fill up to N_STATIONS
    pick, seen = [], set()
    for w in n.index:
        state = str(st.loc[w, "name"])[:2]
        if state not in seen and n[w] >= 15:
            pick.append(w); seen.add(state)
    pick = pick[:N_STATIONS]
    jobs = []
    for w in pick:
        gl = pd.DataFrame({"g": granules(st.loc[w, "lon"], st.loc[w, "lat"])})
        p_ = gl.g.str.split("_", expand=True)            # NISAR_L2_PR_GCOV_<cycle>_<track>_<A/D>_<frame>_...
        gl["track"], gl["dirn"], gl["frame"] = p_[5], p_[6], p_[7]
        top = gl.groupby(["track", "dirn"]).frame.agg(lambda x: x.value_counts().index[0])
        gl = gl[gl.frame == gl.apply(lambda r: top[(r.track, r.dirn)], axis=1)]
        for g in gl.g:
            jobs.append((int(w), float(st.loc[w, "lat"]), float(st.loc[w, "lon"]), g))
    print(len(pick), "stations,", len(jobs), "granule windows to read", flush=True)
    rows = []
    with Pool(12) as pool:
        for k, r in enumerate(pool.imap_unordered(read_one, jobs), 1):
            rows.append(r)
            if k % 20 == 0:
                print(k, "done |", pd.Series([x["status"][:12] for x in rows]).value_counts().to_dict(), flush=True)
    idx = pd.DataFrame(rows); idx.to_csv(os.path.join(PROC, "nisar_windows_index.csv"), index=False)
    print(idx.status.str[:14].value_counts().to_string()); print("MB read:", idx.get("mb", pd.Series(dtype=float)).sum())
