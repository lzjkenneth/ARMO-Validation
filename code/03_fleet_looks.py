"""Step 3 - real look times of the existing field-scale radar fleet over every USCRN station (metadata only, no login):
  Sentinel-1 IW GRD : Microsoft Planetary Computer STAC (collection sentinel-1-grd)
  NISAR L-band      : NASA CMR granule search, NISAR_L2_GCOV_BETA_V1, _PROVISIONAL_V1 and the validated NISAR_L2_GCOV_V1
  ALOS-2 L-band     : NASA CMR, ALOS2_L1_PSR2 (the ALOS-2 ScanSAR archive held at ASF), to document what L-band data exists
Adjacent frames of the same pass (< 1 h apart) count as one look.
Output: data/processed/looks_s1.csv, looks_nisar.csv, looks_alos2.csv"""
import os, json, time, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROC = os.path.join(ROOT, "data", "processed")
H = {"User-Agent": "ES4304-student-validation", "Content-Type": "application/json"}
T0, T1 = "2025-03-15T00:00:00Z", "2026-11-01T00:00:00Z"

def http(url, body=None):
    err = None
    for k in range(8):
        try:
            hdr = H if body else {"User-Agent": H["User-Agent"]}          # CMR rejects GETs that carry a JSON content type
            req = urllib.request.Request(url, data=json.dumps(body).encode() if body else None, headers=hdr)
            return json.loads(urllib.request.urlopen(req, timeout=180).read())
        except Exception as e:
            err = e; time.sleep(5 * (k + 1))
    raise RuntimeError(f"{url} -> {err}")

def s1(lon, lat):
    out, url, body = [], "https://planetarycomputer.microsoft.com/api/stac/v1/search", {
        "collections": ["sentinel-1-grd"], "intersects": {"type": "Point", "coordinates": [lon, lat]},
        "datetime": f"{T0}/{T1}", "limit": 250}
    while True:
        j = http(url, body)
        for f in j["features"]:
            pr = f["properties"]
            out.append(dict(time=pr["datetime"], platform=pr.get("platform"), orbit=pr.get("sat:orbit_state"),
                            rel_orbit=pr.get("sat:relative_orbit")))
        nxt = next((l for l in j.get("links", []) if l.get("rel") == "next"), None)
        if not nxt: return out
        body = nxt.get("body", body)

def nisar(lon, lat, names=("NISAR_L2_GCOV_BETA_V1", "NISAR_L2_GCOV_PROVISIONAL_V1", "NISAR_L2_GCOV_V1")):
    out = []
    for sn in names:
        q = urllib.parse.urlencode({"short_name": sn, "point": f"{lon},{lat}", "page_size": 1000, "temporal": f"{T0},{T1}"})
        for e in http("https://cmr.earthdata.nasa.gov/search/granules.json?" + q)["feed"]["entry"]:
            parts = e["title"].split("_")
            out.append(dict(time=e["time_start"], product=sn, track=parts[5] if len(parts) > 6 else "", direction=parts[6] if len(parts) > 6 else ""))
    return out

def one_per_pass(a):
    a = a.sort_values(["wban", "time"])
    new = (a.wban != a.wban.shift()) | (a.time.diff() > pd.Timedelta(hours=1))
    return a[new].reset_index(drop=True)

if __name__ == "__main__":
    st = pd.read_csv(os.path.join(PROC, "station_seasons.csv")); st = st[st.included].drop_duplicates("wban")
    CACHE = os.path.join(ROOT, "data", "raw", "looks_cache"); os.makedirs(CACHE, exist_ok=True)
    def fetch(r):
        f = os.path.join(CACHE, f"{r.wban}.json")
        if os.path.exists(f): return tuple(json.load(open(f)))
        res = ([dict(wban=r.wban, **d) for d in s1(r.lon, r.lat)], [dict(wban=r.wban, **d) for d in nisar(r.lon, r.lat)],
               [dict(wban=r.wban, **d) for d in nisar(r.lon, r.lat, names=("ALOS2_L1_PSR2",))])
        json.dump(res, open(f, "w")); return res
    with ThreadPoolExecutor(3) as ex: res = list(ex.map(fetch, st.itertuples()))
    out = {}
    for k, name in enumerate(("s1", "nisar", "alos2")):
        d = pd.DataFrame([x for r in res for x in r[k]], columns=["wban", "time"] if not any(r[k] for r in res) else None)
        d["time"] = pd.to_datetime(d.time, utc=True, format="ISO8601"); n0 = len(d)
        d = one_per_pass(d) if len(d) else d
        d.to_csv(os.path.join(PROC, f"looks_{name}.csv"), index=False); out[name] = d
        print(f"{name}: {n0} frames/granules -> {len(d)} looks at {d.wban.nunique() if len(d) else 0} of {len(st)} stations;",
              "by month", d.groupby(d.time.dt.strftime("%Y-%m")).size().to_dict() if len(d) else {})
