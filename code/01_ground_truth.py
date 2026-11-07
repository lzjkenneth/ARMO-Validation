"""Step 1 - ground truth. Download NOAA USCRN hourly02 station files for 2025 and 2026 and build a quality-controlled
hourly table of 5 cm soil moisture and rain. Written independently of the team's code, following the definitions in
their report (Section 3.1): missing sentinels to NaN, 5 cm soil moisture outside [0, 0.6] to NaN, frozen or unknown soil
(5 cm soil temperature not above 1 C) to NaN, warm season 1 Apr-31 Oct, a station-season counts if >= 70 % of its hours
have valid soil moisture.
Source: https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/ (open, no login)."""
import os, re, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw", "uscrn"); PROC = os.path.join(ROOT, "data", "processed")
BASE = "https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/"
H = {"User-Agent": "ES4304-student-validation"}
# hourly02 column layout (HEADERS.txt): 38 fields
COLS = ["WBANNO", "UTC_DATE", "UTC_TIME", "LST_DATE", "LST_TIME", "CRX_VN", "LONGITUDE", "LATITUDE", "T_CALC", "T_HR_AVG",
        "T_MAX", "T_MIN", "P_CALC", "SOLARAD", "SOLARAD_FLAG", "SOLARAD_MAX", "SOLARAD_MAX_FLAG", "SOLARAD_MIN",
        "SOLARAD_MIN_FLAG", "SUR_TEMP_TYPE", "SUR_TEMP", "SUR_TEMP_FLAG", "SUR_TEMP_MAX", "SUR_TEMP_MAX_FLAG", "SUR_TEMP_MIN",
        "SUR_TEMP_MIN_FLAG", "RH_HR_AVG", "RH_HR_AVG_FLAG", "SOIL_MOISTURE_5", "SOIL_MOISTURE_10", "SOIL_MOISTURE_20",
        "SOIL_MOISTURE_50", "SOIL_MOISTURE_100", "SOIL_TEMP_5", "SOIL_TEMP_10", "SOIL_TEMP_20", "SOIL_TEMP_50", "SOIL_TEMP_100"]

def get(url, tries=4):
    for k in range(tries):
        try: return urllib.request.urlopen(urllib.request.Request(url, headers=H), timeout=120).read()
        except Exception:
            if k == tries - 1: raise
            time.sleep(3 * (k + 1))

def download():
    for yr in (2025, 2026):
        os.makedirs(os.path.join(RAW, str(yr)), exist_ok=True)
        names = sorted(set(re.findall(r'href="(CRNH0203-%d-[^"]+\.txt)"' % yr, get(BASE + f"{yr}/").decode())))
        todo = [n for n in names if not os.path.exists(os.path.join(RAW, str(yr), n))]
        print(yr, len(names), "files,", len(todo), "to fetch", flush=True)
        def one(n):
            open(os.path.join(RAW, str(yr), n), "wb").write(get(BASE + f"{yr}/{n}"))
        with ThreadPoolExecutor(6) as ex: list(ex.map(one, todo))

def parse_file(path):
    rows = []
    with open(path, "rb") as f:
        for ln in f.read().replace(b"\x00", b"").decode("ascii", "ignore").splitlines():
            p = ln.split()
            if len(p) == 38 and p[0].isdigit(): rows.append(p)
    if not rows: return None
    d = pd.DataFrame(rows, columns=COLS)
    out = pd.DataFrame({"wban": d.WBANNO.astype(int),
                        "time": pd.to_datetime(d.UTC_DATE + d.UTC_TIME.str.zfill(4), format="%Y%m%d%H%M", utc=True),
                        "lon": pd.to_numeric(d.LONGITUDE), "lat": pd.to_numeric(d.LATITUDE),
                        "rain": pd.to_numeric(d.P_CALC), "sm5": pd.to_numeric(d.SOIL_MOISTURE_5),
                        "st5": pd.to_numeric(d.SOIL_TEMP_5)})
    for c in ("rain", "sm5", "st5"): out.loc[out[c] <= -99, c] = np.nan
    out.loc[(out.sm5 < 0) | (out.sm5 > 0.6), "sm5"] = np.nan
    out.loc[~(out.st5 > 1.0), "sm5"] = np.nan          # frozen or unknown soil temperature
    out["name"] = re.sub(r"CRNH0203-20\d\d-|\.txt", "", os.path.basename(path))
    return out

def build():
    files = sorted([os.path.join(RAW, y, f) for y in ("2025", "2026") for f in os.listdir(os.path.join(RAW, y))])
    with ThreadPoolExecutor(8) as ex: parts = [p for p in ex.map(parse_file, files) if p is not None]
    df = pd.concat(parts, ignore_index=True).sort_values(["wban", "time"]).drop_duplicates(["wban", "time"])
    df = df[(df.lat > 24) & (df.lat < 50) & (df.lon > -125) & (df.lon < -66)]      # contiguous US
    os.makedirs(PROC, exist_ok=True)
    df.to_parquet(os.path.join(PROC, "uscrn_hourly.parquet"))
    rows = []
    for yr in (2025, 2026):
        w = df[(df.time >= f"{yr}-04-01") & (df.time < f"{yr}-11-01")]
        g = w.groupby("wban").agg(name=("name", "first"), lat=("lat", "median"), lon=("lon", "median"), sm_valid=("sm5", "count"))
        g["season"] = yr; g["valid_frac"] = g.sm_valid / (214 * 24); rows.append(g.reset_index())
    st = pd.concat(rows, ignore_index=True); st["included"] = st.valid_frac >= 0.7
    st.to_csv(os.path.join(PROC, "station_seasons.csv"), index=False)
    print("hourly rows", len(df), "| stations", df.wban.nunique(), "| station-seasons included", int(st.included.sum()),
          "| stations included in >= 1 season", st[st.included].wban.nunique())

if __name__ == "__main__":
    if "--build-only" not in sys.argv: download()
    build()
