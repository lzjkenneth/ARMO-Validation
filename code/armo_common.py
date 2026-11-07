"""Shared helpers for the capture, waiting-time and design-sweep scripts.

Definitions used everywhere
  look      : one satellite pass that images the station (frames < 1 h apart merged)
  captured  : at least one look inside the event's wetting window [ws, we]
  wait      : hours from the window opening (ws) to the first look at or after it ("waiting time"); if that first look comes
              after the window has closed, the event is missed
  mean gap  : for one station, the number of days in the period divided by the number of looks (equals the average gap
              between consecutive looks); "mean of gaps" pools every consecutive gap instead
Cohort: storms starting 17 Jun-31 Oct 2026 (the continuous NISAR record), as in the team's report."""
import os, importlib.util
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROC = os.path.join(ROOT, "data", "processed")
FIG = os.path.join(ROOT, "figures")
_spec = importlib.util.spec_from_file_location("orb", os.path.join(ROOT, "code", "04_orbit_sim.py"))
orb = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(orb)
ERA0, ERA1 = pd.Timestamp("2026-06-17", tz="UTC"), pd.Timestamp("2026-11-01", tz="UTC")
SIM0, SIM1 = "2026-06-10", "2026-11-25"           # simulate a margin either side so late windows can still be looked at
PHASES = [0, 45, 90, 135, 180, 225, 270, 315]
EPOCH = pd.Timestamp("1970-01-01", tz="UTC")
DEFS = {"half": "50 % of rise (report)", "half_p3": "50 %, 3 h below to end", "efold": "1/e of rise", "abs": "anomaly >= 0.03",
        "q25": "25 % of rise", "q75": "75 % of rise"}


def hrs(s):
    return ((pd.to_datetime(s, utc=True) - EPOCH) / pd.Timedelta(hours=1)).to_numpy(dtype=float)


def load_events(cohort="nisar_era"):
    ev = pd.read_parquet(os.path.join(PROC, "events.parquet"))
    if cohort == "nisar_era":
        ev = ev[(ev.season == 2026) & (ev.t_start >= ERA0) & (ev.t_start < ERA1)]
    ev = ev.reset_index(drop=True)
    for d in DEFS:
        ev[f"ws_{d}_h"], ev[f"we_{d}_h"] = hrs(ev[f"ws_{d}"]), hrs(ev[f"we_{d}"])
    return ev


def stations_for(ev):
    st = pd.read_csv(os.path.join(PROC, "station_seasons.csv")).drop_duplicates("wban").set_index("wban")
    return st.loc[sorted(ev.wban.unique()), ["name", "lat", "lon"]].reset_index()


def fleet_looks():
    by = lambda d: {w: np.sort(hrs(g.time)) for w, g in d.groupby("wban")}
    S1 = by(pd.read_csv(os.path.join(PROC, "looks_s1.csv")))
    NI = by(pd.read_csv(os.path.join(PROC, "looks_nisar.csv")))
    FL = {w: np.sort(np.r_[S1.get(w, []), NI.get(w, [])]) for w in set(S1) | set(NI)}
    return {"Sentinel-1": S1, "NISAR": NI, "Sentinel-1 + NISAR": FL}


_cache = {}
def sat_looks(stations, phase, **kw):
    """Look times for one simulated satellite at one starting phase (cached)."""
    key = (tuple(stations.wban), phase, tuple(sorted(kw.items())))
    if key not in _cache:
        L = orb.looks(stations, SIM0, SIM1, phase_deg=phase, **kw)
        _cache[key] = {w: np.sort(L[i]) for i, w in enumerate(stations.wban)}
    return _cache[key]


def constellation(stations, n_sats=2, phase=0.0, planes=None, **kw):
    """n_sats equally spaced in one plane (or a list of (phase_offset, kw_override) in `planes`)."""
    parts = planes if planes is not None else [((360.0 / n_sats) * k, {}) for k in range(n_sats)]
    looks = {}
    for off, over in parts:
        L = sat_looks(stations, (phase + off) % 360, **{**kw, **over})
        for w, o in L.items():
            looks[w] = np.r_[looks.get(w, []), o]
    return {w: np.sort(o) for w, o in looks.items()}


def merge(*sets):
    keys = set().union(*[s.keys() for s in sets])
    return {w: np.sort(np.concatenate([s.get(w, np.array([])) for s in sets])) for w in keys}


def score(ev, looks_by, d="half"):
    """captured (bool) and wait (h) per event for window definition d."""
    cap = np.zeros(len(ev), bool); wait = np.full(len(ev), np.nan)
    WS, WE = ev[f"ws_{d}_h"].to_numpy(), ev[f"we_{d}_h"].to_numpy()
    for w, idx in ev.groupby("wban").indices.items():
        obs = looks_by.get(w, np.array([]))
        if obs.size == 0:
            continue
        ws, we = WS[idx], WE[idx]
        i = np.searchsorted(obs, ws, side="left")
        f = np.where(i < len(obs), obs[np.minimum(i, len(obs) - 1)], np.nan)
        cap[idx] = np.isfinite(f) & (f <= we); wait[idx] = f - ws
    return cap, wait


def gap_table(looks_by, wbans, t0=ERA0, t1=ERA1):
    a, b = (t0 - EPOCH) / pd.Timedelta(hours=1), (t1 - EPOCH) / pd.Timedelta(hours=1)
    rows = []
    for w in wbans:
        o = looks_by.get(w, np.array([])); o = o[(o >= a) & (o < b)]
        g = np.diff(o) / 24.0
        rows.append(dict(wban=w, n_looks=len(o), mean_gap_d=(b - a) / 24.0 / max(len(o), 1),
                         max_gap_d=g.max() if g.size else (b - a) / 24.0, gaps=g))
    return pd.DataFrame(rows)


def gap_summary(gt):
    allg = np.concatenate([g for g in gt.gaps if len(g)]) if len(gt) else np.array([])
    return dict(mean_gap_d=gt.mean_gap_d.mean(), median_station_gap_d=gt.mean_gap_d.median(),
                p90_station_gap_d=gt.mean_gap_d.quantile(0.9), mean_of_gaps_d=allg.mean() if allg.size else np.nan,
                share_stations_lt_1p4=(gt.mean_gap_d < 1.4).mean(), mean_max_gap_d=gt.max_gap_d.mean(),
                worst_gap_d=gt.max_gap_d.max())


def boot_ci(ev, values, n=2000, seed=42):
    """Station-cluster bootstrap 95 % CI of the event mean of `values`."""
    rng = np.random.default_rng(seed)
    st = ev.wban.to_numpy(); u, inv = np.unique(st, return_inverse=True)
    s = np.bincount(inv, weights=values, minlength=len(u)); c = np.bincount(inv, minlength=len(u)).astype(float)
    k = rng.integers(0, len(u), (n, len(u)))
    return np.percentile(s[k].sum(1) / c[k].sum(1), [2.5, 97.5])


def phase_average(ev, make_looks, d="half", base=None, phases=PHASES, per_phase=False):
    """Capture averaged over starting phases (per event) and the waits. make_looks(phase) -> looks dict; base: looks merged
    in (the real fleet). With per_phase=True the waits come back as a (phases x events) array so their distribution can be
    pooled over phases; averaging waits per event first would narrow the distribution and flatter ARMO."""
    caps, waits = [], []
    for ph in phases:
        L = make_looks(ph)
        if base is not None:
            L = merge(L, base)
        c, w = score(ev, L, d); caps.append(c); waits.append(w)
    if per_phase:
        return np.mean(caps, axis=0), np.array(waits)
    return np.mean(caps, axis=0), np.nanmean(waits, axis=0)


def wait_stats(w):
    """Summary of waiting times (hours); w may be 1-D (one real record) or phases x events (pooled)."""
    w = np.asarray(w, float).ravel()
    return dict(wait_mean_h=np.nanmean(w), wait_median_h=np.nanmedian(w), wait_p90_h=np.nanpercentile(w, 90),
                wait_le_24h=np.nanmean(w <= 24), wait_le_48h=np.nanmean(w <= 48))
