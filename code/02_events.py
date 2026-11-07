"""Step 2 - wetting events from hourly 5 cm soil moisture (independent re-implementation of the report's definition,
Section 3.2), plus two alternative definitions to answer "how do you decide what a wetting event is?":

  main (half-amplitude): storm = hours with >= 0.2 mm rain, merged when < 24 h apart; event needs >= 5 mm and a 5 cm rise
       >= 0.03 m3/m3 above the 3-hour pre-storm mean; the window runs while the anomaly stays >= 50 % of the rise, and ends
       at the next storm (censored) or after 30 days.
  abs  : the window runs while the anomaly stays >= 0.03 m3/m3 (a fixed physical threshold).
  efold: the window runs while the anomaly stays >= 1/e of the rise (37 %); with an exponential dry-down this is exactly
       the dry-down time constant tau, the quantity soil-moisture papers report.
  half_p3: the report's 50 % rule, but the window only ends after 3 hours in a row below the line, so one noisy hourly
       reading cannot end it (robustness variant).
  q25, q75: the window runs while the anomaly stays >= 25 % or >= 75 % of the rise (sweep of the 50 % choice itself).
For every event two dry-down curves are fitted from the peak (up to 7 days or the next storm):
  floor fixed at the pre-storm level: SM = pre + A exp(-(t - t_peak) / tau_pre)   (matches the 50 % rule: halves after
       tau_pre * ln 2)
  free floor (literature style):      SM = c + a exp(-(t - t_peak) / tau_free)
Output: data/processed/events.parquet"""
import os, warnings
import numpy as np, pandas as pd
from scipy.optimize import curve_fit

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROC = os.path.join(ROOT, "data", "processed")
WET_MM, STORM_GAP_H, P_MIN, A_MIN, FRAC = 0.2, 24, 5.0, 0.03, 0.5
warnings.simplefilter("ignore")

def series(g, yr):
    idx = pd.date_range(f"{yr}-03-25", f"{yr}-11-30", freq="h", tz="UTC")
    s = g.set_index("time")[["rain", "sm5"]].reindex(idx)
    s["sm5"] = s.sm5.interpolate(limit=3, limit_area="inside")       # gaps up to 3 h, as in the report
    s["rain"] = s.rain.fillna(0.0)
    return s

def drydown(t_h, y, pre, A):
    """Two exponential dry-down fits from the peak. Returns (tau_pre_h, c, a, tau_free_h); nan where a fit fails."""
    ok = np.isfinite(y)
    if ok.sum() < 12: return (np.nan,) * 4
    t, y = t_h[ok], y[ok]
    try:
        p1, _ = curve_fit(lambda t, tau: pre + A * np.exp(-t / tau), t, y, p0=(48.0,), bounds=([1.0], [5000.0]), maxfev=4000)
        tau_pre = p1[0]
    except Exception:
        tau_pre = np.nan
    try:
        p2, _ = curve_fit(lambda t, c, a, tau: c + a * np.exp(-t / tau), t, y, p0=(y[-1], y[0] - y[-1], 24.0),
                          bounds=([-0.1, 0.0, 1.0], [0.7, 0.7, 2000.0]), maxfev=4000)
        free = tuple(p2) if p2[1] > 0.005 else (np.nan,) * 3
    except Exception:
        free = (np.nan,) * 3
    return (tau_pre,) + free

def window(sm, a, ip, lim, thr, persist=1):
    """Window [ws, we): opens at the first hour after onset at or above thr, closes at the first hour after the peak that
    starts a run of `persist` hours below thr (missing hours never close it). Censored if it never closes before lim."""
    seg = np.nan_to_num(sm[a:ip + 1], nan=-1.0)
    up = np.flatnonzero(seg >= thr)
    if up.size == 0: return None
    ws = a + int(up[0])
    below = np.nan_to_num(sm[ip:lim], nan=np.inf) < thr
    if persist > 1 and below.size >= persist:
        run = np.convolve(below.astype(int), np.ones(persist, int), "valid") == persist
        dn = np.flatnonzero(run)
    else:
        dn = np.flatnonzero(below)
    return ws, ip + (int(dn[0]) if dn.size else len(below)), dn.size == 0

def detect(s, yr, wban):
    sm, p, t = s.sm5.to_numpy(), s.rain.to_numpy(), s.index
    wet = np.flatnonzero(p >= WET_MM)
    if wet.size == 0: return []
    br = np.flatnonzero(np.diff(wet) > STORM_GAP_H)
    starts, ends = np.r_[wet[0], wet[br + 1]], np.r_[wet[br], wet[-1]]
    tot = np.array([p[a:b + 1].sum() for a, b in zip(starts, ends)])
    keep = tot >= P_MIN; starts, ends, tot = starts[keep], ends[keep], tot[keep]
    lo, hi = t.searchsorted(pd.Timestamp(f"{yr}-04-01", tz="UTC")), t.searchsorted(pd.Timestamp(f"{yr}-11-01", tz="UTC"))
    out = []
    for k, (a, b, P) in enumerate(zip(starts, ends, tot)):
        if a < lo or a >= hi or a < 3: continue
        pre = np.nanmean(sm[a - 3:a])
        if not np.isfinite(pre): continue
        nxt = starts[k + 1] if k + 1 < len(starts) else len(sm)
        seg = sm[a:min(b + 24, nxt) + 1]
        if not np.isfinite(seg).any(): continue
        ip = a + int(np.nanargmax(seg)); A = sm[ip] - pre
        if A < A_MIN: continue
        lim = min(nxt, ip + 30 * 24, len(sm))
        w_half = window(sm, a, ip, lim, pre + FRAC * A)
        w_abs = window(sm, a, ip, lim, pre + A_MIN)
        w_e = window(sm, a, ip, lim, pre + A / np.e)
        w_p3 = window(sm, a, ip, lim, pre + FRAC * A, persist=3)
        w_q25 = window(sm, a, ip, lim, pre + 0.25 * A); w_q75 = window(sm, a, ip, lim, pre + 0.75 * A)
        if w_half is None: continue
        fit_end = min(nxt, ip + 7 * 24, len(sm))
        tau_pre, c_free, a_free, tau_free = drydown(np.arange(fit_end - ip, dtype=float), sm[ip:fit_end], pre, A)
        rec = dict(wban=wban, season=yr, t_start=t[a], t_peak=t[ip], rain_mm=P, pre=pre, amp=A, rise_h=ip - a,
                   tau_pre_h=tau_pre, tau_free_h=tau_free, c_free=c_free, a_free=a_free,
                   next_storm=t[nxt] if nxt < len(t) else pd.NaT)
        for name, w in (("half", w_half), ("abs", w_abs), ("efold", w_e), ("half_p3", w_p3), ("q25", w_q25), ("q75", w_q75)):
            if w is None: rec.update({f"ws_{name}": pd.NaT, f"we_{name}": pd.NaT, f"cens_{name}": np.nan}); continue
            rec.update({f"ws_{name}": t[w[0]], f"we_{name}": t[min(w[1], len(t) - 1)], f"cens_{name}": w[2]})
        out.append(rec)
    return out

if __name__ == "__main__":
    df = pd.read_parquet(os.path.join(PROC, "uscrn_hourly.parquet"))
    st = pd.read_csv(os.path.join(PROC, "station_seasons.csv"))
    rows = []
    for r in st[st.included].itertuples():
        g = df[df.wban == r.wban]
        rows += detect(series(g, r.season), r.season, r.wban)
    ev = pd.DataFrame(rows)
    for n in ("half", "abs", "efold", "half_p3", "q25", "q75"): ev[f"W_{n}_h"] = (ev[f"we_{n}"] - ev[f"ws_{n}"]) / pd.Timedelta(hours=1)
    ev.to_parquet(os.path.join(PROC, "events.parquet"))
    print("events", len(ev), "| stations", ev.wban.nunique(), "| median half window", round(ev.W_half_h.median() / 24, 2), "d",
          "| share < 48 h", round((ev.W_half_h < 48).mean(), 3), "| median 3-h-persistent window",
          round(ev.W_half_p3_h.median() / 24, 2), "d | median tau (pre floor)", round(ev.tau_pre_h.median(), 1), "h | median tau (free floor)",
          round(ev.tau_free_h.median(), 1), "h")
