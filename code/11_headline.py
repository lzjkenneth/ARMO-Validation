"""Step 11 - collect every number quoted in README.md and docs/RESULTS.md into data/processed/headline.json, formatted the
way the documents quote them. tests/test_numbers.py checks that each formatted value appears in the documents, so the
write-up cannot drift away from the analysis outputs."""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from armo_common import PROC, load_events

H = {}
def put(key, value, fmt):
    H[key] = dict(value=float(value), text=fmt.format(value))

ev_all = pd.read_parquet(os.path.join(PROC, "events.parquet")); ev = load_events()
put("events_all", len(ev_all), "{:,.0f}"); put("stations_all", ev_all.wban.nunique(), "{:.0f}")
put("median_window_d", ev_all.W_half_h.median() / 24, "{:.2f}"); put("share_lt_48h", 100 * (ev_all.W_half_h < 48).mean(), "{:.0f}")
put("cohort_events", len(ev), "{:.0f}"); put("cohort_stations", ev.wban.nunique(), "{:.0f}")

s = pd.read_csv(os.path.join(PROC, "capture_summary.csv")).set_index("system")
for key, sys_ in (("s1", "Sentinel-1"), ("nisar", "NISAR"), ("fleet", "Sentinel-1 + NISAR"), ("armo", "ARMO D alone"),
                  ("both", "Fleet + ARMO D"), ("armoA", "ARMO A alone"), ("bothA", "Fleet + ARMO A"), ("bothB", "Fleet + ARMO B"),
                  ("bothC", "Fleet + ARMO C"), ("bothE", "Fleet + ARMO E")):
    r = s.loc[sys_]
    put(f"cap_{key}", 100 * r.capture, "{:.1f}"); put(f"ci_{key}_lo", 100 * r.ci_lo, "{:.1f}"); put(f"ci_{key}_hi", 100 * r.ci_hi, "{:.1f}")
    put(f"wait_med_{key}", r.wait_median_h, "{:.0f}"); put(f"wait_mean_{key}", r.wait_mean_h, "{:.0f}")
    put(f"wait_p90_{key}", r.wait_p90_h, "{:.0f}"); put(f"wait24_{key}", 100 * r.wait_le_24h, "{:.0f}")
    put(f"gap_{key}", r.mean_gap_d, "{:.2f}"); put(f"worst_{key}", r.worst_gap_d, "{:.1f}")
per = pd.read_parquet(os.path.join(PROC, "capture_events.parquet"))
from armo_common import boot_ci
gain = (per["cap_Fleet + ARMO D"] - per["cap_Sentinel-1 + NISAR"].astype(float)).to_numpy()
lo, hi = boot_ci(ev, gain)
put("gain", 100 * gain.mean(), "{:.1f}"); put("gain_lo", 100 * lo, "{:.1f}"); put("gain_hi", 100 * hi, "{:.1f}")
fl = per["cap_Sentinel-1 + NISAR"].astype(float)
put("armo_catches_fleet_missed", 100 * per["cap_ARMO D alone"][fl == 0].mean(), "{:.0f}")
W = per.W_half_h / 24
for lo_, hi_, k in ((0, 1, "lt1"), (1, 2, "1to2")):
    m = (W >= lo_) & (W < hi_)
    put(f"fleet_{k}", 100 * fl[m].mean(), "{:.0f}"); put(f"both_{k}", 100 * per["cap_Fleet + ARMO D"][m].mean(), "{:.0f}")
m = W < 0.5
put("share_lt12h", 100 * m.mean(), "{:.0f}"); put("both_lt12h", 100 * per["cap_Fleet + ARMO D"][m].mean(), "{:.0f}")

g = pd.read_csv(os.path.join(PROC, "gap_by_station.csv"))
put("armo_gap_mean", g.armo_gap.mean(), "{:.2f}"); put("armo_gap_south", g[g.lat < 35].armo_gap.mean(), "{:.2f}")
put("armo_gap_north", g[g.lat > 42].armo_gap.mean(), "{:.2f}"); put("armo_share_lt14", 100 * (g.armo_gap < 1.4).mean(), "{:.0f}")
put("s1_gap_median", g["Sentinel-1_gap"].median(), "{:.1f}"); put("nisar_gap_median", g["NISAR_gap"].median(), "{:.1f}")

bd = pd.read_csv(os.path.join(PROC, "capture_by_definition.csv")).pivot_table(index="system", columns="definition", values="capture")
for d in ("half", "half_p3", "efold", "abs"):
    put(f"gain_def_{d}", 100 * (bd.loc["Fleet + ARMO D", d] - bd.loc["Sentinel-1 + NISAR", d]), "{:.0f}")
    put(f"fleet_def_{d}", 100 * bd.loc["Sentinel-1 + NISAR", d], "{:.0f}"); put(f"both_def_{d}", 100 * bd.loc["Fleet + ARMO D", d], "{:.0f}")
ed = pd.read_csv(os.path.join(PROC, "event_definition_summary.csv")).set_index("rule")
put("win_half", ed.loc["50 % of the rise (report)", "median_d"], "{:.2f}"); put("win_p3", ed.loc["50 %, must stay below 3 h", "median_d"], "{:.2f}")
put("win_efold", ed.loc["1/e of the rise (37 %)", "median_d"], "{:.2f}"); put("win_abs", ed.loc["anomaly >= 0.03 m3/m3", "median_d"], "{:.2f}")
put("tau_pre", ed.loc["dry-down tau, floor = pre-storm level", "median_d"], "{:.1f}")
put("tau_ratio", ed.loc["peak-to-50 % time / (tau_pre ln2), uncensored", "median_d"], "{:.2f}")
put("censored_half", 100 * ed.loc["50 % of the rise (report)", "share_censored"], "{:.0f}")

sw = pd.read_csv(os.path.join(PROC, "design_sweeps.csv")).set_index("label")
for n in (1, 2, 3, 4):
    r = sw.loc[f"{n} satellite(s), one plane"]
    put(f"nsat{n}_alone", 100 * r.capture_alone, "{:.1f}"); put(f"nsat{n}_fleet", 100 * r.capture_with_fleet, "{:.1f}")
    put(f"nsat{n}_gap", r.mean_gap_d, "{:.2f}"); put(f"nsat{n}_wait", r.wait_median_alone_h, "{:.0f}")
for w_ in (200, 250, 300, 400, 420, 440, 450):
    r = sw.loc[f"2 sat, usable swath {w_} km"]
    put(f"sw{w_}_alone", 100 * r.capture_alone, "{:.1f}"); put(f"sw{w_}_fleet", 100 * r.capture_with_fleet, "{:.1f}")
    put(f"sw{w_}_share14", 100 * r.share_stations_lt_1p4, "{:.0f}"); put(f"sw{w_}_gap", r.mean_gap_d, "{:.2f}")
for lab, k in (("2 sats, dawn-dusk, asc 18:00 (proposed)", "tod_ddusk"), ("2 sats, asc 22:30 (looks ~10:30/22:30)", "tod_1030"),
               ("2 sats, asc 13:30 (looks ~01:30/13:30)", "tod_1330"), ("2 sats, look side alt", "side_alt"),
               ("2 sats, look side both", "side_both"), ("2 sats, 2nd sat in a 10:30/22:30 plane", "plane2"),
               ("1 sat, alternate sides, 25-45 deg", "one_alt"), ("1 sat, both sides every pass (bound)", "one_both"),
               ("1 sat, 20-50 deg (574 km access)", "one_wide"), ("1 sat at 800 km, 20-50 deg", "one_800")):
    r = sw.loc[lab]
    put(f"{k}_alone", 100 * r.capture_alone, "{:.1f}"); put(f"{k}_fleet", 100 * r.capture_with_fleet, "{:.1f}")
    put(f"{k}_swath", r.swath_km, "{:.0f}"); put(f"{k}_gap", r.mean_gap_d, "{:.2f}")
alt = pd.read_csv(os.path.join(PROC, "altitude_scan.csv")).set_index("h_km")
for h in (560, 620, 640, 720):
    put(f"alt{h}_alone", 100 * alt.loc[h, "capture_alone"], "{:.1f}"); put(f"alt{h}_worst", alt.loc[h, "worst_gap_d"], "{:.0f}")
rp = pd.read_csv(os.path.join(PROC, "repeat_orbit_options.csv")).set_index("orbit")
r6, r0 = rp.loc["89/6 (6-day repeat)"], rp.loc["620 km as proposed (no exact repeat)"]
put("rep6_h", r6.h_km, "{:.0f}"); put("rep6_alone", r6.capture_alone, "{:.1f}"); put("rep6_fleet", r6.capture_with_fleet, "{:.1f}")
put("rep6_worst", r6.worst_gap_d, "{:.1f}"); put("rep6_same", r6.same_geometry_repeat_median_d, "{:.0f}")
put("rep0_same", r0.same_geometry_repeat_median_d, "{:.0f}"); put("rep0_worst", r0.worst_gap_d, "{:.1f}")
lat = pd.read_csv(os.path.join(PROC, "nisar_latency_sample.csv"))
put("lat_gcov_med_h", lat[lat["product"] == "GCOV"].prod_lag_h.median(), "{:.0f}")
put("lat_sme2_med_h", lat[lat["product"] == "SME2"].prod_lag_h.median(), "{:.0f}")
put("lat_gcov_p90_d", lat[lat["product"] == "GCOV"].prod_lag_h.quantile(0.9) / 24, "{:.0f}")
put("lat_sme2_p90_d", lat[lat["product"] == "SME2"].prod_lag_h.quantile(0.9) / 24, "{:.0f}")
if os.path.exists(os.path.join(PROC, "noise_vs_scale.csv")):
    ns = pd.read_csv(os.path.join(PROC, "noise_vs_scale.csv")); ns = ns[ns.kind.isna() & (ns.lag_days == 12)].set_index(["pol", "scale_m"])
    for pol in ("HH", "HV"):
        for sc in (10, 100, 200, 500):
            put(f"noise_{pol}_{sc}", ns.loc[(pol, sc), "median_db"], "{:.2f}")
            put(f"speckle_{sc}", ns.loc[(pol, sc), "speckle_nisar_fit_db"], "{:.2f}")
            if sc >= 100:
                put(f"armo40_{sc}", ns.loc[(pol, sc), "speckle_armo_8b_40mhz_db"], "{:.2f}")
                put(f"armo20_{sc}", ns.loc[(pol, sc), "speckle_armo_11b_20mhz_db"], "{:.2f}")
    put("enl_fit", json.load(open(os.path.join(PROC, "noise_enl.json")))["enl_fit"], "{:.1f}")
    rs = pd.read_csv(os.path.join(PROC, "noise_station_residual.csv")).set_index(["pol", "scale_m"])
    for pol in ("HH", "HV"):
        for sc in (100, 200, 500):
            put(f"resid_{pol}_{sc}", rs.loc[(pol, sc), "resid_sd_db"], "{:.1f}")
    put("resid_HH_200_robust", rs.loc[("HH", 200), "resid_robust_db"], "{:.2f}")
    se = pd.read_csv(os.path.join(PROC, "noise_sensitivity.csv")).set_index(["pol", "scale_m"])
    put("sens_HH_200", se.loc[("HH", 200), "slope_db_per_0p1"], "{:.2f}")
    put("sens_HH_200_lo", se.loc[("HH", 200), "ci_lo"], "{:.2f}"); put("sens_HH_200_hi", se.loc[("HH", 200), "ci_hi"], "{:.2f}")
    put("noise_pairs", len(pd.read_csv(os.path.join(PROC, "noise_pairs.csv"))), "{:,.0f}")
    pr_ = pd.read_csv(os.path.join(PROC, "noise_pairs.csv")); put("noise_stations", pr_.wban.nunique(), "{:.0f}")
    ns_all = pd.read_csv(os.path.join(PROC, "noise_vs_scale.csv"))
    lag = ns_all[ns_all.kind.isna() & (ns_all.pol == "HH") & (ns_all.scale_m == 200)].set_index("lag_days")
    for L in (12, 24, 36):
        put(f"lag_HH200_{L}", lag.loc[L, "median_db"], "{:.2f}"); put(f"lag_n_{L}", lag.loc[L, "n_pairs"], "{:.0f}")
    put("sens_HV_200", se.loc[("HV", 200), "slope_db_per_0p1"], "{:.2f}")
# review checks (code/13_review_checks.py) and the Sentinel-1 interval test (code/12b_s1_noise_by_interval.py)
rc = pd.read_csv(os.path.join(PROC, "review_checks.csv"))
for k, r in rc.iterrows():
    H[f"check_{k}"] = dict(value=float("nan"), text=str(r["value"]), item=r["item"], note=str(r["note"]) if pd.notna(r["note"]) else "")
s1n = pd.read_csv(os.path.join(PROC, "s1_noise_by_lag.csv")).set_index(["pol", "scale_m", "lag_days"])
for L in (6, 12, 18, 24):
    put(f"s1_spatial_VV200_{L}", s1n.loc[("VV", 200, L), "spatial_median_db"], "{:.2f}")
    put(f"s1_resid_VV200_{L}", s1n.loc[("VV", 200, L), "station_resid_sd_db"], "{:.2f}")
    put(f"s1_n_{L}", s1n.loc[("VV", 200, L), "n_stable"], "{:.0f}")
m6 = pd.read_csv(os.path.join(PROC, "s1_noise_6v12_matched.csv")).set_index(["pol", "scale_m"])
put("s1_ratio_VV200", m6.loc[("VV", 200), "median_ratio_12_over_6"], "{:.2f}")
put("s1_ratio_VV200_lo", m6.loc[("VV", 200), "ci_lo"], "{:.2f}"); put("s1_ratio_VV200_hi", m6.loc[("VV", 200), "ci_hi"], "{:.2f}")
put("s1_ratio_stations", m6.loc[("VV", 200), "stations"], "{:.0f}")
json.dump(H, open(os.path.join(PROC, "headline.json"), "w"), indent=1)
for k, v in H.items():
    print(f"{k:28} {v['text']}")
