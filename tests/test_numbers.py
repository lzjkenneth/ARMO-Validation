"""Checks that the documents quote the analysis correctly, and that the orbit model and the event rule behave as
described. Run: python -m pytest -q tests"""
import json, math, os, sys, importlib.util
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H = json.load(open(os.path.join(ROOT, "data", "processed", "headline.json")))
DOCS = {n: open(os.path.join(ROOT, p), encoding="utf-8").read() for n, p in
        (("results", "docs/RESULTS.md"), ("readme", "README.md"))}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, path)); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m


# keys whose formatted value must appear in RESULTS.md
RESULTS_KEYS = """events_all stations_all cohort_events cohort_stations median_window_d share_lt_48h
cap_s1 cap_nisar cap_fleet cap_armo cap_both ci_s1_lo ci_s1_hi ci_nisar_lo ci_nisar_hi ci_fleet_lo ci_fleet_hi ci_armo_lo
ci_armo_hi ci_both_lo ci_both_hi gain gain_lo gain_hi cap_bothA cap_bothB cap_bothC cap_bothE
wait_med_fleet wait_med_armo wait_med_both wait_mean_fleet wait_mean_armo wait_mean_both wait_p90_fleet wait_p90_armo
wait_p90_both wait24_fleet wait24_armo wait24_both gap_fleet gap_armo gap_both worst_fleet worst_armo
armo_gap_mean armo_gap_south armo_gap_north armo_share_lt14 armo_catches_fleet_missed fleet_lt1 both_lt1 fleet_1to2
both_1to2 share_lt12h both_lt12h win_half win_p3 win_efold win_abs tau_pre tau_ratio censored_half
gain_def_half gain_def_half_p3 gain_def_efold gain_def_abs fleet_def_half fleet_def_efold both_def_efold
nsat1_alone nsat1_fleet nsat1_gap nsat1_wait nsat2_alone nsat2_fleet nsat2_gap nsat2_wait nsat3_alone nsat3_fleet
nsat3_gap nsat3_wait nsat4_alone nsat4_fleet nsat4_gap nsat4_wait sw200_alone sw200_fleet sw250_alone sw250_fleet
sw300_alone sw300_fleet sw400_alone sw400_fleet sw400_share14 sw420_alone sw420_share14 sw440_alone sw440_fleet
tod_ddusk_alone tod_1030_alone tod_1330_alone side_alt_alone side_both_alone plane2_alone one_alt_alone one_wide_alone
one_800_alone one_800_swath one_both_alone alt560_alone alt620_alone alt640_alone alt720_alone rep6_alone rep6_fleet
rep6_worst rep6_same rep0_same noise_HH_10 noise_HH_100 noise_HH_200 noise_HH_500 noise_HV_10 noise_HV_100
noise_HV_200 noise_HV_500 speckle_100 speckle_200 armo40_100 armo40_200 armo20_200 enl_fit sens_HH_200 sens_HH_200_lo
sens_HH_200_hi sens_HV_200 resid_HH_100 resid_HH_200 resid_HH_500 resid_HH_200_robust lag_HH200_12 lag_HH200_24
lag_HH200_36 s1_spatial_VV200_6 s1_spatial_VV200_24 s1_resid_VV200_6 s1_resid_VV200_12 s1_resid_VV200_24 s1_n_6 s1_n_12
s1_ratio_VV200 s1_ratio_VV200_lo s1_ratio_VV200_hi""".split()


def test_results_quote_every_headline_number():
    missing = [f"{k}={H[k]['text']}" for k in RESULTS_KEYS if H[k]["text"] not in DOCS["results"]]
    assert not missing, "RESULTS.md does not quote: " + ", ".join(missing)


def test_readme_table_matches():
    for k in ("cap_fleet", "cap_armo", "cap_both", "ci_fleet_lo", "ci_both_hi", "wait_med_fleet", "wait_med_both",
              "wait24_fleet", "wait24_both", "gap_fleet", "gap_armo", "gap_both", "worst_fleet", "gain", "alt560_alone",
              "alt720_alone", "armo_gap_south", "armo_gap_north"):
        assert H[k]["text"] in DOCS["readme"], k


def test_sensitivity_checks_present_in_results():
    for k in ("check_0", "check_1", "check_2"):
        assert H[k]["text"].split(" (")[0] in DOCS["results"], k


def test_orbit_geometry():
    orb = load("orb", "code/04_orbit_sim.py")
    assert abs(orb.sso_inclination(620) - 97.87) < 0.02
    assert abs(orb.access_width_km(620, 25, 45) - 362) < 1.5
    st = pd.DataFrame({"lat": [36.12], "lon": [-97.09]})
    L = orb.looks(st, "2026-07-01", "2026-07-31", phase_deg=0)[0]
    assert 8 <= len(L) <= 14                       # one right-looking satellite: a look every ~2.8 days
    G = orb.looks(st, "2026-07-01", "2026-07-31", phase_deg=0, return_geom=True)[0]
    assert np.all((G[:, 1] >= 25) & (G[:, 1] <= 45))


def test_event_rule_on_a_synthetic_storm():
    ev = load("ev", "code/02_events.py")
    idx = pd.date_range("2026-06-01", "2026-06-20", freq="h", tz="UTC")
    sm = np.full(len(idx), 0.20); rain = np.zeros(len(idx))
    k0 = 100; rain[k0:k0 + 3] = 4.0                 # 12 mm storm
    t = np.arange(len(idx)) - (k0 + 3)
    sm[k0 + 3:] = 0.20 + 0.10 * np.exp(-t[k0 + 3:] / 48.0)   # rise of 0.10, e-folding 48 h
    s = pd.DataFrame({"sm5": sm, "rain": rain}, index=idx)
    out = ev.detect(s, 2026, 1)
    assert len(out) == 1
    e = out[0]
    assert abs(e["amp"] - 0.10) < 1e-6
    half = (e["we_half"] - e["t_peak"]) / pd.Timedelta(hours=1)
    assert abs(half - 48 * math.log(2)) <= 1.0     # 50 % of the rise is reached after tau ln 2
    assert abs(e["tau_pre_h"] - 48) < 1.0
