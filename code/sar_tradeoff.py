#!/usr/bin/env python3
"""ARMO SAR trade-off: swath, resolution, looks, speckle noise, NESZ, power, data and revisit.

Companion to docs/RESOLUTION_AND_FEASIBILITY.md. Pure Python 3 + numpy + matplotlib.

Run:  python code/sar_tradeoff.py
Writes: figures/feas_*.png and data/processed/feas_*.csv (paths relative to the repo root).

Every number here is either (a) a CALCULATION from the formulas below, (b) an ASSUMPTION listed
in ASSUME with its basis, or (c) a published mission value in MISSIONS with its source URL.
Nothing is fitted to the team's results; the team's own numbers (362 km corridor, 1.4 d gap,
2.8 kW peak, 316 kg, 0.38 dB target, 1.45 dB measured) are quoted only for comparison.
"""
import csv
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
DAT = os.path.join(ROOT, "data", "processed")
os.makedirs(FIG, exist_ok=True)
os.makedirs(DAT, exist_ok=True)

ACCESS = "2026-11-07"          # date on which every URL below was read

# ----------------------------------------------------------------------------------------------
# constants
# ----------------------------------------------------------------------------------------------
RE = 6371.0                    # km, mean Earth radius (spherical Earth throughout)
MU = 398600.4418               # km^3/s^2
C = 299792458.0                # m/s
KB = 1.380649e-23              # J/K
T0 = 290.0                     # K, reference noise temperature
DB = 10.0 / math.log(10.0)     # 4.343, converts a natural-log std to dB

# ----------------------------------------------------------------------------------------------
# ARMO as proposed, plus the engineering assumptions used here (value, basis)
# ----------------------------------------------------------------------------------------------
ASSUME = dict(
    h_km=(620.0, "team's orbit"),
    f_hz=(1.2575e9, "team's centre frequency; NISAR and ALOS-2 use 1257.5 MHz too"),
    D_m=(5.9, "team's reflector diameter"),
    eta_min=(25.0, "team's near edge, off-nadir deg"),
    eta_max=(45.0, "team's far edge, off-nadir deg"),
    eta_ap=(0.60, "aperture efficiency of an offset mesh reflector (team's value; typical 0.55-0.70)"),
    k_bw=(1.15, "3 dB beamwidth = k_bw * lambda / D for a tapered circular aperture (1.03 uniform .. 1.2 tapered)"),
    k_prf=(1.15, "PRF = k_prf * Doppler bandwidth: 15 % oversampling for azimuth ambiguities (ASSUMPTION; "
                 "reproduces NISAR's 1650 Hz from its 12 m reflector)"),
    k_win=(1.20, "resolution broadening by the processing window (Hamming-like) (ASSUMPTION)"),
    k_burst=(1.70, "ScanSAR/TOPS azimuth resolution = k_burst * n_beams * stripmap resolution; the ideal lower bound "
                   "is (n_beams + 1), practical systems reach ~2 n_beams; 1.7 reproduces Sentinel-1 IW (3 bursts, 4.3-4.9 m -> 22 m) "
                   "and EW (5 bursts, -> 43 m) from the Sentinel-1 product table (sentiwiki)"),
    guard_us=(5.0, "receive-protection guard each side of the transmit pulse, microseconds (ASSUMPTION)"),
    duty=(0.10, "transmit duty cycle (team's value; solid-state SARs run 5-15 %)"),
    nf_loss_db=(6.0, "receiver noise figure + RF/processing losses (team's value)"),
    edge_db=(3.0, "two-way beam-edge loss with elevation beams overlapped at -1.5 dB (team's value)"),
    scan_margin_db=(3.0, "extra ScanSAR/scalloping margin (team's value)"),
    bw_hz=(20e6, "baseline chirp bandwidth (team); NISAR uses 5/20/40/80 MHz inside 1215-1300 MHz"),
    p_peak_w=(2800.0, "team's revised peak RF power"),
    nesz_req_db=(-25.0, "team's NESZ requirement at the far edge"),
    bits=(4, "BAQ bits per sample (team; SAOCOM also uses 4-bit BAQ)"),
    n_pol=(2, "HH + HV"),
    rx_duty=(0.70, "fraction of the pulse interval spent receiving (team)"),
    oversamp=(1.20, "complex sampling rate / chirp bandwidth (team)"),
    xband_bps=(400e6, "downlink rate of a commercial X-band radio (team)"),
    duty_conus=(0.015, "fraction of time imaging CONUS, from the team's orbit model (Table 3: 1.5 %)"),
    duty_global=(0.30, "what a global-land mapper would need (ASSUMPTION; ALOS-2 is designed for 50 %, Sentinel-1 for ~25-30 min per 99-min orbit)"),
    sspa_eff=(0.35, "DC-to-RF efficiency end to end (team)"),
    radar_elec_w=(150.0, "radar electronics while imaging (team)"),
    bus_w=(150.0, "bus orbit-average power (team)"),
    lat_ref=(40.0, "reference latitude for the revisit formula; USCRN stations lie mostly at 30-48 N"),
    dry_mass_kg=(290.0, "dry mass used only for the propellant estimate (team's 316 kg wet minus ~26 kg propellant)"),
    isp_s=(220.0, "monopropellant hydrazine (team)"),
    hh_db=(-15.0, "typical L-band HH backscatter of moist cropland at 30-45 deg, dB (ASSUMPTION, order of magnitude)"),
    hv_db=(-24.0, "typical L-band HV backscatter of cropland, dB (ASSUMPTION, order of magnitude)"),
    tec_dawn=(10.0, "ionospheric TEC at 06:00 local, TECU (ASSUMPTION: typical quiet mid-latitude value)"),
    tec_pm=(50.0, "ionospheric TEC at 14:00-16:00 local near solar maximum, TECU (ASSUMPTION)"),
    b_field_t=(4.5e-5, "geomagnetic field, tesla (ASSUMPTION, mid-latitude)"),
    cos_b=(0.5, "cosine of the angle between the ray and the field line (ASSUMPTION)"),
)
A = {k: v[0] for k, v in ASSUME.items()}


# ----------------------------------------------------------------------------------------------
# geometry and orbit
# ----------------------------------------------------------------------------------------------
def geom(eta_deg, h_km):
    """Off-nadir angle -> slant range (km), incidence angle (deg), ground range from nadir (km)."""
    r = RE + h_km
    e = np.deg2rad(np.asarray(eta_deg, float))
    R = r * np.cos(e) - np.sqrt(RE ** 2 - (r * np.sin(e)) ** 2)
    inc = np.arcsin(r * np.sin(e) / RE)
    g = RE * (inc - e)
    return R, np.rad2deg(inc), g


def eta_from_inc(inc_deg, h_km):
    """Incidence angle -> off-nadir angle (deg)."""
    return np.rad2deg(np.arcsin(RE / (RE + h_km) * np.sin(np.deg2rad(inc_deg))))


def dg_deta(eta_deg, h_km):
    """Ground range per radian of off-nadir angle (km/rad): how wide one elevation beam is on the ground."""
    r = RE + h_km
    e = np.deg2rad(eta_deg)
    inc = np.arcsin(r * np.sin(e) / RE)
    return r * np.cos(e) / np.cos(inc) - RE


def orbit(h_km):
    r = RE + h_km
    v = math.sqrt(MU / r)                       # km/s inertial
    T = 2 * math.pi * math.sqrt(r ** 3 / MU)    # s
    return dict(v_km_s=v, vg_km_s=v * RE / r, period_min=T / 60, orbits_per_day=86400 / T)


def corridor_width(eta_lo, eta_hi, h_km):
    return float(geom(eta_hi, h_km)[2] - geom(eta_lo, h_km)[2])


def mean_gap_days(width_km, n_sat, lat_deg, h_km):
    """Mean interval between looks at latitude lat for n_sat satellites each imaging a strip width_km wide on
    every ascending and descending pass (one side only). Coverage per day = 2 * n_orb * width per satellite."""
    n_orb = orbit(h_km)["orbits_per_day"]
    return 2 * math.pi * RE * math.cos(math.radians(lat_deg)) / (2 * n_orb * width_km * n_sat)


# ----------------------------------------------------------------------------------------------
# radar: beams, PRF, timing, resolution
# ----------------------------------------------------------------------------------------------
def beamwidth_deg(D_m, lam_m):
    return math.degrees(A["k_bw"] * lam_m / D_m)


def doppler_bw_hz(D_m, lam_m, v_m_s):
    """Doppler bandwidth of the full azimuth beam: 2 v theta_az / lambda."""
    return 2 * v_m_s * math.radians(beamwidth_deg(D_m, lam_m)) / lam_m


def stripmap_az_res_m(D_m):
    """Full-aperture single-look azimuth resolution after windowing: ~D/2 * k_win."""
    return D_m / 2 * A["k_win"]


def ground_range_res_m(bw_hz, inc_deg):
    return A["k_win"] * C / (2 * bw_hz * math.sin(math.radians(inc_deg)))


def echo_window_swath_km(prf_hz, duty, inc_deg):
    """Ground swath that fits between two transmit events at this PRF (single receive beam)."""
    pri = 1.0 / prf_hz
    tau = duty * pri
    usable = pri - 2 * tau - 2 * A["guard_us"] * 1e-6
    dr_slant_km = C * usable / 2 / 1e3
    return dr_slant_km / math.sin(math.radians(inc_deg)), dr_slant_km


def n_beams(eta_lo, eta_hi, bw_deg, crossover_db):
    """Elevation beams needed to tile the corridor with beams crossing at -crossover_db (Gaussian beam)."""
    width = bw_deg * math.sqrt(crossover_db / 3.0)       # full width at -x dB relative to the -3 dB width
    return int(math.ceil((eta_hi - eta_lo) / width))


# ----------------------------------------------------------------------------------------------
# NESZ (SAR radar equation for a distributed target)
# ----------------------------------------------------------------------------------------------
def nesz_db(p_avg_w, eta_deg, h_km, bw_hz, loss_db, a_eff_m2, lam_m, v_m_s):
    """NESZ = 4 (4 pi)^3 R^3 v sin(inc) k T0 F L B / (P_avg G^2 lambda^3 c), G = 4 pi A_eff / lambda^2.
    Equivalent to the team's '16 pi lambda R^3 v sin(inc) k T F L B / (P_avg A_eff^2 c)'."""
    R, inc, _ = geom(eta_deg, h_km)
    R = np.asarray(R) * 1e3
    th = np.deg2rad(inc)
    G = 4 * math.pi * a_eff_m2 / lam_m ** 2
    FL = 10 ** (loss_db / 10)
    nesz = 4 * (4 * math.pi) ** 3 * R ** 3 * v_m_s * np.sin(th) * KB * T0 * FL * bw_hz / (p_avg_w * G ** 2 * lam_m ** 3 * C)
    return 10 * np.log10(nesz)


def p_avg_for_nesz(nesz_target_db, eta_deg, h_km, bw_hz, loss_db, a_eff_m2, lam_m, v_m_s):
    """Average RF power giving nesz_target_db at eta_deg (NESZ is proportional to 1/P_avg)."""
    ref = nesz_db(1.0, eta_deg, h_km, bw_hz, loss_db, a_eff_m2, lam_m, v_m_s)
    return float(10 ** ((ref - nesz_target_db) / 10))


# ----------------------------------------------------------------------------------------------
# speckle
# ----------------------------------------------------------------------------------------------
def trigamma(x):
    """psi'(x) for x > 0: recurrence to x >= 6, then the asymptotic series."""
    shape = np.shape(x)
    x = np.array(x, dtype=float, ndmin=1)
    acc = np.zeros_like(x)
    xx = x.copy()
    while np.any(xx < 6):
        m = xx < 6
        acc[m] += 1.0 / xx[m] ** 2
        xx[m] += 1.0
    s = 1 / xx + 1 / (2 * xx ** 2) + 1 / (6 * xx ** 3) - 1 / (30 * xx ** 5) + 1 / (42 * xx ** 7) - 1 / (30 * xx ** 9)
    return (acc + s).reshape(shape)


def speckle_single_db(n_looks):
    """Std of 10 log10(I) for an N-look intensity (gamma law): DB * sqrt(psi'(N)). N=1 gives 5.57 dB; large N -> 4.34/sqrt(N)."""
    return DB * np.sqrt(trigamma(n_looks))


def speckle_change_db(n_looks):
    """Std of the dB difference of two independent N-look images: sqrt(2) times the single-image value (6.14/sqrt(N) for large N)."""
    return math.sqrt(2) * speckle_single_db(n_looks)


def looks_per_cell(cell_m, az_res_m, rg_res_m):
    return (cell_m / az_res_m) * (cell_m / rg_res_m)


def looks_for_change_noise(target_db):
    """Smallest N with speckle_change_db(N) <= target (search)."""
    n = np.arange(1, 20000)
    ok = speckle_change_db(n) <= target_db
    return int(n[ok][0])


def noise_with_thermal_db(n_looks, sigma0_db, nesz_db_val):
    """Approximate std (dB) of a noise-subtracted multi-look estimate: 6.14/sqrt(N) * (1 + NESZ/sigma0) (large N)."""
    x = 10 ** ((nesz_db_val - sigma0_db) / 10)
    return 6.14 / math.sqrt(n_looks) * (1 + x)


# ----------------------------------------------------------------------------------------------
# data, power, delta-v, eclipse, Faraday
# ----------------------------------------------------------------------------------------------
def data_rate_bps(bw_hz, n_channels=1):
    return A["n_pol"] * n_channels * 2 * A["bits"] * A["oversamp"] * bw_hz * A["rx_duty"]


def hohmann_dv_m_s(h1_km, h2_km):
    r1, r2 = RE + h1_km, RE + h2_km
    a = (r1 + r2) / 2
    v1 = math.sqrt(MU / r1)
    v2 = math.sqrt(MU / r2)
    vp = math.sqrt(MU * (2 / r1 - 1 / a))
    va = math.sqrt(MU * (2 / r2 - 1 / a))
    return 1e3 * (abs(vp - v1) + abs(v2 - va))


def perigee_lower_dv_m_s(h_km, hp_km):
    """Single burn at apogee (the circular orbit) to lower perigee to hp_km."""
    ra, rp = RE + h_km, RE + hp_km
    a = (ra + rp) / 2
    return 1e3 * (math.sqrt(MU / ra) - math.sqrt(MU * (2 / ra - 1 / a)))


def beta_angle_deg(day_of_year, inc_deg, ltan_h):
    """Sun-orbit beta angle for a sun-synchronous orbit with local time of ascending node ltan_h."""
    dec = math.radians(23.44) * math.sin(2 * math.pi * (day_of_year - 81) / 365.25)
    dra = math.radians((ltan_h - 12) * 15)        # RAAN minus solar right ascension
    i = math.radians(inc_deg)
    s = math.cos(dec) * math.sin(i) * math.sin(dra) + math.sin(dec) * math.cos(i)
    return math.degrees(math.asin(max(-1, min(1, s))))


def eclipse_fraction(beta_deg, h_km):
    r = RE + h_km
    cb = math.cos(math.radians(beta_deg))
    arg = math.sqrt(1 - (RE / r) ** 2) / cb if cb > 0 else 2
    if arg >= 1:
        return 0.0
    return math.acos(arg) / math.pi


def faraday_deg(f_hz, tec_tecu):
    """One-way Faraday rotation: 2.36e4 * B * cos(angle) * TEC / f^2 (SI units)."""
    return math.degrees(2.36e4 * A["b_field_t"] * A["cos_b"] * tec_tecu * 1e16 / f_hz ** 2)


# ----------------------------------------------------------------------------------------------
# published mission values (all read on ACCESS; 'unverified' marks values I could not read on a source page)
# ----------------------------------------------------------------------------------------------
SRC = {
    "esa_s1": "https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-1/Facts_and_figures",
    "sentiwiki_mission": "https://sentiwiki.copernicus.eu/web/s1-mission",
    "sentiwiki_products": "https://sentiwiki.copernicus.eu/web/s1-products",
    "eop_s1": "https://www.eoportal.org/satellite-missions/copernicus-sentinel-1",
    "nisar_quick": "https://nisar.jpl.nasa.gov/mission/quick-facts/",
    "nisar_sat": "https://science.nasa.gov/mission/nisar/about-the-satellite/",
    "nisar_numbers": "https://science.nasa.gov/blogs/nisar/2025/07/30/nasa-isros-nisar-mission-by-the-numbers/",
    "eop_nisar": "https://www.eoportal.org/satellite-missions/nisar",
    "jaxa_alos2": "https://www.eorc.jaxa.jp/ALOS/en/alos-2/a2_about_e.htm",
    "eop_alos2": "https://www.eoportal.org/satellite-missions/alos-2",
    "jaxa_alos4": "https://www.eorc.jaxa.jp/ALOS/en/alos-4/a4_about_e.htm",
    "jaxa_alos4_news": "https://www.satnavi.jaxa.jp/en/news/2024/09/03/9688/index.html",
    "eop_alos4": "https://www.eoportal.org/satellite-missions/alos-4",
    "eop_saocom": "https://www.eoportal.org/satellite-missions/saocom",
    "eop_rosel": "https://www.eoportal.org/satellite-missions/rose-l",
    "esa_rosel_contract": "https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Contract_signed_for_new_Copernicus_ROSE-L_mission",
    "eop_smap": "https://www.eoportal.org/satellite-missions/smap",
    "nsidc_smap_e": "https://nsidc.org/data/spl3smp_e/versions/6",
    "eop_novasar": "https://www.eoportal.org/satellite-missions/novasar-1",
    "eop_capella": "https://www.eoportal.org/satellite-missions/capella-x-sar",
    "esa_harmony": "https://www.esa.int/Applications/Observing_the_Earth/FutureEO/Contract_secures_build_for_ESA_s_Harmony_mission",
    "csa_rcm": "https://www.asc-csa.gc.ca/eng/satellites/radarsat/faq.asp",
    "evo1_faq": "https://essp.larc.nasa.gov/EVO-1/faq.html",
    "firefly": "https://fireflyspace.com/alpha/",
    "fcc_fr": "https://www.federalregister.gov/documents/2024/08/09/2024-17093/space-innovation-mitigation-of-orbital-debris-in-the-new-space-age",
    "oscar": "https://space.oscar.wmo.int/variables/view/soil_moisture_at_surface",
}

# name, band, alt km, repeat d per sat, fleet revisit note, wide-mode swath km, wide-mode res m, res kind,
# mass kg, peak RF kW, duty note, cost as published, source keys
MISSIONS = [
    ("Sentinel-1 (C/D)", "C 5.405 GHz", 693, 12, "6 d with two satellites", 250, 20, "single-look az (IW 5x20 m)", 2185,
     None, "SAR 25-30 min per orbit", "unverified", "esa_s1; sentiwiki_mission"),
    ("NISAR", "L 1.2575 GHz (+S)", 747, 12, "~6 d asc+desc average", 240, 7, "single-look az (3-10 m modes)", 2800,
     None, "26 Tbit/day plan; 80 TB/day products", "unverified", "nisar_quick; nisar_sat; nisar_numbers; eop_nisar"),
    ("SAOCOM-1A/1B", "L 1.275 GHz", 620, 16, "8 d with two satellites", 350, 50, "TOPSAR wide product (<50 m)", 3000,
     6.7, "15 % (about 15 min per orbit)", "unverified", "eop_saocom"),
    ("ALOS-2", "L 1.2575 GHz", 628, 14, "one satellite", 350, 100, "ScanSAR product (490 km at 60 m)", 2120,
     6.1, "50 % (49 min per orbit)", "unverified", "jaxa_alos2; eop_alos2"),
    ("ALOS-4", "L", 628, 14, "one satellite", 200, 3, "stripmap single-look (DBF); ScanSAR 700 km at 25 m", 3000,
     None, "n/a", "unverified", "jaxa_alos4; jaxa_alos4_news; eop_alos4"),
    ("ROSE-L (2028)", "L 1.2575 GHz", 693, 12, "6 d global land, 3 d Europe (2 sats)", 260, 10, "5-10 m geohazards; 5-50 m modes", 2223,
     None, ">35 min per orbit; 12.5 TB/day", "EUR 482 M, 2 sats, build", "eop_rosel; esa_rosel_contract"),
    ("SMAP radiometer", "L 1.41 GHz", 685, 8, "~3 d global", 1000, 36000, "radiometer; 9 km enhanced grid", 944,
     None, "continuous (spinning)", "unverified", "eop_smap; nsidc_smap_e"),
    ("NovaSAR-1", "S", 580, None, "one satellite", 150, 30, "ScanSAR product (20 m at 50-100 km)", 430,
     None, ">= 2 min per orbit (~2 %)", "GBP 21 M UK share only", "eop_novasar"),
    ("Capella Acadia", "X 9.4-9.9 GHz", 505, None, "fleet, tasked", 10, 1.2, "stripmap 100 x 10 km", 197,
     None, "n/a", "unverified", "eop_capella"),
    ("ARMO (proposed)", "L 1.2575 GHz", 620, None, "1.4 d mean gap (2 sats, simulated)", 362, 200, "product 100-200 m", 316,
     2.8, "1.5 % over CONUS", "USD 115 M request / 150 M estimate", "team report"),
]


# ----------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------
def write_csv(name, header, rows):
    path = os.path.join(DAT, name)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(["" if v is None else (f"{v:.6g}" if isinstance(v, float) else v) for v in r])
    print("wrote", os.path.relpath(path, ROOT))


PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.titlesize": 10, "axes.titleweight": "bold", "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF, "legend.frameon": False, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False})


# ----------------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------------
def main():
    h = A["h_km"]
    lam = C / A["f_hz"]
    D = A["D_m"]
    orb = orbit(h)
    v_ms = orb["v_km_s"] * 1e3
    a_phys = math.pi * D ** 2 / 4
    a_eff = A["eta_ap"] * a_phys
    loss_full = A["nf_loss_db"] + A["edge_db"] + A["scan_margin_db"]       # 12 dB, the team's full margin set
    loss_lean = A["nf_loss_db"]                                             # 6 dB, no mode margins
    p_avg_team = A["p_peak_w"] * A["duty"]                                  # 280 W

    # ---------------- geometry of the corridor ----------------
    R_lo, inc_lo, g_lo = geom(A["eta_min"], h)
    R_hi, inc_hi, g_hi = geom(A["eta_max"], h)
    corridor = float(g_hi - g_lo)
    print(f"\nOrbit {h:.0f} km: v = {orb['v_km_s']:.3f} km/s, ground speed {orb['vg_km_s']:.3f} km/s, "
          f"period {orb['period_min']:.1f} min, {orb['orbits_per_day']:.2f} orbits/day")
    print(f"Corridor {A['eta_min']:.0f}-{A['eta_max']:.0f} deg off-nadir: incidence {inc_lo:.1f}-{inc_hi:.1f} deg, "
          f"slant range {R_lo:.0f}-{R_hi:.0f} km, ground range {g_lo:.0f}-{g_hi:.0f} km -> width {corridor:.0f} km")

    # ---------------- beams, PRF, timing ----------------
    bw3 = beamwidth_deg(D, lam)
    bdop = doppler_bw_hz(D, lam, v_ms)
    prf = A["k_prf"] * bdop
    nb3 = n_beams(A["eta_min"], A["eta_max"], bw3, 3.0)
    nb15 = n_beams(A["eta_min"], A["eta_max"], bw3, 1.5)
    foot = {e: float(dg_deta(e, h) * math.radians(bw3)) for e in (25, 35, 45)}
    sw_mid, dr_mid = echo_window_swath_km(prf, A["duty"], float(geom(35, h)[1]))
    sw_near, _ = echo_window_swath_km(prf, A["duty"], float(inc_lo))
    sw_far, _ = echo_window_swath_km(prf, A["duty"], float(inc_hi))
    print(f"5.9 m dish: 3 dB beamwidth {bw3:.2f} deg; one elevation beam covers {foot[25]:.0f}/{foot[35]:.0f}/{foot[45]:.0f} km "
          f"at 25/35/45 deg; {nb3} beams at -3 dB crossover, {nb15} at -1.5 dB")
    print(f"Doppler bandwidth {bdop:.0f} Hz -> PRF {prf:.0f} Hz -> echo window fits {dr_mid:.0f} km of slant range = "
          f"{sw_near:.0f}/{sw_mid:.0f}/{sw_far:.0f} km of ground swath at the near/mid/far edge")
    # NISAR check of the PRF rule
    nis = orbit(747.0)
    print(f"check: the same PRF rule on NISAR's 12 m reflector gives {A['k_prf'] * doppler_bw_hz(12.0, lam, nis['v_km_s'] * 1e3):.0f} Hz "
          f"(published 1650 Hz)")

    # ---------------- mode table ----------------
    az_strip = stripmap_az_res_m(D)
    inc_mid = float(geom(35, h)[1])
    modes = []   # name, swath, n_b, az_res, nesz_penalty_db, n_channels
    modes.append(("Stripmap - one beam", sw_mid, 1, az_strip, 0.0, 1))
    modes.append((f"ScanSAR/TOPS - {nb3} beams (-3 dB crossover)", corridor, nb3, A["k_burst"] * nb3 * az_strip, 0.0, 1))
    modes.append((f"ScanSAR/TOPS - {nb15} beams (-1.5 dB crossover; team)", corridor, nb15, A["k_burst"] * nb15 * az_strip, 0.0, 1))
    modes.append((f"SweepSAR/DBF - {nb15} receive beams", corridor, nb15, az_strip, 10 * math.log10(nb15), nb15))
    rows = []
    for name, sw, nb, az, pen, nch in modes:
        rg20 = ground_range_res_m(20e6, inc_mid)
        rg40 = ground_range_res_m(40e6, inc_mid)
        rg80 = ground_range_res_m(80e6, inc_mid)
        n100 = looks_per_cell(100, az, rg20)
        n200 = looks_per_cell(200, az, rg20)
        n200_40 = looks_per_cell(200, az, rg40)
        n200_80 = looks_per_cell(200, az, rg80)
        nesz_far = float(nesz_db(p_avg_team, A["eta_max"], h, A["bw_hz"], loss_full, a_eff, lam, v_ms)) + pen
        p_for_req = p_avg_for_nesz(A["nesz_req_db"], A["eta_max"], h, A["bw_hz"], loss_full, a_eff, lam, v_ms) * 10 ** (pen / 10)
        rate = data_rate_bps(A["bw_hz"], nch)
        rows.append([name, round(sw), nb, round(az, 1), round(rg20, 1), round(rg40, 1), round(rg80, 1),
                     round(n100), round(n200), round(n200_40), round(n200_80),
                     round(float(speckle_change_db(max(n100, 1))), 2), round(float(speckle_change_db(max(n200, 1))), 2),
                     round(float(speckle_change_db(max(n200_40, 1))), 2), round(float(speckle_change_db(max(n200_80, 1))), 2),
                     round(nesz_far, 1), round(p_for_req), round(p_for_req / A["duty"] / 1e3, 2),
                     round(rate / 1e6), round(rate * A["duty_conus"] * 86400 / 8e9, 1),
                     round(rate * A["duty_conus"] * 86400 / A["xband_bps"] / 60, 1)])
    write_csv("feas_modes.csv",
              ["mode", "swath_km", "elevation_beams", "single_look_az_m", "single_look_ground_range_m_20MHz",
               "single_look_ground_range_m_40MHz", "single_look_ground_range_m_80MHz", "looks_100m_20MHz", "looks_200m_20MHz",
               "looks_200m_40MHz", "looks_200m_80MHz", "speckle_change_db_100m_20MHz", "speckle_change_db_200m_20MHz",
               "speckle_change_db_200m_40MHz", "speckle_change_db_200m_80MHz",
               "nesz_far_edge_db_at_280W_avg_12dB_losses", "p_avg_W_for_minus25dB", "p_peak_kW_for_minus25dB_10pct_duty",
               "data_rate_Mbit_s_20MHz", "GB_per_day_CONUS_1p5pct", "downlink_min_per_day_400Mbit_s"], rows)
    print("\nMode table (20 MHz, mid-corridor incidence %.0f deg):" % inc_mid)
    for r in rows:
        print(f"  {r[0]:<48s} swath {r[1]:>4} km  az {r[3]:>5} m  rg {r[4]:>4} m  looks/200 m {r[8]:>4}  "
              f"speckle change {r[12]:.2f} dB  NESZ far {r[15]:.1f} dB  P_peak for -25 dB {r[17]:.2f} kW")

    # ---------------- speckle vs product cell ----------------
    cells = np.array([10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 750, 1000], float)
    rg20_40 = ground_range_res_m(20e6, 40.0)
    rg80_40 = ground_range_res_m(80e6, 40.0)
    series = {
        f"ARMO ScanSAR {nb15} beams 20 MHz": (A["k_burst"] * nb15 * az_strip, rg20_40),
        f"ARMO ScanSAR {nb15} beams 80 MHz": (A["k_burst"] * nb15 * az_strip, rg80_40),
        "ARMO stripmap or SweepSAR 20 MHz": (az_strip, rg20_40),
        "NISAR 20 MHz (7 m az)": (7.0, ground_range_res_m(20e6, 40.0)),
        "Sentinel-1 IW (22 m az; 2.7-3.5 m slant range)": (22.0, 3.1 / math.sin(math.radians(39.0)) * 1.0),
    }
    rows = []
    for c_ in cells:
        row = [c_]
        for k, (az, rg) in series.items():
            n = looks_per_cell(c_, az, rg)
            row += [round(n, 1), round(float(speckle_change_db(max(n, 1.0))), 3)]
        rows.append(row)
    hdr = ["product_cell_m"]
    for k in series:
        hdr += [f"looks | {k}", f"speckle_change_db | {k}"]
    write_csv("feas_speckle.csv", hdr, rows)
    n_req = looks_for_change_noise(0.38)
    n_team = looks_for_change_noise(1.45)
    nis_box = looks_per_cell(600, 7.0, ground_range_res_m(20e6, 40.0))
    print(f"\nSpeckle-only change noise: N=1 {float(speckle_change_db(1)):.2f} dB; N=4 {float(speckle_change_db(4)):.2f} dB; "
          f"N=100 {float(speckle_change_db(100)):.2f} dB; 0.38 dB needs N >= {n_req}; 1.45 dB corresponds to N = {n_team}")
    print(f"The team's NISAR HH test averages 3x3 cells of 200 m (600 m box): ~{nis_box:.0f} looks at 20 MHz -> speckle "
          f"{float(speckle_change_db(nis_box)):.2f} dB, versus 1.45 dB measured")

    # thermal-noise comparison ScanSAR vs SweepSAR at equal average power (200 m cell, 20 MHz)
    n_scan = looks_per_cell(200, A["k_burst"] * nb15 * az_strip, rg20_40)
    n_sweep = looks_per_cell(200, az_strip, rg20_40)
    nesz_scan = float(nesz_db(p_avg_team, A["eta_max"], h, A["bw_hz"], loss_full, a_eff, lam, v_ms))
    nesz_sweep = nesz_scan + 10 * math.log10(nb15)
    cmp_rows = []
    for pol, s0 in (("HH", A["hh_db"]), ("HV", A["hv_db"])):
        cmp_rows.append([pol, s0, round(n_scan), round(nesz_scan, 1), round(noise_with_thermal_db(n_scan, s0, nesz_scan), 2),
                         round(n_sweep), round(nesz_sweep, 1), round(noise_with_thermal_db(n_sweep, s0, nesz_sweep), 2)])
    write_csv("feas_scan_vs_sweep.csv", ["pol", "sigma0_db_assumed", "scansar_looks_200m", "scansar_nesz_far_db",
                                         "scansar_noise_db_200m", "sweepsar_looks_200m", "sweepsar_nesz_far_db",
                                         "sweepsar_noise_db_200m"], cmp_rows)
    print("ScanSAR vs SweepSAR at the same 280 W average (200 m cell, speckle + thermal, far edge):")
    for r in cmp_rows:
        print(f"  {r[0]}: ScanSAR {r[4]:.2f} dB ({r[2]} looks, NESZ {r[3]} dB) | SweepSAR {r[7]:.2f} dB ({r[5]} looks, NESZ {r[6]} dB)")

    # ---------------- NESZ and power ----------------
    etas = np.arange(25.0, 45.01, 1.0)
    p_lean = p_avg_for_nesz(A["nesz_req_db"], A["eta_max"], h, A["bw_hz"], loss_lean, a_eff, lam, v_ms)
    p_full = p_avg_for_nesz(A["nesz_req_db"], A["eta_max"], h, A["bw_hz"], loss_full, a_eff, lam, v_ms)
    p_full_28 = p_avg_for_nesz(-28.0, A["eta_max"], h, A["bw_hz"], loss_full, a_eff, lam, v_ms)
    p_first = p_avg_for_nesz(A["nesz_req_db"], A["eta_max"], h, A["bw_hz"], A["nf_loss_db"] + A["scan_margin_db"], a_phys, lam, v_ms)
    p_full_40 = p_avg_for_nesz(A["nesz_req_db"], A["eta_max"], h, 40e6, loss_full, a_eff, lam, v_ms)
    p_full_7m = p_avg_for_nesz(A["nesz_req_db"], A["eta_max"], h, A["bw_hz"], loss_full, A["eta_ap"] * math.pi * 7.0 ** 2 / 4, lam, v_ms)
    recon = [
        ["100 % aperture efficiency, 9 dB (NF+losses 6 + ScanSAR 3), 20 MHz  [team's first model, 498 W peak]", round(p_first), round(p_first / A["duty"])],
        ["60 % aperture efficiency, 6 dB NF+losses only, 20 MHz", round(p_lean), round(p_lean / A["duty"])],
        ["60 % aperture efficiency, 12 dB (6 + 3 beam edge + 3 ScanSAR), 20 MHz  [team's revised model, 2.8 kW peak]", round(p_full), round(p_full / A["duty"])],
        ["same, but NESZ -28 dB (SAOCOM / ROSE-L class)", round(p_full_28), round(p_full_28 / A["duty"])],
        ["same as the revised model, but 40 MHz bandwidth", round(p_full_40), round(p_full_40 / A["duty"])],
        ["same as the revised model, but a 7.0 m reflector", round(p_full_7m), round(p_full_7m / A["duty"])],
        [f"SweepSAR with {nb15} beams, revised-model losses, 20 MHz", round(p_full * nb15), round(p_full * nb15 / A["duty"])],
    ]
    write_csv("feas_power_reconciliation.csv", ["case (NESZ -25 dB at 45 deg off-nadir unless stated)", "average_rf_W", "peak_rf_W_at_10pct_duty"], recon)
    print(f"\nAverage RF power for NESZ -25 dB at 45 deg off-nadir, 20 MHz: {p_first:.0f} W with a 100 % aperture and 9 dB "
          f"(team's first model: 498 W peak), {p_lean:.0f} W with 60 % efficiency and 6 dB NF+losses, {p_full:.0f} W with the "
          f"full 12 dB (team's 2.8 kW peak at 10 % duty = 280 W); -28 dB needs {p_full_28:.0f} W; 40 MHz needs {p_full_40:.0f} W; "
          f"a 7 m reflector needs {p_full_7m:.0f} W")
    nesz_rows = []
    for e in etas:
        R, inc, g = geom(e, h)
        n_scan_e = float(nesz_db(p_avg_team, e, h, A["bw_hz"], loss_full, a_eff, lam, v_ms))
        n_lean_e = float(nesz_db(p_avg_team, e, h, A["bw_hz"], loss_lean, a_eff, lam, v_ms))
        n_sweep_e = n_scan_e + 10 * math.log10(nb15)
        n_40 = float(nesz_db(p_avg_team, e, h, 40e6, loss_full, a_eff, lam, v_ms))
        nesz_rows.append([e, round(float(inc), 1), round(float(R)), round(float(g - g_lo)), round(n_scan_e, 1), round(n_lean_e, 1),
                          round(n_sweep_e, 1), round(n_40, 1)])
    write_csv("feas_nesz.csv", ["off_nadir_deg", "incidence_deg", "slant_range_km", "ground_range_from_near_edge_km",
                                "nesz_db_scansar_280W_12dB", "nesz_db_scansar_280W_6dB", f"nesz_db_sweepsar_{nb15}beams_280W_12dB",
                                "nesz_db_scansar_280W_12dB_40MHz"], nesz_rows)

    # ---------------- swath vs revisit ----------------
    lat = A["lat_ref"]
    gap_rows = []
    for w in (50, 100, 150, 200, 220, 250, 300, 350, 362, 400, 500, 575, 724):
        gap_rows.append([w, round(mean_gap_days(w, 1, lat, h), 2), round(mean_gap_days(w, 2, lat, h), 2),
                         round(mean_gap_days(w, 2, 0.0, h), 2), round(mean_gap_days(w, 2, 55.0, h), 2)])
    write_csv("feas_swath_gap.csv", ["usable_swath_km", "mean_gap_d_1sat_40N", "mean_gap_d_2sat_40N", "mean_gap_d_2sat_equator",
                                     "mean_gap_d_2sat_55N"], gap_rows)
    w80 = 2 * math.pi * RE * math.cos(math.radians(lat)) / (2 * orb["orbits_per_day"] * 2 * 1.46)
    print(f"\nMean gap at {lat:.0f} N (formula): 1 sat x 362 km {mean_gap_days(362, 1, lat, h):.2f} d (team design A: 2.8 d); "
          f"1 sat x 575 km {mean_gap_days(575, 1, lat, h):.2f} d (team design C: 1.8 d); 2 sats x 362 km "
          f"{mean_gap_days(362, 2, lat, h):.2f} d (team design D: 1.4 d)")
    print(f"Two satellites: 220/250/300 km usable -> {mean_gap_days(220, 2, lat, h):.2f}/{mean_gap_days(250, 2, lat, h):.2f}/"
          f"{mean_gap_days(300, 2, lat, h):.2f} d; the 1.46 d needed for 80 % capture needs {w80:.0f} km imaged on every pass")

    # incidence-angle choices and altitude: how wide is the usable band?
    inc_rows = []
    for hh in (550.0, 620.0, 700.0, 800.0):
        for lo, hi in ((30, 45), (25, 45), (25, 50), (28.0, 51.0), (20, 55)):
            e_lo, e_hi = eta_from_inc(lo, hh), eta_from_inc(hi, hh)
            w = corridor_width(float(e_lo), float(e_hi), hh)
            Rf = float(geom(float(e_hi), hh)[0])
            pen = 10 * math.log10((Rf ** 3 * math.sin(math.radians(hi))) / (float(R_hi) ** 3 * math.sin(math.radians(float(inc_hi)))))
            inc_rows.append([hh, lo, hi, round(float(e_lo), 1), round(float(e_hi), 1), round(w), round(mean_gap_days(w, 2, lat, hh), 2),
                             round(mean_gap_days(w, 1, lat, hh), 2), round(Rf), round(pen, 1)])
    write_csv("feas_incidence_altitude.csv", ["altitude_km", "incidence_min_deg", "incidence_max_deg", "off_nadir_min_deg",
                                              "off_nadir_max_deg", "usable_width_km", "mean_gap_d_2sat_40N", "mean_gap_d_1sat_40N",
                                              "far_edge_slant_range_km", "far_edge_nesz_penalty_db_vs_ARMO_45deg"], inc_rows)
    print("Usable width by incidence range (620 km): " + "; ".join(
        f"{r[1]:.0f}-{r[2]:.0f} deg = {r[5]} km" for r in inc_rows if r[0] == 620.0))

    # one satellite matching two: altitude at which the 25-45 deg corridor is 2 x 362 km wide
    hs = np.arange(620, 1600, 5.0)
    widths = np.array([corridor_width(A["eta_min"], A["eta_max"], float(x)) for x in hs])
    h_one = float(hs[np.argmin(np.abs(widths - 2 * corridor))])
    R_one = float(geom(A["eta_max"], h_one)[0])
    pen_one = 10 * math.log10(R_one ** 3 / float(R_hi) ** 3)
    print(f"One satellite with the same 25-45 deg corridor needs ~{h_one:.0f} km altitude for a {2 * corridor:.0f} km corridor: "
          f"slant range {R_one:.0f} km, NESZ penalty {pen_one:.1f} dB, Hohmann from 620 km {hohmann_dv_m_s(620, h_one):.0f} m/s")

    # ---------------- budgets: delta-v, power, eclipse, Faraday ----------------
    dv_raise = hohmann_dv_m_s(500, h)
    dv_deorbit_300 = perigee_lower_dv_m_s(h, 300)
    dv_to_500 = hohmann_dv_m_s(h, 500)
    dv_phase, dv_drag = 10.0, 10.0      # ASSUMPTIONS: slow 180 deg phasing via a small altitude offset; drag make-up over 3 yr at 620 km
    dv_total = dv_raise + dv_phase + dv_drag + dv_deorbit_300
    prop = A["dry_mass_kg"] * (math.exp(dv_total / (9.80665 * A["isp_s"])) - 1)
    dc_img = p_avg_team / A["sspa_eff"] + A["radar_elec_w"]
    rate20 = data_rate_bps(A["bw_hz"])
    dl_min_conus = rate20 * A["duty_conus"] * 86400 / A["xband_bps"] / 60
    orb_avg_conus = A["bus_w"] + dc_img * A["duty_conus"] + 60.0 * dl_min_conus * 60 / 86400
    orb_avg_global = A["bus_w"] + dc_img * A["duty_global"] + 60.0 * (dl_min_conus * A["duty_global"] / A["duty_conus"]) * 60 / 86400
    array = lambda w: w * 1.3 / (1361 * 0.28 * 0.85 * 0.85)
    days = np.arange(1, 366)
    betas = np.array([beta_angle_deg(d, 97.87, 18.0) for d in days])
    ecl = np.array([eclipse_fraction(b, h) for b in betas])
    ecl_days = int((ecl > 0).sum())
    ecl_max_min = float(ecl.max() * orb["period_min"])
    far_dawn, far_pm = faraday_deg(A["f_hz"], A["tec_dawn"]), faraday_deg(A["f_hz"], A["tec_pm"])
    budget = [
        ["delta-v: raise 500 -> 620 km (rideshare insertion), m/s", round(dv_raise, 1), "Hohmann, calculation"],
        ["delta-v: 180 deg phasing of the second satellite, m/s", dv_phase, "ASSUMPTION (slow drift, small altitude offset)"],
        ["delta-v: drag make-up at 620 km over 3 years, m/s", dv_drag, "ASSUMPTION (order of magnitude)"],
        ["delta-v: deorbit burn, 620 km circular -> 300 km perigee, m/s", round(dv_deorbit_300, 1), "calculation (single burn)"],
        ["delta-v: lower to a 500 km circular orbit instead, m/s", round(dv_to_500, 1), "Hohmann, calculation (cheaper, decay time not computed)"],
        ["delta-v total with the 300 km deorbit, m/s", round(dv_total, 1), "sum; team assumed 120 m/s"],
        ["propellant for that total at Isp 220 s on 290 kg dry, kg", round(prop, 1), "rocket equation"],
        ["average RF power while imaging, W", round(p_avg_team), "2.8 kW peak x 10 % duty"],
        ["DC power while imaging, W", round(dc_img), "35 % DC-to-RF + 150 W electronics"],
        ["orbit-average power, CONUS duty 1.5 %, W", round(orb_avg_conus), "calculation with the team's bus 150 W"],
        ["solar array for that, m2", round(array(orb_avg_conus), 2), "28 % cells, 0.85 packing, 0.85 EOL, 1.3 margin"],
        ["orbit-average power if global-land duty 30 %, W", round(orb_avg_global), "calculation"],
        ["solar array for that, m2", round(array(orb_avg_global), 2), "same factors"],
        ["raw data rate while imaging, 20 MHz, HH+HV, Mbit/s", round(rate20 / 1e6), "team's formula"],
        ["data per day, CONUS duty 1.5 %, GB", round(rate20 * A["duty_conus"] * 86400 / 8e9, 1), "calculation"],
        ["downlink minutes per day at 400 Mbit/s", round(dl_min_conus, 1), "calculation"],
        ["data per day if 40 MHz and global-land duty 30 %, GB", round(2 * rate20 * A["duty_global"] * 86400 / 8e9), "calculation"],
        ["days per year with eclipses (LTAN 18:00, 620 km)", ecl_days, "beta-angle calculation"],
        ["longest eclipse per orbit, min", round(ecl_max_min, 1), "beta-angle calculation"],
        ["Faraday rotation at 06:00 local (TEC 10 TECU), deg one-way", round(far_dawn, 1), "calculation, assumed TEC"],
        ["Faraday rotation mid-afternoon near solar maximum (TEC 50 TECU), deg", round(far_pm, 1), "calculation, assumed TEC"],
        ["orbital period at 400 / 620 / 800 km, min", f"{orbit(400)['period_min']:.1f} / {orbit(620)['period_min']:.1f} / {orbit(800)['period_min']:.1f}", "calculation"],
        ["delta-v: deorbit burn from 800 km to a 300 km perigee, m/s", round(perigee_lower_dv_m_s(800, 300), 1), "calculation"],
        ["elevation beams for 20-55 deg incidence (-3 dB / -1.5 dB crossover)", f"{n_beams(float(eta_from_inc(20, h)), float(eta_from_inc(55, h)), bw3, 3.0)} / {n_beams(float(eta_from_inc(20, h)), float(eta_from_inc(55, h)), bw3, 1.5)}", "calculation"],
        ["elevation beams for 20-50 deg off-nadir, design C (-3 dB / -1.5 dB crossover)", f"{n_beams(20.0, 50.0, bw3, 3.0)} / {n_beams(20.0, 50.0, bw3, 1.5)}", "calculation"],
        ["Firefly Alpha capacity to 500 km SSO, kg", 630, SRC["firefly"]],
        ["two ARMO satellites + 40 kg adapter, kg", 2 * 316 + 40, "team"],
    ]
    write_csv("feas_budget.csv", ["item", "value", "basis"], budget)
    print(f"\nDelta-v {dv_total:.0f} m/s total -> {prop:.0f} kg propellant; imaging DC {dc_img:.0f} W; orbit-average "
          f"{orb_avg_conus:.0f} W (CONUS) / {orb_avg_global:.0f} W (global); eclipses on {ecl_days} days/yr, max {ecl_max_min:.0f} min; "
          f"Faraday {far_dawn:.1f} deg (dawn) vs {far_pm:.1f} deg (afternoon, solar max)")

    # ---------------- mission table ----------------
    mrows = []
    for m in MISSIONS:
        mrows.append(list(m[:12]) + ["; ".join(SRC.get(k.strip(), k.strip()) for k in m[12].split(";")), ACCESS])
    write_csv("feas_missions.csv", ["mission", "band", "altitude_km", "repeat_days_per_satellite", "fleet_revisit", "wide_mode_swath_km",
                                    "wide_mode_resolution_m", "resolution_kind", "mass_kg", "peak_rf_kw", "duty_cycle",
                                    "cost_as_published", "sources", "accessed"], mrows)

    # ==========================================================================================
    # figures
    # ==========================================================================================
    # (a) swath vs single-look azimuth resolution
    fig, ax = plt.subplots(figsize=(7.6, 5.2), dpi=150)
    nbs = np.arange(1, 13)
    sw_curve = nbs * (corridor / nb3)
    az_curve = np.where(nbs == 1, az_strip, A["k_burst"] * nbs * az_strip)
    ax.plot(sw_curve, az_curve, color=PAL[0], lw=2, label="ARMO ScanSAR/TOPS: beams added one by one")
    ax.plot([sw_mid, corridor], [az_strip, az_strip], color=PAL[1], lw=2, ls="--",
            label=f"ARMO SweepSAR: full resolution, NESZ {10 * math.log10(nb15):.0f} dB worse (or {nb15}x power)")
    ax.scatter([sw_mid], [az_strip], s=60, color=PAL[0], zorder=5)
    ax.annotate("ARMO stripmap (one beam)", (sw_mid, az_strip), xytext=(-112, -40), textcoords="offset points", fontsize=8, color=INK, ha="left",
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
    ax.scatter([corridor], [A["k_burst"] * nb15 * az_strip], s=60, color=PAL[0], zorder=5)
    ax.annotate(f"ARMO ScanSAR {nb15} beams:\n{corridor:.0f} km, {A['k_burst'] * nb15 * az_strip:.0f} m", (corridor, A["k_burst"] * nb15 * az_strip),
                xytext=(-70, 40), textcoords="offset points", fontsize=8, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
    pts_sl = [("S1 SM", 80, 5), ("S1 IW", 250, 22), ("S1 EW", 400, 43), ("NISAR", 240, 7), ("ALOS-2 SM", 50, 3),
              ("ALOS-4 SM (DBF)", 200, 3), ("ROSE-L", 260, 10), ("NovaSAR SM", 20, 6), ("Capella SM", 10, 1.2)]
    pts_pr = [("ALOS-2 ScanSAR", 350, 100), ("ALOS-2 ScanSAR wide", 490, 60), ("ALOS-4 ScanSAR", 700, 25), ("SAOCOM TW", 350, 50),
              ("SAOCOM TN", 150, 30), ("SAOCOM SM", 40, 10), ("NovaSAR ScanSAR", 100, 20), ("SMAP radar (lost 2015)", 1000, 3000)]
    ax.scatter([p[1] for p in pts_sl], [p[2] for p in pts_sl], s=34, color=PAL[2], zorder=4, label="published single-look resolution")
    ax.scatter([p[1] for p in pts_pr], [p[2] for p in pts_pr], s=34, facecolors="none", edgecolors=PAL[2], linewidths=1.5, zorder=4,
               label="published mode/product resolution (multi-looked)")
    offsets = {"S1 IW": (5, -10), "ROSE-L": (-6, 6), "ALOS-4 SM (DBF)": (5, -9), "SAOCOM TW": (6, -9), "S1 EW": (6, 2),
               "ALOS-2 SM": (5, -10), "ALOS-2 ScanSAR": (-82, 4), "ALOS-2 ScanSAR wide": (6, -4), "NISAR": (5, -9), "SMAP radar (lost 2015)": (-95, -12),
               "NovaSAR SM": (-8, -11), "S1 EW": (7, 5), "SAOCOM TW": (4, -11)}
    for name, x, y in pts_sl + pts_pr:
        dx, dy = offsets.get(name, (5, 3))
        ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=7.5, color=INK2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(8, 1500)
    ax.set_ylim(0.8, 6000)
    ax.set_xlabel("imaged swath width (km)")
    ax.set_ylabel("azimuth resolution (m)")
    ax.set_title("Swath against azimuth resolution: ARMO's 5.9 m reflector at 620 km and real missions")
    ax.legend(loc="upper left", fontsize=8)
    fig.text(0.01, 0.012, "ARMO curves: calculation (stripmap 0.6 D; ScanSAR azimuth = 1.7 x beams x stripmap; one beam = 45 km of ground).\n"
             "Points: ESA Sentinel-1 facts and sentiwiki, NASA NISAR, JAXA ALOS-2/-4, eoPortal SAOCOM / ROSE-L / NovaSAR-1 / Capella / SMAP, read " + ACCESS + ".",
             fontsize=6.5, color=INK2)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(os.path.join(FIG, "feas_swath_vs_resolution.png"))
    plt.close(fig)

    # (b) speckle change noise vs product cell
    fig, ax = plt.subplots(figsize=(7.6, 5.0), dpi=150)
    cells_f = np.logspace(1, 3, 60)
    styles = [(PAL[0], "-"), (PAL[0], "--"), (PAL[1], "-"), (PAL[2], "-"), (PAL[3], "-")]
    for (k, (az, rg)), (col, ls) in zip(series.items(), styles):
        n = looks_per_cell(cells_f, az, rg)
        y = speckle_change_db(np.maximum(n, 1.0))
        ax.plot(cells_f, y, color=col, ls=ls, lw=2, label=k)
    ax.axhline(0.38, color=PAL[7], lw=1.2, ls=":")
    ax.text(10.5, 0.40, "team's target 0.38 dB (about 0.038 m3/m3)", fontsize=8, color=PAL[7], va="bottom", backgroundcolor=SURF)
    ax.axhline(1.45, color=INK2, lw=1.2, ls=":")
    ax.text(130, 1.5, "measured no-change noise, NISAR HH: 1.45 dB", fontsize=8, color=INK2, va="bottom", backgroundcolor=SURF)
    ax.scatter([600], [float(speckle_change_db(nis_box))], s=50, color=PAL[2], zorder=5)
    ax.annotate("team's NISAR HH test box\n(3 x 3 cells of 200 m): speckle alone\n= %.2f dB, measured 1.45 dB" % float(speckle_change_db(nis_box)),
                (600, float(speckle_change_db(nis_box))), xytext=(-150, -34), textcoords="offset points", fontsize=7.5, color=INK2,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10, 1000)
    ax.set_ylim(0.05, 10)
    ax.set_xlabel("product cell size (m)")
    ax.set_ylabel("speckle-only std of a backscatter change (dB)")
    ax.set_title("Speckle noise of a two-date change against product cell size")
    ax.legend(loc="lower left", fontsize=7.5)
    fig.text(0.01, 0.012, "Calculation: std = 4.34 sqrt(2 psi'(N)) dB for N independent looks per cell (gamma speckle);\n"
             "N = cell^2 / (azimuth resolution x ground-range resolution) at 40 deg incidence. Published resolutions for NISAR and Sentinel-1.",
             fontsize=6.5, color=INK2)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(os.path.join(FIG, "feas_speckle_vs_cell.png"))
    plt.close(fig)

    # (c) NESZ across the corridor and peak power against usable swath
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 4.6), dpi=150)
    incs = np.array([float(geom(e, h)[1]) for e in etas])
    ax1.plot(incs, [float(nesz_db(p_avg_team, e, h, A["bw_hz"], loss_full, a_eff, lam, v_ms)) for e in etas], color=PAL[0], lw=2,
             label="ScanSAR/TOPS, 280 W average RF, 12 dB NF + losses + margins")
    ax1.plot(incs, [float(nesz_db(p_avg_team, e, h, A["bw_hz"], loss_lean, a_eff, lam, v_ms)) for e in etas], color=PAL[0], lw=2, ls="--",
             label="same, 6 dB NF + losses only (no mode margins)")
    ax1.plot(incs, [float(nesz_db(p_avg_team, e, h, A["bw_hz"], loss_full, a_eff, lam, v_ms)) + 10 * math.log10(nb15) for e in etas],
             color=PAL[1], lw=2, label=f"SweepSAR, {nb15} receive beams, same 280 W")
    ax1.axhline(-25, color=PAL[7], lw=1.2, ls=":")
    ax1.text(incs[0] + 0.3, -24.7, "team's requirement -25 dB", fontsize=8, color=PAL[7])
    ax1.axhline(-28, color=INK2, lw=1.0, ls=":")
    ax1.text(incs[0] + 0.3, -27.7, "SAOCOM/ROSE-L class -28 dB", fontsize=8, color=INK2)
    ax1.set_xlabel("incidence angle across the corridor (deg)")
    ax1.set_ylabel("NESZ (dB)")
    ax1.set_title("Noise floor across the 25-45 deg corridor")
    ax1.legend(fontsize=7, loc="lower right")
    ax1.set_ylim(-42, -12)
    sw_axis = np.arange(50, 701, 10.0)
    nb_axis = np.maximum(1, sw_axis / (corridor / nb15))
    ax2.plot(sw_axis, np.full_like(sw_axis, p_full / A["duty"] / 1e3), color=PAL[0], lw=2, label="ScanSAR/TOPS, NESZ -25 dB")
    ax2.plot(sw_axis, np.full_like(sw_axis, p_full_28 / A["duty"] / 1e3), color=PAL[0], lw=2, ls="--", label="ScanSAR/TOPS, NESZ -28 dB")
    ax2.plot(sw_axis, nb_axis * p_full / A["duty"] / 1e3, color=PAL[1], lw=2, label="SweepSAR, NESZ -25 dB (power grows with beams)")
    ax2.scatter([corridor], [A["p_peak_w"] / 1e3], s=60, color=INK, zorder=5)
    ax2.annotate("team: 2.8 kW at 362 km", (corridor, A["p_peak_w"] / 1e3), xytext=(6, 6), textcoords="offset points", fontsize=8, color=INK2)
    ax2.axvspan(220, 300, color=GRID, alpha=0.8, lw=0)
    ax2.text(222, 0.62, "usable 220-300 km\n(30-45 deg incidence)", fontsize=7.5, color=INK2)
    ax2.set_xlabel("usable imaged swath (km)")
    ax2.set_ylabel("peak RF power at 10 % duty (kW)")
    ax2.set_title("Peak RF power for the far edge (5.9 m, 20 MHz)")
    ax2.set_yscale("log")
    ax2.set_ylim(0.5, 60)
    ax2.set_yticks([0.5, 1, 2, 5, 10, 20, 50])
    ax2.set_yticklabels(["0.5", "1", "2", "5", "10", "20", "50"])
    ax2.legend(fontsize=7, loc="upper left")
    fig.text(0.01, 0.012, "Calculation: NESZ = 4 (4 pi)^3 R^3 v sin(inc) k T0 F L B / (P_avg G^2 lambda^3 c), G = 4 pi (0.6 A) / lambda^2, T0 = 290 K,\n"
             "B = 20 MHz, duty 10 %, 620 km. SweepSAR: transmit gain spread over all receive beams (10 log10 of the beam count).", fontsize=6.5, color=INK2)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(FIG, "feas_nesz_power.png"))
    plt.close(fig)

    # (d) mission comparison
    names = [m[0] for m in MISSIONS]
    fig, axes = plt.subplots(1, 4, figsize=(12.8, 5.0), dpi=150, sharey=True)
    y = np.arange(len(names))[::-1]
    rep = [m[3] if m[3] else np.nan for m in MISSIONS]
    sw = [m[5] for m in MISSIONS]
    res = [m[6] for m in MISSIONS]
    axes[0].barh(y, rep, color=PAL[0], height=0.5)
    axes[0].set_title("exact repeat per satellite (d)")
    for yi, m, r in zip(y, MISSIONS, rep):
        axes[0].text((r if r == r else 0) + 0.6, yi, m[4], fontsize=6.5, color=INK2, va="center")
    axes[0].set_xlim(0, 30)
    axes[1].barh(y, sw, color=PAL[1], height=0.5)
    axes[1].set_xscale("log")
    axes[1].set_xlim(5, 3000)
    axes[1].set_title("widest routine swath (km)")
    axes[2].barh(y, res, color=PAL[2], height=0.5)
    axes[2].set_xscale("log")
    axes[2].set_xlim(0.5, 100000)
    axes[2].set_title("resolution of that mode (m)")
    for yi, m in zip(y, MISSIONS):
        axes[2].text(m[6] * 1.4, yi, m[7], fontsize=6, color=INK2, va="center")
    axes[3].barh(y, [m[8] for m in MISSIONS], color=PAL[3], height=0.5)
    axes[3].set_xscale("log")
    axes[3].set_xlim(30, 10000)
    axes[3].set_title("mass (kg) and published cost")
    for yi, m in zip(y, MISSIONS):
        axes[3].text(m[8] * 1.3, yi, m[11], fontsize=6, color=INK2, va="center")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(names, fontsize=8)
    for ax_ in axes:
        ax_.grid(axis="y", visible=False)
    fig.suptitle("Soil-moisture-relevant radar missions compared with ARMO (values as published; see feas_missions.csv for sources)",
                 fontsize=10, fontweight="bold")
    fig.text(0.01, 0.01, "Sources read " + ACCESS + ": ESA, Copernicus sentiwiki, NASA/JPL, JAXA, eoPortal (SAOCOM, ROSE-L, SMAP, NovaSAR-1, Capella). "
             "'unverified' = no cost figure found on an agency page. ARMO values are the team's proposal.", fontsize=6.5, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(os.path.join(FIG, "feas_mission_comparison.png"))
    plt.close(fig)
    print("\nfigures written to", os.path.relpath(FIG, ROOT))


if __name__ == "__main__":
    main()
