"""Mass, delta-v, launch and budget numbers for the ARMO report (8 Nov 2026).

Inputs are the team's realistic sizing (code/06_mission_design.py, realistic_sizing):
  dry mass current best estimate 230 kg, 30 % growth margin, hydrazine Isp 220 s.
New here: a rideshare launch (SpaceX published price list, fetched 8 Nov 2026 from
api-rideshare.spacex.com/api/v2/price: US$ 7,000/kg in 2026 and 7,500/kg in 2027; plate limits
quarter 300 kg, half 452 kg, full 831 kg; listed SSO flights 500-600 km) and a delta-v budget for
raising from a 500 km insertion (worst case) to the 620 km operating orbit, drag make-up and a
controlled deorbit to a 300 km perigee. (Set H = 620.0 for the exact-repeat alternative: 213 m/s, 330 kg wet.)
Run: python launch_mass_budget.py  (prints everything quoted in the report)
"""
import math

MU, RE = 398600.4418, 6371.0
G0 = 9.80665


def v_circ(h):
    return math.sqrt(MU / (RE + h))


def hohmann(h1, h2):
    r1, r2 = RE + h1, RE + h2
    a = (r1 + r2) / 2
    dv1 = math.sqrt(MU * (2 / r1 - 1 / a)) - math.sqrt(MU / r1)
    dv2 = math.sqrt(MU / r2) - math.sqrt(MU * (2 / r2 - 1 / a))
    return 1e3 * (abs(dv1) + abs(dv2))


def deorbit(h, hp=300.0):
    r1, r2 = RE + h, RE + hp
    a = (r1 + r2) / 2
    return 1e3 * (math.sqrt(MU / r1) - math.sqrt(MU * (2 / r1 - 1 / a)))


def drift_deg_per_orbit(h_low, h_high):
    n1 = math.sqrt(MU / (RE + h_low) ** 3); n2 = math.sqrt(MU / (RE + h_high) ** 3)
    T = 2 * math.pi / n2
    return math.degrees((n1 - n2) * T)


def sso_inclination(h):
    J2, R_EQ = 1.08262668e-3, 6378.137
    a = R_EQ + h; rho = 2 * math.pi / (365.2422 * 86400)
    return math.degrees(math.acos(-(2 * rho * a ** 3.5) / (3 * J2 * R_EQ ** 2 * math.sqrt(MU))))


def raise_with_trim(h1, h2):
    """Hohmann raise with the sun-synchronous inclination change done in the second burn."""
    r1, r2 = RE + h1, RE + h2; a = (r1 + r2) / 2
    v1 = math.sqrt(MU / r1); vp = math.sqrt(MU * (2 / r1 - 1 / a)); va = math.sqrt(MU * (2 / r2 - 1 / a)); v2 = math.sqrt(MU / r2)
    di = math.radians(sso_inclination(h2) - sso_inclination(h1))
    return 1e3 * ((vp - v1) + math.sqrt(va ** 2 + v2 ** 2 - 2 * va * v2 * math.cos(di))), math.degrees(di)


H = 620.0
dry_cbe = 230.0
dry = dry_cbe * 1.30
raise_dv, trim_deg = raise_with_trim(500.0, H)
drag_dv = 10.0
deorb_dv = deorbit(H)
sub = raise_dv + drag_dv + deorb_dv
total = 1.10 * sub
prop = dry * (math.exp(total / (G0 * 220.0)) - 1)
wet = dry + prop
print(f"raise 500->{H:.0f} km with {trim_deg:.2f} deg inclination trim: {raise_dv:.1f} m/s (raise alone {hohmann(500.0, H):.1f}); drag {drag_dv:.0f}; deorbit to 300 km perigee: {deorb_dv:.1f}; sum {sub:.1f}; +10 % = {total:.0f} m/s")
print(f"dry {dry:.0f} kg, propellant {prop:.1f} kg, wet {wet:.0f} kg, pair {2*wet:.0f} kg")
print(f"half plate 452 kg -> margin {452 - wet:.0f} kg per satellite")
d = drift_deg_per_orbit(500.0, H)
print(f"phasing: satellite B waits at 500 km while A raises: {d:.1f} deg per orbit -> 180 deg in {180/d:.0f} orbits ({180/d*2*math.pi*math.sqrt((RE+H)**3/MU)/86400:.1f} days)")

# launch price at the published per-kg rate, escalated at the list's +500 US$/kg per year to a 2030 launch
for yr, rate in ((2026, 7000), (2027, 7500), (2030, 9000)):
    print(f"{yr}: {rate} US$/kg x {2*wet:.0f} kg = US$ {rate*2*wet/1e6:.2f} M")
sep_2026 = 629_600  # 24-inch separation system, 2026 list (US$)
print(f"two 24-inch separation systems at the 2026 list: US$ {2*sep_2026/1e6:.2f} M")
launch_budget = 9000 * 2 * wet / 1e6 + 2 * sep_2026 * 1.2 / 1e6
print(f"launch line incl. separation systems escalated 20 %: US$ {launch_budget:.1f} M -> budget 8.0")

# budget: team allocation with the launch line changed; freed money to reserves
lines = [("Two spacecraft and radar payloads", 50.0, "AD"), ("Integration, test and qualification", 8.0, "AD"),
         ("Launch (rideshare, both satellites)", 8.0, "AD"), ("Ground systems and data processing", 5.0, "AD"),
         ("Project management, systems engineering, mission assurance", 7.0, "AD"),
         ("Science team, algorithms and validation", 3.0, "mixed"), ("Mission operations and open data (1-year prime)", 5.0, "E")]
base = sum(x[1] for x in lines)
res = 115.0 - base
sci_ad, sci_e = 1.5, 1.5
ad = sum(x[1] for x in lines if x[2] == "AD") + sci_ad
e = 5.0 + sci_e
res_e = 1.5
res_ad = res - res_e
print(f"base {base:.1f}, reserves {res:.1f} ({100*res/115:.1f} % of total)")
print(f"A-D cost excl. reserves {ad:.1f}; A-D reserves {res_ad:.1f} = {100*res_ad/ad:.0f} % (AO minimum 25 %); E {e:.1f}, E reserves {res_e:.1f} = {100*res_e/e:.0f} %")

# threshold mission: one satellite of the baseline radar (design A), scaled with the team's first/second-unit split (45 + 25 of 70)
fs1 = 50.0 * 45 / 70
it1 = 8.0 * 45 / 70
launch1 = 8.0 / 2
other = 5.0 + 7.0 + 3.0 + 5.0
base1 = fs1 + it1 + launch1 + other
ad1 = fs1 + it1 + launch1 + 5.0 + 7.0 + sci_ad
tot1 = base1 + 0.35 * ad1 + res_e
print(f"threshold mission (1 satellite): base {base1:.1f}, with 35 % A-D reserves {tot1:.1f} US$ M")
# higher-price threshold: the first unit at the team's earlier concept figure (45 of 70)
fsh = 45.0; adh = fsh + 0.16 * fsh + 8.0 / 2 + 5.0 + 0.14 * fsh + sci_ad
toth = adh + (5.0 + sci_e) + 0.35 * adh + res_e
print(f"threshold mission at the earlier US$ 70 M estimate (first unit 45): {toth:.1f} US$ M")

# commercial X-band tasking cost for 120 stations
for days in (137, 214):
    looks = days / 1.42
    lo, hi = looks * 120 * 1320 / 1e6, looks * 120 * 1385 / 1e6
    print(f"X-band tasking, {days} days: {looks:.0f} looks per station -> US$ {lo:.1f}-{hi:.1f} M")
