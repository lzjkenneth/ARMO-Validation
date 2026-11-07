"""Step 4 - independent orbit simulator for SAR look times over the USCRN stations.

Physics: circular sun-synchronous orbit, J2 secular rates (node precession, argument-of-latitude rate), Earth rotation
(GMST), WGS84 station positions. A side-looking SAR "looks" at a station at zero Doppler (the moment the station is abeam,
i.e. the line of sight is perpendicular to the velocity) if the station's off-nadir angle is inside the access range on the
looking side. Positions every 10 s. Each design is averaged over N_PHASE starting positions along the orbit.

The geometric access width of an off-nadir range [t1, t2] at altitude h is reported with Earth curvature:
  incidence eta = asin((Re + h) / Re * sin(theta)), central angle gamma = eta - theta, ground range = Re * gamma.
This is ACCESS, not a validated imaging swath (see docs/RESOLUTION_AND_FEASIBILITY.md).
"""
import os, json, math
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROC = os.path.join(ROOT, "data", "processed")
MU, RE, J2, RE_MEAN = 398600.4418, 6378.137, 1.08262668e-3, 6371.0088
OMEGA_E = 7.2921150e-5                      # Earth rotation rad/s
RHO_SSO = 2 * math.pi / (365.2422 * 86400)  # required node precession rate for sun-synchronism
J2000 = pd.Timestamp("2000-01-01 12:00", tz="UTC")

def sso_inclination(h_km):
    a = RE + h_km
    c = -(2 * RHO_SSO * a ** 3.5) / (3 * J2 * RE ** 2 * math.sqrt(MU))
    return math.degrees(math.acos(c))

def ground_range_km(h_km, theta_deg):
    th = np.radians(theta_deg); eta = np.arcsin(np.clip((RE_MEAN + h_km) / RE_MEAN * np.sin(th), -1, 1))
    return RE_MEAN * (eta - th)

def access_width_km(h_km, t1, t2):
    return float(ground_range_km(h_km, t2) - ground_range_km(h_km, t1))

def stations_ecef(lat, lon):
    a, f = 6378.137, 1 / 298.257223563; e2 = f * (2 - f)
    la, lo = np.radians(lat), np.radians(lon); N = a / np.sqrt(1 - e2 * np.sin(la) ** 2)
    return np.c_[N * np.cos(la) * np.cos(lo), N * np.cos(la) * np.sin(lo), N * (1 - e2) * np.sin(la)]

def _state(t, a, inc, raan0, raan_dot, u0, u_dot, gmst0):
    """ECEF position (km) and unit velocity direction for a circular J2-precessing orbit at times t (s since t0)."""
    u = u0 + u_dot * t; raan = raan0 + raan_dot * t; th = gmst0 + OMEGA_E * t
    cu, su, cr, sr = np.cos(u), np.sin(u), np.cos(raan), np.sin(raan); ci, si = math.cos(inc), math.sin(inc)
    x = a * (cr * cu - sr * su * ci); y = a * (sr * cu + cr * su * ci); z = a * (su * si)
    vx = -cr * su - sr * cu * ci; vy = -sr * su + cr * cu * ci; vz = cu * si
    ce, se = np.cos(th), np.sin(th)
    R = np.stack([ce * x + se * y, -se * x + ce * y, z], axis=-1)
    V = np.stack([ce * vx + se * vy, -se * vx + ce * vy, vz], axis=-1)   # Earth-rotation term ignored (geometry only)
    return R, V


def looks(stations, t0, t1, h_km=620.0, ltan_h=18.0, phase_deg=0.0, theta=(25.0, 45.0), side="right", step_s=60.0,
          inc_deg=None, return_geom=False):
    """Look times (hours since 1970, UTC) per station index for one satellite.
    side: 'right', 'left', 'alt' (left and right on alternate orbits) or 'both' (either side on every pass: an upper bound).
    Method: d.V (line of sight dotted with velocity) on a coarse time grid; each + to - sign change is a closest approach;
    its time is refined with two secant steps on the exact orbit, then the off-nadir angle and side are tested there."""
    a = RE + h_km; n = math.sqrt(MU / a ** 3)
    inc = math.radians(inc_deg if inc_deg is not None else sso_inclination(h_km))
    k = 0.75 * J2 * (RE / a) ** 2
    raan_dot = -1.5 * n * J2 * (RE / a) ** 2 * math.cos(inc)
    u_dot = n * (1 + k * (6 - 8 * math.sin(inc) ** 2))
    T0 = pd.Timestamp(t0, tz="UTC"); total = (pd.Timestamp(t1, tz="UTC") - T0).total_seconds()
    d0 = (T0 - J2000).total_seconds() / 86400.0
    sun_ra = math.radians((280.460 + 0.9856474 * d0) % 360)
    raan0 = sun_ra + math.radians((ltan_h - 12.0) * 15.0)
    gmst0 = math.radians((280.46061837 + 360.98564736629 * d0) % 360)
    u0 = math.radians(phase_deg)
    P = (a, inc, raan0, raan_dot, u0, u_dot, gmst0)
    S = stations_ecef(stations.lat.to_numpy(), stations.lon.to_numpy())
    tgrid = np.arange(0.0, total + step_s, step_s)
    ks, ss = [], []
    chunk = 4000
    for c0 in range(0, len(tgrid) - 1, chunk):
        tt = tgrid[c0:min(c0 + chunk + 1, len(tgrid))]
        R, V = _state(tt, *P)
        dot = S @ V.T - np.einsum("tk,tk->t", R, V)[None, :]           # (ns, nt): (S - R).V
        s_i, t_i = np.nonzero((dot[:, :-1] > 0) & (dot[:, 1:] <= 0))
        f = dot[s_i, t_i] / (dot[s_i, t_i] - dot[s_i, t_i + 1])
        ks.append(tt[t_i] + f * step_s); ss.append(s_i)
    tc = np.concatenate(ks); si_ = np.concatenate(ss)
    for _ in range(2):                                                  # secant refinement on the exact orbit
        R1, V1 = _state(tc, *P); R2, V2 = _state(tc + 1.0, *P)
        f1 = np.einsum("nk,nk->n", S[si_] - R1, V1); f2 = np.einsum("nk,nk->n", S[si_] - R2, V2)
        tc = tc - f1 / np.where(f2 - f1 == 0, -1e-9, f2 - f1)
    R, V = _state(tc, *P); los = S[si_] - R; rn = np.linalg.norm(R, axis=1)
    cosang = np.einsum("nk,nk->n", los, -R) / (np.linalg.norm(los, axis=1) * rn)
    ok = (cosang >= math.cos(math.radians(theta[1]))) & (cosang <= math.cos(math.radians(theta[0])))
    right = np.einsum("nk,nk->n", los, np.cross(V, R / rn[:, None])) > 0
    if side == "right": ok &= right
    elif side == "left": ok &= ~right
    elif side == "alt":
        orbit_no = np.floor((u0 + u_dot * tc) / (2 * math.pi)).astype(int)
        ok &= right == (orbit_no % 2 == 0)
    base_h = (T0 - pd.Timestamp("1970-01-01", tz="UTC")).total_seconds() / 3600.0
    if return_geom:                 # per look: time (h since 1970), off-nadir angle (deg), ascending (1) / descending (0)
        ang = np.degrees(np.arccos(np.clip(cosang, -1, 1))); asc = (V[:, 2] > 0).astype(float)
        out = [[] for _ in range(len(S))]
        for s_, t_, a_, d_ in zip(si_[ok], tc[ok], ang[ok], asc[ok]):
            out[s_].append((t_ / 3600.0 + base_h, a_, d_))
        return [np.array(sorted(o)).reshape(-1, 3) for o in out]
    out = [[] for _ in range(len(S))]
    for s_, t_ in zip(si_[ok], tc[ok]):
        out[s_].append(t_)
    return [np.sort(np.array(o)) / 3600.0 + base_h for o in out]


def theta_for_swath(h_km, swath_km, centre_deg=35.0):
    """Off-nadir range (t1, t2) whose ground width is swath_km, centred in ground range on the off-nadir angle centre_deg."""
    from scipy.optimize import brentq
    g_c = ground_range_km(h_km, centre_deg)
    g1, g2 = g_c - swath_km / 2, g_c + swath_km / 2
    inv = lambda g: brentq(lambda th: ground_range_km(h_km, th) - g, 0.01, 70.0)
    return inv(max(g1, 0.1)), inv(g2)


if __name__ == "__main__":
    for h in (500, 560, 620, 700):
        print(h, "km: SSO inclination", round(sso_inclination(h), 2), "deg | access 25-45 deg", round(access_width_km(h, 25, 45)), "km",
              "| 20-50 deg", round(access_width_km(h, 20, 50)), "km")
