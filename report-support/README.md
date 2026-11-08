# Report support (8 November 2026)

Scripts behind new numbers and one figure in the revised ARMO report of 8 November 2026. They were written during the
report revision and copied here so the numbers can be re-run from this repository. Run them from the repository root
with the packages in `../requirements.txt`. Checked on 8 November 2026: in a fresh clone these copies give
byte-identical CSV files and figures to the ones used for the report.

| Script | What it reproduces | Inputs | Output |
|---|---|---|---|
| `launch_mass_budget.py` | Mass, delta-v, launch and budget numbers: 213 m/s of delta-v (raise from a 500 km insertion to 613 km with the 0.44° sun-synchronous trim, drag make-up, deorbit to a 300 km perigee, plus 10 %), 330 kg wet mass per satellite (299 kg dry), the rideshare launch cost at SpaceX's published per-kg price, the budget lines and reserves, a one-satellite threshold mission and the cost of commercial X-band tasking | None: constants and their sources are in the file | Printed |
| `gap_by_orbit.py` | Station-by-station mean revisit gap of design D (two satellites 180° apart in one dawn-dusk plane) at 620 km, at five exact-repeat altitudes (599.9-623.5 km) and with two wider access ranges, averaged over 12 starting positions; also capture alone and with today's fleet | The committed files in `../data/processed/` (events, station list, Sentinel-1 and NISAR look times) | `data/gap_by_orbit.csv`, `data/gap_by_orbit_stations.csv` (committed) |
| `fig_gap_latitude.py` | Figure: station mean gap against latitude for the proposed 620 km orbit, the 613 km exact-repeat orbit and the 24-46° access range, with the 1.5-day Level-1 line | `data/gap_by_orbit_stations.csv`, `../data/processed/gap_by_station.csv` | `figures/val_gap_by_latitude.png` |
| `fig08_design_d.py` | Report Figure 8: design D on 22 August 2026 over the contiguous United States (access swaths, ground tracks, scale bar, north arrow) and the Stillwater, Oklahoma timeline of 5 cm soil moisture, wetting windows and Sentinel-1, NISAR and ARMO looks | `../data/processed/` (events, look times and `uscrn_hourly.parquet`, which `code/01_ground_truth.py` writes and which is not committed) and the Natural Earth states file (below) | `figures/fig08_design_d.png` |

From `data/gap_by_orbit.csv`: at the proposed 620 km with 25-45° off-nadir access, the mean gap over the 107 stations
is 1.42 days (1.54 south of 35° N, 1.35 north of it), the worst station 1.66 days, and 29 stations are above 1.52
days. The 613 km exact-repeat orbit has a 1.47-day mean but a 2.34-day worst station. Widening access to 24-46° gives
a 1.28-day mean, a 1.49-day worst station and 79.3 % capture for ARMO alone.

## Run

```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
curl -L -o report-support/data/ne_50m_admin_1_states_provinces.geojson \
  https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_1_states_provinces.geojson
.venv/bin/python code/01_ground_truth.py           # only needed for fig08: USCRN hourly files (about 600 MB)
.venv/bin/python report-support/launch_mass_budget.py
.venv/bin/python report-support/gap_by_orbit.py     # under a minute
.venv/bin/python report-support/fig_gap_latitude.py
.venv/bin/python report-support/fig08_design_d.py
```

The Natural Earth file (public domain, 2.3 MB) and the generated figures are not committed (see `.gitignore`).

