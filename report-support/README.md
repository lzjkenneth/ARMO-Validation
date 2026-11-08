# Report support

Scripts behind the new numbers and figures in the final ARMO report (8 November 2026), so they can be re-run from this
repository. Run them from the repository root with the packages in `../requirements.txt`. In a clean clone they
reproduce the report's numbers and figure files exactly.

| Script | What it reproduces | Inputs | Output |
|---|---|---|---|
| `launch_mass_budget.py` | Mass, delta-v, launch and budget numbers: 221 m/s of delta-v (raise from a 500 km rideshare insertion to the 620 km orbit with the 0.46° sun-synchronous inclination trim, drag make-up, deorbit to a 300 km perigee, plus 10 %), 331 kg wet mass per satellite (299 kg dry), the rideshare launch cost at SpaceX's published per-kg price, the budget lines and reserves, the one-satellite threshold mission and the cost of commercial X-band tasking. Set `H = 613.0` for the exact-repeat alternative (213 m/s, 330 kg) | None: constants and their sources are in the file | Printed |
| `gap_by_orbit.py` | Station-by-station mean revisit gap of design D (two satellites 180° apart in one dawn–dusk plane) at 620 km, at five exact-repeat altitudes (599.9–623.5 km) and with two wider access ranges, averaged over 12 starting positions; also capture alone and with today's fleet | `../data/processed/` (events, station list, Sentinel-1 and NISAR look times) | `data/gap_by_orbit.csv`, `data/gap_by_orbit_stations.csv` (committed) |
| `fig_gap_latitude.py` | Station mean gap against latitude for the 620 km orbit, the 613 km exact repeat and the 24–46° access range, with the 1.5-day line | `data/gap_by_orbit_stations.csv`, `../data/processed/gap_by_station.csv` | `figures/val_gap_by_latitude.png` |
| `fig04_workflow_compact.py` | Report Figure 4: the analysis workflow (phases, objectives, validation steps) | None | `figures/fig04_workflow_compact.png` |
| `fig_print_redraws.py` | Report Figures 6, 9 and 10 (detection progression, capture gains, usable-strip sensitivity), drawn at their printed widths from the team's own output tables | `data/team_generated/` (four tables from [zheniannn/ARMO](https://github.com/zheniannn/ARMO), MIT licence) | `figures/fig06_print.png`, `fig09_print.png`, `fig10_print.png` |
| `fig08_design_d.py` | Report Figure 8: design D on 22 August 2026 over the contiguous United States (access swaths, ground tracks, scale bar, north arrow) and the Stillwater, Oklahoma timeline of 5 cm soil moisture, wetting windows and Sentinel-1, NISAR and ARMO looks | `../data/processed/` (events, look times and `uscrn_hourly.parquet`, which `code/01_ground_truth.py` writes) and the Natural Earth states file (below) | `figures/fig08_design_d.png` |

From `data/gap_by_orbit.csv`: at 620 km with 25–45° off-nadir access the mean gap over the 107 stations is 1.42 days
(1.35 north of 35° N, 1.54 south), the worst station 1.66 days averaged over 12 starting positions (1.64 days in
`../data/processed/gap_by_station.csv`, which the report quotes), and with today's fleet no station averages more
than 1.22 days. The 613 km exact repeat has a 1.47-day mean but leaves four Gulf-coast stations at 2.0–2.3 days.
Widening the access to 24–46° gives a 1.28-day mean, a 1.49-day worst station and 79.3 % capture for ARMO alone.

The dry mass, growth margin and propellant specific impulse in `launch_mass_budget.py` come from the team's sizing in
[zheniannn/ARMO](https://github.com/zheniannn/ARMO) (`code/06_mission_design.py`); the launch price comes from SpaceX's
published rideshare price list (read 8 November 2026).

## Run

```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
curl -L -o report-support/data/ne_50m_admin_1_states_provinces.geojson \
  https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_1_states_provinces.geojson
.venv/bin/python code/01_ground_truth.py           # only needed for fig08: USCRN hourly files (about 600 MB)
.venv/bin/python report-support/launch_mass_budget.py
.venv/bin/python report-support/gap_by_orbit.py     # under a minute
.venv/bin/python report-support/fig_gap_latitude.py
.venv/bin/python report-support/fig04_workflow_compact.py
.venv/bin/python report-support/fig_print_redraws.py
.venv/bin/python report-support/fig08_design_d.py
```

The Natural Earth file (public domain, 2.3 MB) and the generated figures are not committed (see `.gitignore`).
