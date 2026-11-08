# ARMO validation: a companion analysis for the ARMO mission proposal

**ARMO** (Agricultural Rapid Moisture Observatory) is a two-satellite L-band SAR mission proposed by ES4304 Group 2
(Satellite Remote Sensing, NTU, 2026) to observe soil wetting within days of rain. This repository is a companion to the
team's ARMO analysis ([zheniannn/ARMO](https://github.com/zheniannn/ARMO), built by Tang Zhen Ian). It implements the
capture analysis a second time from the definitions in the team's report, as a cross-check, and extends it with
further tests: waiting time from wetting to the first look, latitude, design sweeps, orbit repeat options and measured
L-band noise. It also holds the scripts behind the final report's engineering, budget and revisit numbers
([report-support/](report-support/README.md)). Everything uses open data: NOAA USCRN hourly ground truth, real
Sentinel-1 and NISAR acquisition records, and NISAR and Sentinel-1 images.

The cross-check was made on 7 November 2026, after the group presented ARMO on 5 November, and the report-support
scripts were added on 8 November. How the group arrived at ARMO, from dozens of candidate problems and pilot analyses,
is summarised in [docs/IDEA_SEARCH.md](docs/IDEA_SEARCH.md).

## What it checks

The report's first three objectives, on the same data and definitions. O1 uses all 2025–2026 events; O2 and O3 use the
matched cohort of storms from 17 June to 31 October 2026 (the continuous NISAR record): 982 events at 107 stations
here, 979 in the team's run.

| Report objective | This cross-check | Team |
|---|---|---|
| O1: wetting windows are short | 3,343 events at 120 stations; median window 1.79 days; 55 % shorter than 48 h | 3,340 events; 1.8 days; 55 % |
| O2: today's satellites miss about half of the events | Sentinel-1 30.3 %, NISAR 37.4 %, either 51.3 % | 30.2, 37.9 and 51.6 % |
| O3: how often a satellite must look, and what ARMO adds | 80 % capture needs a look every 1.45 days; ARMO alone 77.0 %, today + ARMO 82.2 % | 1.46 days; 76.1 % and 81.5 % |

The headline numbers reproduce within about one percentage point. Details, confidence intervals and the method
differences (for example, gap filling of up to 3 hours) are in [docs/RESULTS.md](docs/RESULTS.md), sections 1, 2 and 6.

## Main results

| | Today (Sentinel-1 + NISAR) | ARMO alone (2 satellites) | Today + ARMO |
|---|---|---|---|
| Wetting events with a look inside the window | 51.3 % (95 % CI 47.3–55.6) | 77.0 % (73.7–80.3) | 82.2 % (78.9–85.4) |
| Median wait from wetting to first look | 39 h | 19 h | 13 h |
| Events looked at within 24 h | 35 % | 60 % | 74 % |
| Mean gap between looks | 2.81 days | 1.42 days | 0.93 days |
| Worst gap at any station | 20.0 days | 4.5 days | 4.5 days |

982 wetting events at 107 USCRN stations, 17 June–31 October 2026. Today's fleet uses real look times; ARMO is
simulated (620 km dawn–dusk orbit, 25–45° off-nadir, 8 starting positions).

![Capture and waiting time](figures/fig_capture_and_wait.png)

- **The team's numbers reproduce** within one point (team: 51.6 %, 76.1 %, 81.5 %, gain +29.9 points; here +30.9,
  95 % CI 28.0–33.6).
- **Two satellites are needed.** Event by event, two satellites beat the best single satellite by 8.4 points alone
  (95 % CI 7.1–9.7). One satellite matches two only by imaging about twice the ground per pass.
- **The gain holds** under every wetting-window rule tried (+24 to +31 points) and without NISAR (Sentinel-1 + ARMO
  79.4 %).
- **The revisit gap varies with latitude:** 1.42 days on average, 1.54 in the south and 1.27 in the north; with
  today's fleet every station is below 1.4 days.
- **Orbit repeat cycle matters.** Altitudes near short repeat cycles capture far less (560 km: 60.3 %; 720 km: 48.3 %).
  An exact 6-day repeat at 613 km keeps capture and gives a same-geometry image every 3 days instead of every 45.
- **Noise.** Measured NISAR pairs show change noise above about 50 m is real scene change rather than speckle (about
  1.6 dB at the station from 100 to 500 m); ARMO's burst-mode speckle is 0.95 dB at 200 m and 0.38 dB at 500 m.

## Documents

| File | What it holds |
|---|---|
| [docs/RESULTS.md](docs/RESULTS.md) | All results with tables and figures: capture, waiting time, latitude, event definition, design sweeps, orbit repeat options, measured noise, scope and assumptions |
| [docs/DATA_ACCESS.md](docs/DATA_ACCESS.md) | Which radar records are openly available, and NISAR's measured gaps (2026 hole, revisit, delivery delay) |
| [docs/IDEA_SEARCH.md](docs/IDEA_SEARCH.md) | How the group arrived at ARMO: concept notes, idea bank, ground-truth checks and pilots (October–November 2026) |
| [report-support/](report-support/README.md) | Scripts behind the final report's mass, delta-v, launch, budget and station-by-station revisit numbers, and Figures 4, 6, 8, 9 and 10 |

## Run it

Python 3.11+ with `numpy pandas pyarrow scipy matplotlib pyproj h5py rasterio` (`python3 -m venv .venv` then
`.venv/bin/pip install -r requirements.txt`). Steps 1–6 and 9–11 need no login.
Step 7 needs a free NASA Earthdata login in a netrc file (`EARTHDATA_NETRC`, default `~/.netrc`); it is never
printed or stored in the repository.

```
python code/01_ground_truth.py         # USCRN hourly files 2025-2026 (about 600 MB) -> quality-controlled table
python code/02_events.py               # wetting events under six window rules, dry-down fits
python code/03_fleet_looks.py          # real Sentinel-1, NISAR and ALOS-2 look times at every station
python code/05_capture.py              # capture, waiting time and gaps for today, ARMO designs A-E and both
python code/06_event_definition.py     # the 50 % rule checked three ways
python code/09_design_sweeps.py        # satellites, swath, altitude, time of day, look side, one-satellite options
python code/09b_altitude_resonance.py  # altitude in 5 km steps
python code/07_nisar_windows.py        # NISAR GCOV windows at 16 stations (range reads, about 5 GB transferred)
python code/08_noise_vs_resolution.py  # change noise against resolution, measured
python code/12_s1_pair_noise.py        # Sentinel-1 RTC windows; then code/12b_s1_noise_by_interval.py
python code/13_review_checks.py        # paired and sensitivity numbers; code/13b_usable_strip.py
python code/10_figures.py && python code/11_headline.py && python -m pytest -q tests
```

`code/04_orbit_sim.py` is the orbit model (J2 sun-synchronous, zero-Doppler looks); `code/sar_tradeoff.py` makes the
radar trade tables and `figures/feas_*.png`. Every number quoted in the documents is collected in
`data/processed/headline.json`, and `tests/test_numbers.py` checks the documents against it. Raw downloads go to
`data/raw/` and are not committed. The report-support scripts run the same way; see
[report-support/README.md](report-support/README.md).

## Data sources

| Data | Source | Access |
|---|---|---|
| USCRN hourly02 soil moisture and rain | NOAA NCEI, https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/ | Open |
| Sentinel-1 look times, RTC backscatter | Microsoft Planetary Computer (`sentinel-1-grd`, `sentinel-1-rtc`), CC BY 4.0 | Open |
| NISAR look times; GCOV backscatter | NASA CMR; ASF (`NISAR_L2_GCOV_PROVISIONAL_V1`) | CMR open; ASF free Earthdata login |
| ALOS-2 ScanSAR look times | NASA CMR (`ALOS2_L1_PSR2`) | Open metadata |

## Related repository

[zheniannn/ARMO](https://github.com/zheniannn/ARMO): the team's ARMO analysis, figures and method notes (Tang Zhen Ian).

## Licence

Code in this repository is released under the MIT License ([LICENSE](LICENSE)). Data downloaded or derived from NOAA
USCRN, Copernicus Sentinel-1 (via Microsoft Planetary Computer), NASA NISAR and ALOS-2 metadata keep their providers'
terms (see Data sources). The four tables in `report-support/data/team_generated/` come from the team's MIT-licensed
repository [zheniannn/ARMO](https://github.com/zheniannn/ARMO).

## Evidence labels

**OBSERVED** (ground truth, real acquisitions), **SIMULATED** (orbit model), **MEASURED** (real radar images),
**ASSUMED** (stated assumption). Main limits: one season (June–October 2026), few cropland stations, capture is not
detection, and the ARMO runs assume the whole swath is imaged on every pass. See `docs/RESULTS.md` section 10.
