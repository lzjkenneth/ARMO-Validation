# Radar data access

Sentinel-1 (C-band) and NISAR (L-band) look times, backscatter and the NISAR soil-moisture product are all open, and
all are used here. The other L-band missions (ALOS-2, SAOCOM, ALOS-4) are not practical comparisons for US stations in
2025–2026, and ROSE-L has not launched. Everything below was checked on 7 Nov 2026.

## Mission by mission

| Mission (band) | What we can get | How | Used here | Coverage of the 121 USCRN stations |
|---|---|---|---|---|
| Sentinel-1 A/C/D (C-band) | Every acquisition time, GRD and RTC backscatter (10 m) | Microsoft Planetary Computer STAC, no login | Look times (`03_fleet_looks.py`) | 14,151 looks, Mar 2025–Oct 2026. In 17 Jun–31 Oct 2026 most US looks came from Sentinel-1D (1,864) and Sentinel-1C (996); Sentinel-1A only 189 |
| NISAR (L-band, 1.239 GHz, 40 MHz) | Acquisition times; GCOV backscatter (10 m, HH + HV); SME2 soil moisture (200 m) | NASA CMR search (no login); ASF download with a free Earthdata login | Look times; 10 m backscatter windows at 16 stations (`07_nisar_windows.py`) | 6,313 looks Oct 2025–Oct 2026, but a hole from February to mid-June 2026 (below) |
| ALOS-2 PALSAR-2 (L-band) | ScanSAR Level 1 held at ASF (`ALOS2_L1_PSR2`) | NASA CMR / ASF | Look times only | 65 looks at 18 stations, all in March–May 2025, none after. Not usable for 2026 |
| SAOCOM-1A/1B (L-band) | Only via ESA's third-party mission scheme or CONAE, by project application | Application, not open download | Not used | Not checked per station (no open catalogue) |
| ALOS-4 PALSAR-3 (L-band) | Not in NASA's catalogue; JAXA distribution | Not open | Not used | — |
| ROSE-L (L-band) | Not launched | — | — | — |
| SMAP (L-band radiometer) | 9 km and 36 km soil moisture | NSIDC with Earthdata login | Not re-run here (the team's repo covers it) | Daily-ish, but 9 km |

## NISAR's gaps, measured from the catalogue

All from real catalogue records, not from NISAR's design numbers.

1. **A four-month hole.** Products for most US stations stop in February 2026 and restart on 17 June 2026. From
   February to May only 46 of 121 stations had any NISAR look (824 looks in total, mostly a few West Coast sites).
   The only continuous record is 17 June–31 October 2026, which is why the matched cohort starts on 17 June.
2. **Slow revisit, and uneven.** In 17 June–31 October each station got 18–40 NISAR looks (median 28), a mean gap
   of 4.9 days at the median station. Sentinel-1 is slower still over the US (median 6.5 days).
3. **Same two times of day.** Every NISAR look falls at about 03:00–06:00 or 18:00–21:00 local solar time, and
   Sentinel-1 is also near 06:00/18:00. More dawn-dusk satellites cannot add new times of day. Our sweep shows that
   time of day hardly matters for capture, though (see RESULTS.md section 7).
4. **Late delivery.** Over four sample stations (Jun–Oct 2026) the provisional GCOV backscatter reached the
   catalogue a median 74 hours after the pass (90th percentile about 14 days). The SME2 soil-moisture product took a
   median 105 hours (90th percentile about 19 days). The median wetting window is 1.8 days, so even when NISAR looks
   inside a window, its data usually arrive after the window has closed. These are provisional products in NISAR's first
   year, so the latency may improve. `data/processed/nisar_latency_sample.csv`
5. **HH-led soil moisture at 200 m.** The SME2 product is 200 m, and its quality flags mark many retrievals "not
   recommended" in dense vegetation (team repo, METHOD_NOTES section 5).

## Summary

- Over the US in 2026 the openly available field-scale radar records are Sentinel-1 and NISAR; SAOCOM and ALOS-2/4
  data are not openly available for these stations.
- NISAR's provisional products reach the catalogue a median 3–4 days after the pass, longer than the median wetting
  window (1.8 days), which is why a delivery target of 24 h matters for this application.
- NISAR's nominal 12-day repeat corresponds to a measured 4.9-day mean gap per station (ascending, descending and
  overlapping tracks), with a four-month gap in 2026 products and two fixed times of day.
