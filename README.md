# IMERG/ERA5 precipitation calibration for RUSLE R-factor

## Purpose
Repository for station-based calibration of precipitation time series and raster stacks, focused on:
- realistic annual water balance,
- preserved sub-daily intensity structure for `R` and `I30`,
- robust behavior in extreme wet/dry years.

Primary operational IMERG mode is now **`v5_year_anchor`**.

## Current main method (IMERG)
Main raster workflow (`src/apply_qm_to_rasters.py`, dataset `imerg`) applies:
1. Seasonal 3-hour soft QM (station-zone based, nearest-station Voronoi).
2. Daily soft guard.
3. Year-specific annual anchor (`station_year / sat_year` ratio for full years).
4. Annual transfer mapping.
5. Final annual sanity envelope guard.

This is the best validated setup for historical reconstruction in this project.

## Quick start
```bash
pip install -r requirements.txt
```

### 1) Build station calibration tables (`*_calib.csv`)
```bash
python src/main.py --dataset imerg
python src/main.py --dataset era5land
```

Outputs:
- `output/calib_imerg/`
- `output/calib_era5land/`

### 2) Apply IMERG calibration to annual raster archive (main pipeline)
```bash
python src/apply_qm_to_rasters.py \
  --dataset imerg \
  --zip "D:\Cache\Yandex.Disk\РНФ25-28\Осадки\IMERG_RFACTOR_ANNUAL-20260223T120530Z-1-001.zip" \
  --calib output/calib_imerg
```

Default IMERG output (main):
- `output/imerg_rfactor_calib_v5_year_anchor/`

### 3) Build full evidence pack (tables + plots + report)
```bash
python src/build_imerg_evidence_pack.py
```

Main evidence report:
- `docs/evidence_imerg_qc_main.md`

## Full pipeline (recommended order)
1. Prepare/update source station and satellite tables in `data/`.
2. Run `src/main.py --dataset imerg`.
3. Run `src/apply_qm_to_rasters.py --dataset imerg ...`.
4. Run `src/build_imerg_evidence_pack.py`.
5. Review:
   - `docs/kazan_summary_v2_v3_v4_v5.csv`
   - `docs/kazan_v2_v3_v4_v5_full_years.csv`
   - `docs/aoi_annual_quantiles_v5_year_anchor.csv`
   - `docs/evidence_imerg_qc_main.md`

## Historical experiments (for traceability)
- `v2`: temporal mismatch + KNN fix + annual sanity.
- `v3_soft`: adaptive soft-QM, softer daily/annual guards.
- `v4_transfer`: annual raw-to-station transfer.
- `v5_year_anchor` (main): year-specific annual anchor + transfer + sanity.

Comparison tables/plots are generated automatically by `src/build_imerg_evidence_pack.py`.

## Modes and when to use
- `v5_year_anchor` (main): best for **historical reconstruction** where same-year station annual totals are available.
- `v4_transfer` / `v3_soft`: safer for **strict out-of-sample** scenarios (no same-year station anchor).

## Important paths
- Station calibration metrics:
  - `output/calib_imerg/validation_metrics_imerg.csv`
- Main IMERG rasters:
  - `output/imerg_rfactor_calib_v5_year_anchor/`
- Evidence docs:
  - `docs/README.md`
  - `docs/evidence_imerg_qc_main.md`
  - `docs/methodology_imerg_main.md`

## Script map
- `src/main.py`: builds per-station calibrated tables and validation metrics.
- `src/apply_qm_to_rasters.py`: applies calibration to raster archives.
- `src/build_imerg_evidence_pack.py`: regenerates comparison tables, graphics, and final evidence report.
- `src/analyze_results.py`: boxplots/scatter for station-level validation metrics.
