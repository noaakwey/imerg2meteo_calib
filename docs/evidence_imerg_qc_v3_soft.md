# IMERG: доказательная база качества калибровки (v3_soft)

## Что изменено
- В `src/apply_qm_to_rasters.py` внедрена мягкая (soft) схема применения QM к растру:
  - station-specific `blend_alpha` (`raw + alpha * (qm - raw)`),
  - мягкий суточный guard (ограниченное масштабирование),
  - мягкий годовой sanity-guard с ограничением по обучающему envelope.
- Пересчитаны годовые растра IMERG за 2001-2025:
  - `output/imerg_rfactor_calib_v3_soft/IMERG_V07_P30min_mmh_YYYY_calib_qm.tif`

## Контроль 1: Казань (WMO 27595), годовые суммы
- Полные годы для сравнения: `2001-2024` (`N=24`, 2025 исключен как неполный).
- Сводка PBIAS (растр vs станция):
  - `v2`: mean/median `-8.43% / -7.44%`
  - `v3_soft`: mean/median `-4.24% / -4.40%`
- Сводка `|PBIAS|`:
  - `v2`: mean/median `9.48% / 7.46%`
  - `v3_soft`: mean/median `7.31% / 4.95%`
- Худший `|PBIAS|`:
  - `v2`: `32.15%`
  - `v3_soft`: `26.78%`
- Улучшение по `|PBIAS|`: `18/24` полных лет.

## Контроль 2: AOI envelope (годовые суммы по пикселям)
- По всем годам `2001-2025` для `v3_soft`:
  - `ann_p99`: `510.84 .. 861.99` мм/год
  - `ann_max`: `556.60 .. 929.08` мм/год
  - отрицательных пикселей: `0`
  - NaN пикселей: `0`
- Относительно `v2` (по годам):
  - средний сдвиг `ann_p99`: `+4.89%`
  - средний сдвиг `ann_max`: `+3.71%`
- Интерпретация: поле стало немного «влажнее», но осталось в физически правдоподобном диапазоне для региона (нет кратного раздува «как в тропиках»).

## Контроль 3: station-level метрики на калибровочных CSV
- Файл: `output/calib_imerg/validation_metrics_imerg.csv` (`202` станции).
- Daily KGE mean: `0.3409 -> 0.3424`
- Monthly KGE mean: `0.4797 -> 0.5104`
- Median `|daily PBIAS|`: `8.59% -> 6.92%`
- Улучшено станций:
  - daily KGE: `47/202`
  - monthly KGE: `81/202`
  - `|daily PBIAS|`: `102/202`

## Контроль 4: station-level годовой баланс (2001-2024)
- Агрегировано по `200` станциям (`docs/station_annual_pbias_summary_v3_soft.csv`):
  - median `|annual PBIAS|` по станциям: `9.79% -> 8.45%`
  - mean `|annual PBIAS|` по станциям: `24.44% -> 18.42%`
  - станций с median `|annual PBIAS| <= 10%`: `103 -> 126`
  - станций с mean annual overbias `>20%`: `37 -> 27`

## Артефакты проверки
- Таблицы:
  - `docs/kazan_annual_balance_v3_soft.csv`
  - `docs/kazan_v2_vs_v3_soft_full_years.csv`
  - `docs/aoi_annual_quantiles_v3_soft.csv`
  - `docs/aoi_v2_vs_v3_soft_delta.csv`
  - `docs/station_annual_pbias_summary_v3_soft.csv`
- Графики:
  - `docs/graphics/evidence_kazan_v2_vs_v3_soft_series.png`
  - `docs/graphics/evidence_kazan_v2_vs_v3_soft_pbias.png`
  - `docs/graphics/evidence_aoi_annual_envelope_v3_soft.png`
  - `docs/graphics/evidence_station_annual_abs_pbias_v3_soft.png`
