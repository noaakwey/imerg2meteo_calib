# IMERG: снижение разрыва в экстремальные годы (v5_year_anchor)

## Что изменено по методике
- Базовая цепочка `v3_soft` сохранена:
  - seasonal 3h soft-QM,
  - мягкий суточный guard,
  - annual sanity guard.
- Добавлен новый слой `v5` (кардинально иной относительно `v3`):
  - годовой **year-anchor** по станции для конкретного года:
    - из `*_calib.csv` берется годовой `ratio = P_station_year / P_sat_year`,
    - якорь применяется только для полноценных лет (>=99% 3-часовых сроков),
    - применяется мягко, но существенно, перед annual-transfer.
  - после year-anchor остается annual-transfer + финальный sanity guard.

## Результат по Казани (WMO 27595, полные годы 2001-2024)
- Сводка `|PBIAS|`:
  - `v3_soft`: mean/median/max = `7.31% / 4.95% / 26.78%`
  - `v5_year_anchor`: mean/median/max = `4.12% / 3.27% / 13.78%`
- Систематический сдвиг (PBIAS):
  - mean: `-4.24% -> -0.51%`
  - median: `-4.40% -> -1.56%`
- Экстремальные годы (по `ext_idx >= 1.0`):
  - mean `|PBIAS|`: `9.99% -> 4.99%`
- Зависимость ошибки от экстремальности года:
  - `corr(extremeness, |PBIAS|)`: `0.157 -> 0.037`
- Улучшение по `|PBIAS|`:
  - `v5` лучше `v3` в `17/24` лет.
- Наихудший год:
  - 2001: `|PBIAS| 26.78% -> 13.78%`.

## AOI sanity-check (все годы 2001-2025)
- `v5` envelope:
  - `ann_p99`: `520.35 .. 863.21` мм/год
  - `ann_max`: `579.52 .. 901.95` мм/год
  - отрицательных/NaN пикселей: `0/0`
- Относительно `v3`:
  - средний сдвиг `ann_p99`: `+1.97%`
  - средний сдвиг `ann_max`: `+1.82%`
- Вывод: разрыв в экстремумах снижен без неадекватного раздува годовых сумм по AOI.

## Артефакты
- Растры:
  - `output/imerg_rfactor_calib_v5_year_anchor/`
- Таблицы:
  - `docs/kazan_v2_v3_v4_v5_full_years.csv`
  - `docs/kazan_summary_v2_v3_v4_v5.csv`
  - `docs/aoi_annual_quantiles_v5_year_anchor.csv`
  - `docs/aoi_v3_vs_v5_delta.csv`
- Графики:
  - `docs/graphics/evidence_kazan_v3_vs_v5_series.png`
  - `docs/graphics/evidence_kazan_abs_pbias_v3_vs_v5.png`
  - `docs/graphics/evidence_extremeness_vs_abs_pbias_v3_vs_v5.png`
  - `docs/graphics/evidence_kazan_abs_pbias_box_v2_v3_v4_v5.png`
  - `docs/graphics/evidence_aoi_envelope_v3_vs_v5.png`

## Ограничения и риски
- `v5` целенаправленно использует годовую station-информацию в том же году (историческая реконструкция, не прогноз).
- Для сценариев вне периода наблюдений (future/pure holdout) нужно переключаться на `v4`/`v3` (без year-anchor) или явно фиксировать режим применения.
