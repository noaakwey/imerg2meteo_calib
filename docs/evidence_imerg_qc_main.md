# IMERG: доказательная база эволюции v2 -> v3 -> v4 -> v5 (основной отчет)

## 1. Исходная проблема
- В годовых калиброванных растрах сохранялся систематический разрыв со станциями.
- В экстремальные годы (очень влажные/сухие) ошибка росла сильнее всего.
- Цель: уменьшить годовой разрыв и ослабить зависимость ошибки от экстремальности года, не выходя за физически правдоподобные пределы по AOI.

## 2. Эксперименты
| Версия | Что изменили | Ожидаемый эффект |
|---|---|---|
| v2 | фиксация mismatch по срокам + фиксация KNN-координат + annual sanity envelope | убрать грубое завышение |
| v3_soft | адаптивный soft seasonal QM + мягкие daily/annual guard | снизить перекоррекцию |
| v4_transfer | annual raw->station transfer | улучшить межгодовую амплитуду |
| v5_year_anchor | годовой anchor по station/sat ratio + transfer + sanity | целенаправленно сжать разрыв в экстремумах |

## 3. Сравнение по Казани (полные годы 2001-2024, N=24)
- `v2`: mean/median/max `|PBIAS|` = `9.48% / 7.46% / 32.15%`; `corr(ext, |PBIAS|)=0.233`
- `v3_soft`: mean/median/max `|PBIAS|` = `7.31% / 4.95% / 26.78%`; `corr(ext, |PBIAS|)=0.157`
- `v4_transfer`: mean/median/max `|PBIAS|` = `7.05% / 4.93% / 26.78%`; `corr(ext, |PBIAS|)=0.105`
- `v5_year_anchor`: mean/median/max `|PBIAS|` = `4.12% / 3.27% / 13.78%`; `corr(ext, |PBIAS|)=0.037`

## 4. Почему v5 выбран основным
- Минимальный средний `|PBIAS|` среди всех версий (`4.12%`).
- Минимальная ошибка в худшем году (`13.78%` против `26.78%` в v3).
- Самая слабая связь ошибки с экстремальностью (`0.037`), т.е. разрыв в экстремумах реально сжат.

## 5. Sanity-check по AOI (v5, 2001-2025)
- `ann_p99`: `520.35 .. 863.21` мм/год
- `ann_max`: `579.52 .. 901.95` мм/год
- отрицательных пикселей: `0`
- NaN пикселей: `0`
- относительно v3_soft: средний сдвиг `p99 = +1.97%`, `max = +1.82%` (умеренно, без нефизичного раздува)

## 6. Station-level годовой баланс (точечные CSV, 2001-2024)
- Median station `|annual PBIAS|`: `9.79% -> 8.45%`
- Mean station `|annual PBIAS|`: `24.44% -> 18.42%`

## 7. Артефакты
- Основные таблицы:
  - `docs/kazan_v2_v3_v4_v5_full_years.csv`
  - `docs/kazan_summary_v2_v3_v4_v5.csv`
  - `docs/aoi_annual_quantiles_v5_year_anchor.csv`
  - `docs/aoi_v3_vs_v5_delta.csv`
  - `docs/station_annual_pbias_summary_imerg.csv`
- Основные графики:
  - `docs/graphics/evidence_kazan_v3_vs_v5_series.png`
  - `docs/graphics/evidence_kazan_abs_pbias_v3_vs_v5.png`
  - `docs/graphics/evidence_extremeness_vs_abs_pbias_v3_vs_v5.png`
  - `docs/graphics/evidence_kazan_abs_pbias_box_v2_v3_v4_v5.png`
  - `docs/graphics/evidence_experiment_progression_mean_abs_pbias.png`
  - `docs/graphics/evidence_aoi_envelope_v3_vs_v5.png`
- Архивные материалы (трассируемость):
  - `docs/evidence_imerg_qc.md`
  - `docs/evidence_imerg_qc_v3_soft.md`
  - `docs/evidence_imerg_qc_v5_year_anchor.md`
  - `docs/methodology_imerg_v2.md`

## 8. Ограничение применения
- `v5_year_anchor` — лучший режим для исторической реконструкции, когда доступны годовые station-данные того же года.
- Для строгого out-of-sample режима (без годового якоря того же года) использовать `v4_transfer` или `v3_soft`.

