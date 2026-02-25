# Методика IMERG (основная): переход к v5_year_anchor

## 1. Что было не так на старте
После первых версий калибровки оставался годовой разрыв между растром и станциями, причем:
- чем экстремальнее год по осадкам, тем больше ошибка годовой суммы.

Цель методики: одновременно улучшить среднюю точность и уменьшить «разбег» именно в экстремальных годах.

## 2. Эволюция решений
| Версия | Идея | Ограничение |
|---|---|---|
| v2 | устранение temporal mismatch + фиксация KNN + annual sanity | остаточный разрыв в экстремумах |
| v3_soft | мягкий seasonal QM + мягкие guards | экстремальные годы все еще хуже средних |
| v4_transfer | годовой transfer raw->station | частично помогает, но не полностью |
| v5_year_anchor | годовой anchor по ratio station/sat + transfer + sanity | лучшее качество в рамках текущих данных |

## 3. Основной алгоритм v5
Для каждого года:
1. 30-мин слоты агрегируются в 3-часовые суммы.
2. Применяется seasonal soft-QM по зонам ближайшей станции.
3. Выполняется пропорциональная дезагрегация обратно в слоты.
4. Применяется daily soft-guard.
5. Применяется годовой year-anchor (ratio station/sat для этого года, только для полных лет).
6. Применяется annual transfer.
7. Применяется финальный annual sanity-guard.

## 4. Обоснование выбора v5 как main
Контроль Казань, полные годы 2001-2024:
- `v3_soft`: mean/median/max `|PBIAS|` = `7.31% / 4.95% / 26.78%`
- `v5_year_anchor`: mean/median/max `|PBIAS|` = `4.12% / 3.27% / 13.78%`
- связь ошибки с экстремальностью года:
  - `corr(extremeness, |PBIAS|)`: `0.157 -> 0.037`

Вывод: v5 лучше контролирует и «обычные», и экстремальные годы.

## 5. Проверка физической правдоподобности по AOI
Для `v5` (2001-2025):
- `ann_p99`: `520.35 .. 863.21` мм/год
- `ann_max`: `579.52 .. 901.95` мм/год
- `neg_cnt=0`, `nan_cnt=0`

Относительно `v3_soft` сдвиг умеренный, без нефизичного раздува.

## 6. Режимы применения
- `v5_year_anchor` (основной): историческая реконструкция, когда station-данные этого же года доступны.
- `v4_transfer` / `v3_soft`: строгий out-of-sample режим без year-anchor.

## 7. Полный воспроизводимый пайплайн
```bash
python src/main.py --dataset imerg
python src/apply_qm_to_rasters.py --dataset imerg --zip "D:\Cache\Yandex.Disk\РНФ25-28\Осадки\IMERG_RFACTOR_ANNUAL-20260223T120530Z-1-001.zip" --calib output/calib_imerg
python src/build_imerg_evidence_pack.py
```

## 8. Где смотреть результаты
- `docs/evidence_imerg_qc_main.md`
- `docs/kazan_summary_v2_v3_v4_v5.csv`
- `docs/kazan_v2_v3_v4_v5_full_years.csv`
- `docs/aoi_annual_quantiles_v5_year_anchor.csv`
- `docs/aoi_v3_vs_v5_delta.csv`

