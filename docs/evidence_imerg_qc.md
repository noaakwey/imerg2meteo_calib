# IMERG: доказательная база качества калибровки (v2)

## Контекст
- Цель: подтвердить, что новая схема калибровки убирает завышение годовых сумм и держит физически правдоподобный баланс.
- Контрольная станция: Казань (WMO 27595).

## Ключевые численные результаты
- Legacy-сетка сроков в station-файле: `[0, 3, 6, 12, 15]`.
- Среднегодовая IMERG сумма (Казань, full 8 сроков): `534.0` мм/год.
- Среднегодовая IMERG сумма (Казань, legacy-сроки): `328.3` мм/год.
- Коэффициент mismatch full/legacy: `1.63x`.
- Полные годы в сравнении растр/станция: `2001-2024` (`N=24`).
- PBIAS old raster (mean/median): `79.15% / 73.12%`.
- PBIAS new raster (mean/median): `-8.43% / -7.44%`.
- Диапазон PBIAS new raster: `-32.15% .. 7.02%`.
- KNN после фикса координат: уникальные станции `1 -> 8`; медианная дистанция `187.73 -> 64.48` км.
- AOI p99 годовых сумм (v2): `462.7 .. 842.1` мм/год; AOI max: `519.3 .. 906.4` мм/год.

## Баланс 2001 по Казани
- Станция: `691.7` мм
- Сырой IMERG (zip): `498.8` мм
- Старый калиброванный растр: `854.5` мм
- Новый калиброванный растр: `469.3` мм

## Snapshot station-level метрик (output/calib_imerg)
- Median daily KGE raw/corr: `0.558` / `0.298`
- Median monthly KGE raw/corr: `0.726` / `0.521`
- Median daily PBIAS raw/corr: `1.2%` / `15.0%`

## Бюджет ошибок и неопределенности
- Темпоральная репрезентативность: station-CSV хранят только «влажные» сроки; пропуски нужно восстанавливать как нули на полной 3ч-сетке.
- Пространственная репрезентативность: nearest-station назначение вносит локальный шум, особенно в разреженной сети.
- Хвосты распределения: сезонный QM может раздувать объем без суточного/годового guard.
- Нестационарность: перенос, обученный на 2001-2015, требует пост-валидации для поздних лет.

## Графики
- [Комбинации сроков по станциям](graphics/evidence_station_hour_combos.png)
- [Mismatch по Казани](graphics/evidence_kazan_annual_series.png)
- [KNN fix диагностика](graphics/evidence_knn_assignment_fix.png)
- [Баланс 2001 по Казани](graphics/evidence_kazan_2001_balance.png)
- [Envelope отношений (Kazan)](graphics/evidence_kazan_ratio_envelope.png)
- [Snapshot метрик](graphics/evidence_validation_snapshot.png)
- [Станция vs old/new растр](graphics/evidence_kazan_old_vs_new_series.png)
- [Сравнение годового PBIAS](graphics/evidence_kazan_pbias_compare.png)
- [AOI envelope годовых сумм](graphics/evidence_aoi_annual_envelope_v2.png)

## Таблицы
- `docs/kazan_annual_balance_v2.csv`
- `docs/kazan_old_vs_new_full_years.csv`
- `docs/aoi_annual_quantiles_v2.csv`
