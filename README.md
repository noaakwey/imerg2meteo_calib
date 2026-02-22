# Dual-Source Precipitation Calibration & Temporal Downscaling

## Описание

Калибровка осадков IMERG V07 (NASA) и ERA5-Land (ECMWF) по данным **202 наземных метеостанций** Европейской России. Цель — расчёт R-фактора RUSLE и $I_{30}$ с точным воспроизведением экстремальных ливней.

## Результаты (медианы, кросс-валидация 2016–2021)

| Источник | Период | KGE (мес.) | KGE (сут.) | PBIAS |
|----------|--------|-----------|-----------|-------|
| IMERG — сырой | 2001– | 0.03 | −0.04 | −68.8% |
| **IMERG — калиброванный** | **2001–2025** | **0.68** | **0.49** | **−1.4%** |
| ERA5 — сырой | 1966– | 0.03 | −0.03 | −66.1% |
| **ERA5 — калиброванный** | **1966–2025** | **0.72** | **0.58** | **−6.9%** |

## Методы калибровки

### IMERG — Full-Distribution QM (1000 квантилей)
- Маппинг всего CDF (включая нули), тренировка 2001–2015, сезонное разбиение.
- Volume Scaling: `mean(station) / mean(IMERG_corr)` по обучающей выборке.
- Линейная экстраполяция правого хвоста (выше 99.9-го квантиля).

### ERA5-Land — Moving Window Full-Distribution QM
- Для каждого года *Y*: окно обучения `[Y−15, Y+15]` (~30 лет), сезонное разбиение.
- Full-Distribution QM (1000 квантилей) + Volume Scaling.
- Устраняет нестационарность 50-летнего ряда; обеспечивает ретроспективный охват 1966–2025.

## Быстрый старт

```bash
# Установка зависимостей
pip install -r requirements.txt

# Калибровка IMERG (2001–2025)
python src/main.py --dataset imerg

# Калибровка ERA5-Land (1966–2025)
python src/main.py --dataset era5land

# Анализ и графики
python src/analyze_results.py --dataset imerg
python src/analyze_results.py --dataset era5land
```

## Применение калибровки к растрам IMERG (R-фактор)

Скрипт `apply_qm_to_rasters.py` применяет QM-калибровку к **архиву квартальных GeoTIFF** IMERG, автоматически группируя кварталы по годам. Результат — годовые 1-часовые растры в `output/imerg_rfactor_calib/`.

```bash
# Обработать все годы (2001–2025)
python src/apply_qm_to_rasters.py

# Обработать один год
python src/apply_qm_to_rasters.py --year 2010

# Нестандартные пути
python src/apply_qm_to_rasters.py \
    --zip  "path/to/IMERG_RFACTOR_ANNUAL.zip" \
    --out  "path/to/output/" \
    --calib "path/to/calib_imerg/" \
    --meteo "path/to/meteo/"
```

**Алгоритм:**
1. Распаковка архива, группировка `YYYY_Q1…Q4` → год.
2. KNN-сопоставление пикселей растра → ближайшая метеостанция (cKDTree).
3. Full-Distribution QM (1000 квантилей) + Volume Scaling по сезонам (DJF / MAM / JJA / SON).
4. Запись годовых GeoTIFF (float32, LZW, tiled 256×256) с временны́ми метками бэндов.

> Скрипт идемпотентен: пропускает годы, для которых уже существует выходной файл.

## Структура проекта

```
src/
  main.py                      # Калибровка по точечным данным станций
  qm_calibration.py            # Full-Distribution QM, Moving Window QM
  apply_qm_to_rasters.py       # Пакетное применение QM к архиву растров IMERG
  batch_calibrate_rasters.py   # (устарел) Прежний прототип обработки растров
  data_loader.py               # Загрузка метеоданных и спутниковых рядов
  validation.py                # Метрики KGE, PBIAS
  analyze_results.py           # Графики и сводная статистика

data/
  meteo/срочные данные_осадки/   # CSV по станциям (cp866, разделитель ;)
  imerg/IMERG_STATIONS_3H/       # IMERG_V07_P3H_mm_YYYY_*.csv
  era5land/ERA5LAND_STATIONS_3H/ # ERA5Land_P3H_mm_YYYY.csv

output/
  calib_imerg/                   # *_calib.csv + validation_metrics_imerg.csv
  calib_era5land/                # *_calib.csv + validation_metrics_era5land.csv
  imerg_rfactor_calib/           # IMERG_V07_P1h_mm_YYYY_calib_qm.tif (годовые)
```

Данные не включены в репозиторий (`data/` в `.gitignore`).

## Подробная методика
[docs/methodology.md](docs/methodology.md)
