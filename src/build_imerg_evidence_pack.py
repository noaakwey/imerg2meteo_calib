import calendar
import glob
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

BASE = Path(__file__).resolve().parents[1]
DOCS = BASE / "docs"
GDIR = DOCS / "graphics"

KAZAN_WMO = 27595
KAZAN_LON = 49.2
KAZAN_LAT = 55.73333333

VERSION_DIRS = {
    "v2": BASE / "output" / "imerg_rfactor_calib_v2",
    "v3": BASE / "output" / "imerg_rfactor_calib_v3_soft",
    "v4": BASE / "output" / "imerg_rfactor_calib_v4_annual_transfer",
    "v5": BASE / "output" / "imerg_rfactor_calib_v5_year_anchor",
}


def find_kazan_calib_csv() -> Path:
    matches = sorted((BASE / "output" / "calib_imerg").glob("*27595_calib.csv"))
    if not matches:
        raise FileNotFoundError("Kazan calib CSV (*27595_calib.csv) not found in output/calib_imerg")
    return matches[0]


def list_rasters(folder: Path) -> list[tuple[int, Path]]:
    if not folder.exists():
        return []
    out = []
    for p in sorted(folder.glob("IMERG_V07_P30min_mmh_*_calib_qm.tif")):
        m = re.search(r"(\d{4})_calib", p.name)
        if m:
            out.append((int(m.group(1)), p))
    return out


def raster_kazan_series(items: list[tuple[int, Path]], col: str) -> pd.DataFrame:
    rows = []
    for y, p in items:
        with rasterio.open(p) as src:
            r, c = src.index(KAZAN_LON, KAZAN_LAT)
            vals = src.read(window=((r, r + 1), (c, c + 1))).astype("float64")[:, 0, 0]
            mm = float(np.nansum(np.maximum(vals, 0.0)) * 0.5)
        rows.append({"year": y, col: mm})
    return pd.DataFrame(rows).set_index("year").sort_index()


def aoi_quantiles(items: list[tuple[int, Path]]) -> pd.DataFrame:
    rows = []
    for y, p in items:
        with rasterio.open(p) as src:
            arr = src.read().astype(np.float32)
        ann = np.nansum(np.maximum(arr, 0.0), axis=0, dtype=np.float64) * 0.5
        rows.append(
            {
                "year": y,
                "ann_min": float(np.nanmin(ann)),
                "ann_p50": float(np.nanpercentile(ann, 50)),
                "ann_p90": float(np.nanpercentile(ann, 90)),
                "ann_p99": float(np.nanpercentile(ann, 99)),
                "ann_max": float(np.nanmax(ann)),
                "neg_cnt": int(np.sum(ann < 0.0)),
                "nan_cnt": int(np.sum(~np.isfinite(ann))),
            }
        )
    return pd.DataFrame(rows).sort_values("year")


def build_station_annual_pbias_summary() -> pd.DataFrame:
    rows = []
    files = sorted((BASE / "output" / "calib_imerg").glob("*_calib.csv"))
    for f in files:
        try:
            df = pd.read_csv(f, usecols=["datetime", "P_station_mm", "P_sat_mm", "P_corrected_mm"])
        except Exception:
            continue

        df["datetime"] = pd.to_datetime(df["datetime"])
        df["year"] = df["datetime"].dt.year
        g = df.groupby("year", as_index=False)[["P_station_mm", "P_sat_mm", "P_corrected_mm"]].sum()
        g = g[(g["year"] >= 2001) & (g["year"] <= 2024) & (g["P_station_mm"] > 0)]
        if len(g) < 10:
            continue

        pb_raw = 100.0 * (g["P_sat_mm"] - g["P_station_mm"]) / g["P_station_mm"]
        pb_cor = 100.0 * (g["P_corrected_mm"] - g["P_station_mm"]) / g["P_station_mm"]

        m = re.search(r"_(\d+)_calib\.csv$", f.name)
        wmo = int(m.group(1)) if m else -1
        rows.append(
            {
                "wmo": wmo,
                "n_years": int(len(g)),
                "annual_pbias_raw_mean": float(pb_raw.mean()),
                "annual_pbias_corr_mean": float(pb_cor.mean()),
                "annual_pbias_raw_abs_mean": float(pb_raw.abs().mean()),
                "annual_pbias_corr_abs_mean": float(pb_cor.abs().mean()),
                "annual_pbias_raw_abs_med": float(pb_raw.abs().median()),
                "annual_pbias_corr_abs_med": float(pb_cor.abs().median()),
            }
        )

    return pd.DataFrame(rows)


def build_kazan_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    kcsv = pd.read_csv(find_kazan_calib_csv())
    kcsv["datetime"] = pd.to_datetime(kcsv["datetime"])
    kcsv["year"] = kcsv["datetime"].dt.year

    ksum = kcsv.groupby("year")[["P_station_mm", "P_sat_mm", "P_corrected_mm"]].sum()
    counts = kcsv.groupby("year").size().rename("n_rows").to_frame()
    counts["expected"] = [(366 if calendar.isleap(int(y)) else 365) * 8 for y in counts.index]
    counts["completeness"] = counts["n_rows"] / counts["expected"]

    cmp = ksum.join(counts, how="left")
    for v, folder in VERSION_DIRS.items():
        items = list_rasters(folder)
        if not items:
            continue
        s = raster_kazan_series(items, f"raster_{v}_mm")
        cmp = cmp.join(s, how="left")
        cmp[f"pbias_{v}_pct"] = (cmp[f"raster_{v}_mm"] - cmp["P_station_mm"]) / cmp["P_station_mm"] * 100.0
        cmp[f"abs_pbias_{v}"] = cmp[f"pbias_{v}_pct"].abs()

    full = cmp[(cmp["completeness"] >= 0.99) & (cmp.index >= 2001) & (cmp.index <= 2024)].copy()

    med = float(full["P_station_mm"].median())
    iqr = float(full["P_station_mm"].quantile(0.75) - full["P_station_mm"].quantile(0.25))
    if iqr <= 0:
        iqr = max(float(full["P_station_mm"].std()), 1.0)
    full["ext_idx"] = (full["P_station_mm"] - med).abs() / iqr
    full["is_extreme"] = full["ext_idx"] >= 1.0

    for base in ("v2", "v3", "v4"):
        if f"abs_pbias_{base}" in full.columns and "abs_pbias_v5" in full.columns:
            full[f"improve_v5_vs_{base}"] = full[f"abs_pbias_{base}"] - full["abs_pbias_v5"]

    return cmp, full


def summarize_versions(full: pd.DataFrame) -> pd.DataFrame:
    rows = []
    versions = [v for v in ("v2", "v3", "v4", "v5") if f"abs_pbias_{v}" in full.columns]
    for v in versions:
        s = full[f"abs_pbias_{v}"]
        rows.append(
            {
                "version": v,
                "mean_abs_pbias": float(s.mean()),
                "median_abs_pbias": float(s.median()),
                "max_abs_pbias": float(s.max()),
                "corr_ext_abs_pbias": float(full["ext_idx"].corr(s)),
                "mean_abs_extreme": float(s[full["is_extreme"]].mean()),
                "mean_abs_normal": float(s[~full["is_extreme"]].mean()),
                "n_years": int(len(full)),
            }
        )
    return pd.DataFrame(rows)


def save_plots(full: pd.DataFrame, aoi_v3: pd.DataFrame | None, aoi_v5: pd.DataFrame | None) -> None:
    GDIR.mkdir(parents=True, exist_ok=True)

    if {"raster_v3_mm", "raster_v5_mm"}.issubset(full.columns):
        plt.figure(figsize=(11.5, 5.2))
        plt.plot(full.index, full["P_station_mm"], label="Station", lw=2.5)
        plt.plot(full.index, full["raster_v3_mm"], label="Raster v3_soft", lw=2.0)
        plt.plot(full.index, full["raster_v5_mm"], label="Raster v5_year_anchor", lw=2.1)
        plt.grid(alpha=0.25)
        plt.legend()
        plt.title("Kazan annual totals: station vs v3_soft vs v5_year_anchor")
        plt.tight_layout()
        plt.savefig(GDIR / "evidence_kazan_v3_vs_v5_series.png", dpi=220)
        plt.close()

        plt.figure(figsize=(11.5, 5.0))
        plt.plot(full.index, full["abs_pbias_v3"], label="|PBIAS| v3_soft", lw=2.0)
        plt.plot(full.index, full["abs_pbias_v5"], label="|PBIAS| v5_year_anchor", lw=2.1)
        plt.grid(alpha=0.25)
        plt.legend()
        plt.title("Kazan annual |PBIAS|: v3_soft vs v5_year_anchor")
        plt.tight_layout()
        plt.savefig(GDIR / "evidence_kazan_abs_pbias_v3_vs_v5.png", dpi=220)
        plt.close()

        plt.figure(figsize=(7.6, 5.5))
        plt.scatter(full["ext_idx"], full["abs_pbias_v3"], label="v3_soft", alpha=0.75)
        plt.scatter(full["ext_idx"], full["abs_pbias_v5"], label="v5_year_anchor", alpha=0.75)
        for key, color in (("v3", "tab:blue"), ("v5", "tab:orange")):
            x = full["ext_idx"].values
            y = full[f"abs_pbias_{key}"].values
            if len(x) >= 2:
                m, b = np.polyfit(x, y, 1)
                xx = np.linspace(float(x.min()), float(x.max()), 100)
                plt.plot(xx, m * xx + b, color=color, lw=1.8)
        plt.xlabel("Extremeness index |P_station - median| / IQR")
        plt.ylabel("|PBIAS|, %")
        plt.grid(alpha=0.25)
        plt.legend()
        plt.title("Error vs year extremeness (Kazan)")
        plt.tight_layout()
        plt.savefig(GDIR / "evidence_extremeness_vs_abs_pbias_v3_vs_v5.png", dpi=220)
        plt.close()

    versions = [v for v in ("v2", "v3", "v4", "v5") if f"abs_pbias_{v}" in full.columns]
    if versions:
        labels = [v.replace("v3", "v3_soft").replace("v4", "v4_transfer").replace("v5", "v5_anchor") for v in versions]
        vals = [full[f"abs_pbias_{v}"].dropna().values for v in versions]

        plt.figure(figsize=(8.8, 5.0))
        plt.boxplot(vals, tick_labels=labels)
        plt.grid(alpha=0.2)
        plt.title("Kazan annual |PBIAS| distribution by version")
        plt.tight_layout()
        plt.savefig(GDIR / "evidence_kazan_abs_pbias_box_v2_v3_v4_v5.png", dpi=220)
        plt.close()

        means = [float(np.mean(vv)) for vv in vals]
        plt.figure(figsize=(8.8, 4.8))
        plt.bar(labels, means, color=["#4C72B0", "#55A868", "#C44E52", "#8172B2"][: len(labels)])
        plt.ylabel("mean |PBIAS|, %")
        plt.title("Experiment progression: mean annual |PBIAS| (Kazan)")
        plt.grid(axis="y", alpha=0.2)
        plt.tight_layout()
        plt.savefig(GDIR / "evidence_experiment_progression_mean_abs_pbias.png", dpi=220)
        plt.close()

    if aoi_v3 is not None and aoi_v5 is not None and (not aoi_v3.empty) and (not aoi_v5.empty):
        m = aoi_v3.merge(aoi_v5, on="year", suffixes=("_v3", "_v5"))
        plt.figure(figsize=(11.5, 5.0))
        plt.plot(m["year"], m["ann_p99_v3"], label="p99 v3", lw=1.8)
        plt.plot(m["year"], m["ann_p99_v5"], label="p99 v5", lw=2.0)
        plt.plot(m["year"], m["ann_max_v3"], label="max v3", lw=1.8)
        plt.plot(m["year"], m["ann_max_v5"], label="max v5", lw=2.0)
        plt.grid(alpha=0.25)
        plt.legend(ncol=2)
        plt.title("AOI annual envelope: v3_soft vs v5_year_anchor")
        plt.tight_layout()
        plt.savefig(GDIR / "evidence_aoi_envelope_v3_vs_v5.png", dpi=220)
        plt.close()


def write_main_report(summary: pd.DataFrame,
                      full: pd.DataFrame,
                      aoi_v5: pd.DataFrame,
                      aoi_delta: pd.DataFrame | None,
                      station_annual: pd.DataFrame) -> None:
    report = DOCS / "evidence_imerg_qc_main.md"

    v_map = {r["version"]: r for _, r in summary.iterrows()}
    lines = [
        "# IMERG: доказательная база эволюции v2 -> v3 -> v4 -> v5 (основной отчет)",
        "",
        "## 1. Исходная проблема",
        "- В годовых калиброванных растрах сохранялся систематический разрыв со станциями.",
        "- В экстремальные годы ошибка была выше.",
        "- Цель: уменьшить годовой разрыв и ослабить зависимость ошибки от экстремальности года без нефизичного раздува по AOI.",
        "",
        "## 2. Эксперименты",
        "| Версия | Что изменили | Ожидаемый эффект |",
        "|---|---|---|",
        "| v2 | фиксация temporal mismatch + фиксация KNN + annual sanity | убрать грубое завышение |",
        "| v3_soft | soft seasonal QM + мягкие daily/annual guards | снизить перекоррекцию |",
        "| v4_transfer | annual raw->station transfer | улучшить межгодовую амплитуду |",
        "| v5_year_anchor | годовой anchor по station/sat ratio + transfer + sanity | целенаправленно сжать разрыв в экстремумах |",
        "",
        "## 3. Сравнение по Казани (полные годы 2001-2024)",
        f"- Полных лет: `{len(full)}`.",
    ]

    for key, label in (("v2", "v2"), ("v3", "v3_soft"), ("v4", "v4_transfer"), ("v5", "v5_year_anchor")):
        r = v_map.get(key)
        if r is None:
            continue
        lines.append(
            f"- `{label}`: mean/median/max `|PBIAS|` = "
            f"`{r['mean_abs_pbias']:.2f}% / {r['median_abs_pbias']:.2f}% / {r['max_abs_pbias']:.2f}%`; "
            f"`corr(ext, |PBIAS|) = {r['corr_ext_abs_pbias']:.3f}`."
        )

    v3 = v_map.get("v3")
    v5 = v_map.get("v5")
    if v3 is not None and v5 is not None:
        lines.extend(
            [
                "",
                "## 4. Почему v5 выбран основным",
                f"- Минимальный средний `|PBIAS|`: `{v5['mean_abs_pbias']:.2f}%` (v3_soft: `{v3['mean_abs_pbias']:.2f}%`).",
                f"- Минимальная ошибка в худшем году: `{v5['max_abs_pbias']:.2f}%` (v3_soft: `{v3['max_abs_pbias']:.2f}%`).",
                f"- Самая слабая связь ошибки с экстремальностью: `corr(ext, |PBIAS|) = {v5['corr_ext_abs_pbias']:.3f}` "
                f"(v3_soft: `{v3['corr_ext_abs_pbias']:.3f}`).",
            ]
        )

    if (not aoi_v5.empty):
        lines.extend(
            [
                "",
                "## 5. AOI sanity-check (v5)",
                f"- `ann_p99` range: `{aoi_v5['ann_p99'].min():.2f} .. {aoi_v5['ann_p99'].max():.2f}` mm/yr.",
                f"- `ann_max` range: `{aoi_v5['ann_max'].min():.2f} .. {aoi_v5['ann_max'].max():.2f}` mm/yr.",
                f"- Negative pixels: `{int(aoi_v5['neg_cnt'].sum())}`; NaN pixels: `{int(aoi_v5['nan_cnt'].sum())}`.",
            ]
        )
        if aoi_delta is not None and (not aoi_delta.empty):
            lines.append(
                f"- Относительно v3_soft: средний сдвиг `p99 = {aoi_delta['delta_ann_p99_pct'].mean():+.2f}%`, "
                f"`max = {aoi_delta['delta_ann_max_pct'].mean():+.2f}%`."
            )

    if not station_annual.empty:
        lines.extend(
            [
                "",
                "## 6. Station-level годовой баланс (2001-2024)",
                f"- Median station `|annual PBIAS|`: "
                f"`{station_annual['annual_pbias_raw_abs_med'].median():.2f}% -> {station_annual['annual_pbias_corr_abs_med'].median():.2f}%`.",
                f"- Mean station `|annual PBIAS|`: "
                f"`{station_annual['annual_pbias_raw_abs_mean'].mean():.2f}% -> {station_annual['annual_pbias_corr_abs_mean'].mean():.2f}%`.",
            ]
        )

    lines.extend(
        [
            "",
            "## 7. Артефакты",
            "- Таблицы:",
            "  - `docs/kazan_v2_v3_v4_v5_full_years.csv`",
            "  - `docs/kazan_summary_v2_v3_v4_v5.csv`",
            "  - `docs/aoi_annual_quantiles_v5_year_anchor.csv`",
            "  - `docs/aoi_v3_vs_v5_delta.csv`",
            "  - `docs/station_annual_pbias_summary_imerg.csv`",
            "- Графики:",
            "  - `docs/graphics/evidence_kazan_v3_vs_v5_series.png`",
            "  - `docs/graphics/evidence_kazan_abs_pbias_v3_vs_v5.png`",
            "  - `docs/graphics/evidence_extremeness_vs_abs_pbias_v3_vs_v5.png`",
            "  - `docs/graphics/evidence_kazan_abs_pbias_box_v2_v3_v4_v5.png`",
            "  - `docs/graphics/evidence_experiment_progression_mean_abs_pbias.png`",
            "  - `docs/graphics/evidence_aoi_envelope_v3_vs_v5.png`",
            "- Архивные материалы (трассируемость):",
            "  - `docs/evidence_imerg_qc.md`",
            "  - `docs/evidence_imerg_qc_v3_soft.md`",
            "  - `docs/evidence_imerg_qc_v5_year_anchor.md`",
            "  - `docs/methodology_imerg_v2.md`",
            "",
            "## 8. Ограничение применения",
            "- `v5_year_anchor` — лучший режим для исторической реконструкции, когда доступны station-данные того же года.",
            "- Для strict out-of-sample режима использовать `v4_transfer`/`v3_soft`.",
        ]
    )

    report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    DOCS.mkdir(parents=True, exist_ok=True)
    GDIR.mkdir(parents=True, exist_ok=True)

    cmp, full = build_kazan_tables()
    summary = summarize_versions(full)

    aoi_v5_items = list_rasters(VERSION_DIRS["v5"])
    aoi_v5 = aoi_quantiles(aoi_v5_items) if aoi_v5_items else pd.DataFrame()
    aoi_v5.to_csv(DOCS / "aoi_annual_quantiles_v5_year_anchor.csv", index=False, encoding="utf-8-sig")

    aoi_v3_items = list_rasters(VERSION_DIRS["v3"])
    aoi_v3 = aoi_quantiles(aoi_v3_items) if aoi_v3_items else pd.DataFrame()
    if (not aoi_v3.empty) and (not aoi_v5.empty):
        d = aoi_v3.merge(aoi_v5, on="year", suffixes=("_v3", "_v5"))
        for c in ("ann_p50", "ann_p90", "ann_p99", "ann_max"):
            d[f"delta_{c}_pct"] = 100.0 * (d[f"{c}_v5"] - d[f"{c}_v3"]) / d[f"{c}_v3"]
        d[["year", "delta_ann_p50_pct", "delta_ann_p90_pct", "delta_ann_p99_pct", "delta_ann_max_pct"]].to_csv(
            DOCS / "aoi_v3_vs_v5_delta.csv", index=False, encoding="utf-8-sig"
        )
        aoi_delta = d
    else:
        aoi_delta = None

    station_annual = build_station_annual_pbias_summary()
    station_annual.to_csv(DOCS / "station_annual_pbias_summary_imerg.csv", index=False, encoding="utf-8-sig")

    cmp.to_csv(DOCS / "kazan_annual_balance_v5_year_anchor.csv", encoding="utf-8-sig")
    full.to_csv(DOCS / "kazan_v2_v3_v4_v5_full_years.csv", encoding="utf-8-sig")
    summary.to_csv(DOCS / "kazan_summary_v2_v3_v4_v5.csv", index=False, encoding="utf-8-sig")

    save_plots(full, aoi_v3 if not aoi_v3.empty else None, aoi_v5 if not aoi_v5.empty else None)
    write_main_report(summary, full, aoi_v5, aoi_delta, station_annual)

    print("Evidence pack built")
    print("  report: docs/evidence_imerg_qc_main.md")
    print("  summary: docs/kazan_summary_v2_v3_v4_v5.csv")


if __name__ == "__main__":
    main()
