"""LOC95-vs-historical 30-yr moving-average (ΔT, ΔP) shift, overlaid on the raw/all-variants 2043
GCM point cloud.

Standalone analysis (not part of the LOCA-2/CMIP5 bivariate-normal pipeline): for every complete
30-yr window sliding one year at a time across the ingested LOC95 grid's full 1915-2021 daily
record (78 windows, starts 1915-1992), computes the same (ΔT, ΔP) convention the rest of the
package uses for the raw ("point-horizon") precip variant --

    D_tas = mean_annual_tave(window) - mean_annual_tave(historical baseline)
    D_pr  = (mean_annual_pr_total(window) - mean_annual_pr_total(historical baseline))
            / mean_annual_pr_total(historical baseline) * 100

where ``tave = (Tasmax + Tasmin) / 2`` and ``pr_total`` is a year's total precip -- against a FIXED
historical (Livneh/WGEN) baseline of 1989-2018. The rest of the package anchors its baseline on
``base_center=2006`` (i.e. 1992-2021), but the ingested historical grid's record stops at 2018 (see
``historical.py``), three years short of that window; 1989-2018 keeps a full 30-yr baseline rather
than shrinking it to 27 years.

The 78 points are overlaid as an extra scatter layer on the existing "raw + all variants" 2043 GCM
contour figure (``configs/cv-flow-weighted.yaml``'s ``outputs.variants.raw_novaravg``, the same
variant as ``loca2-biv-norm-prob-dt-dp-raw-novaravg-2043_cv-flow-weighted.svg``), so the LOC95
grid's own century-scale internal variability can be compared directly against the GCM ensemble's
projected 2043 spread.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: figures are written to files, never shown

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from .config import load_domain  # noqa: E402
from .distribution import BivariateNormalSurface  # noqa: E402
from .pipeline import LocaPipeline, _mean_tuple, _variant_title  # noqa: E402
from .plotting import ContourPlotter  # noqa: E402
from .precip_extremes import _historical_csv_path, _loc95_csv_path  # noqa: E402

WINDOW = 30
HIST_BASELINE_RANGE = (1989, 2018)  # historical grid's record stops at 2018; shifted back 3 yrs
LOC95_FULL_RANGE = (1915, 2021)     # LOC95 grid's full record
SHIFT_PERIOD = 2043
SHIFT_VARIANT = "raw_novaravg"
SHIFT_LABEL = "LOC95 (30-yr averages)"
SHIFT_COLOR = "blue"  # distinct from the purple/orange/red SSP scatter colors


def _read_daily_with_tave(csv_path: Path) -> pd.DataFrame:
    """Daily Year/'Pr (mm)'/'Tasmax (degC)'/'Tasmin (degC)' CSV -> one row per year: total precip,
    mean daily tave (``(Tasmax + Tasmin) / 2``)."""
    daily = pd.read_csv(csv_path)
    daily = daily.assign(tave=(daily["Tasmax (degC)"] + daily["Tasmin (degC)"]) / 2)
    return daily.groupby("Year", as_index=False).agg(pr_total=("Pr (mm)", "sum"), tave_mean=("tave", "mean"))


def _window_mean(annual: pd.DataFrame, start: int, end: int) -> tuple[float, float]:
    """(mean annual pr total, mean annual tave) across the inclusive ``[start, end]`` year window."""
    sub = annual[(annual["Year"] >= start) & (annual["Year"] <= end)]
    return float(sub["pr_total"].mean()), float(sub["tave_mean"].mean())


def compute_shift(cfg) -> pd.DataFrame:
    """(start_year, D_tas, D_pr) for every complete 30-yr window sliding 1 yr at a time across the
    LOC95 record, vs. the fixed 1989-2018 historical baseline."""
    hist_annual = _read_daily_with_tave(_historical_csv_path(cfg))
    loc_annual = _read_daily_with_tave(_loc95_csv_path(cfg))

    pr_base, tave_base = _window_mean(hist_annual, *HIST_BASELINE_RANGE)

    loc_start, loc_end = LOC95_FULL_RANGE
    rows = []
    for start in range(loc_start, loc_end - WINDOW + 2):
        end = start + WINDOW - 1
        pr_win, tave_win = _window_mean(loc_annual, start, end)
        rows.append({
            "start_year": start,
            "D_tas": tave_win - tave_base,
            "D_pr": (pr_win - pr_base) / pr_base * 100.0,
        })
    return pd.DataFrame(rows)


def _variant_cfg(cfg, variant_name: str):
    """The named ``outputs.variants`` override (e.g. ``raw_novaravg``), applied via
    ``dataclasses.replace`` exactly as ``pipeline.py`` does -- guarantees this matches the
    committed variant figure's config precisely (including its ``filter_gcms`` override)."""
    overrides = cfg.outputs["variants"][variant_name]
    return dataclasses.replace(cfg, **overrides)


def plot_shift_overlay(cfg=None, out_path=None):
    """The raw+novaravg 2043 GCM contour figure, with the LOC95-vs-historical moving-window points
    overlaid as an extra scatter layer + legend entry. Writes ``out_path`` (default:
    ``figures/loca2-biv-norm-prob-dt-dp-raw-novaravg-2043-loc95-shift_<domain>.svg``) and returns a
    dict with the shift points and the written figure path.
    """
    cfg = cfg or load_domain("cv-flow-weighted")
    variant_cfg = _variant_cfg(cfg, SHIFT_VARIANT)
    dist = LocaPipeline(variant_cfg)._compute_distribution()
    plot_cfg = dataclasses.replace(variant_cfg, grid=variant_cfg.plot_grid) if variant_cfg.plot_grid else variant_cfg
    frames = BivariateNormalSurface(plot_cfg).all_periods(dist["gcm_mean"], dist["covariances"])
    frame = next(f for f in frames if int(f["period"].iloc[0]) == SHIFT_PERIOD)
    points = dist["collapsed"][dist["collapsed"]["period"] == SHIFT_PERIOD]

    plotter = ContourPlotter(plot_cfg)
    fig, ax = plotter.filled_contour(
        frame, points, title=_variant_title(SHIFT_PERIOD, variant_cfg),
        mean=_mean_tuple(dist["gcm_mean"], SHIFT_PERIOD), cov=dist["covariances"][SHIFT_PERIOD],
    )

    shift = compute_shift(cfg)
    new_handle = ax.scatter(
        shift["D_pr"], shift["D_tas"], s=22, c=SHIFT_COLOR, marker="D",
        edgecolors="black", linewidths=0.4, zorder=6, label=SHIFT_LABEL,
    )
    legend = ax.get_legend()
    handles = list(legend.legend_handles) if legend else []
    labels = [t.get_text() for t in legend.get_texts()] if legend else []
    ax.legend(
        handles + [new_handle], labels + [SHIFT_LABEL],
        loc="best", fontsize=8, framealpha=0.9, numpoints=2,
    )

    if out_path is None:
        out_path = (
            cfg.paths["figures_dir"]
            / f"loca2-biv-norm-prob-dt-dp-raw-novaravg-{SHIFT_PERIOD}-loc95-shift_{cfg.name}.svg"
        )
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return {"shift": shift, "figure": out_path}
