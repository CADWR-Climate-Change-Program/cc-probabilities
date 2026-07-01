"""Time series of the Pearson correlation coefficient between DeltaT and DeltaP.

A standalone analysis (not part of the bivariate-normal stress-test pipeline): for each of the four
distribution variants ``ccprob`` already produces for a domain -- the default lm/variant-averaged
config plus every override in ``outputs.variants`` (novaravg, raw, raw_novaravg) -- extracts the
Pearson correlation coefficient ``r = cov(DT,DP) / (sd(DT) * sd(DP))`` from that variant's per-period
2x2 covariance (``ClimateDeltas.period_covariances``, the same matrix the crosshair overlay and
``gcm_sigs`` CSV are built from), and plots all four as one time series (period on the x-axis) so how
strongly ΔT/ΔP correlate -- and how much that depends on the averaging/precip-change choice -- can be
compared at a glance. A domain with no ``outputs.variants`` configured just yields a single line (the
default variant).
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: figures are written to files, never shown

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from .pipeline import LocaPipeline

# (pr_token, avg_token) -> line style; lm=teal, raw=orange, varavg=solid, novaravg=dashed
VARIANT_STYLES = {
    ("lm", "varavg"): {"color": "#1b9e77", "linestyle": "-"},
    ("lm", "novaravg"): {"color": "#1b9e77", "linestyle": "--"},
    ("raw", "varavg"): {"color": "#d95f02", "linestyle": "-"},
    ("raw", "novaravg"): {"color": "#d95f02", "linestyle": "--"},
}


def _variant_key(variant_cfg) -> str:
    avg = "varavg" if variant_cfg.filter_nmem_gcms else "novaravg"
    pr = "lm" if variant_cfg.pr_type == "D_pr_lm" else "raw"
    return f"{pr}-{avg}"


def _variant_configs(cfg) -> dict:
    """{'lm-varavg': cfg, 'lm-novaravg': overridden_cfg, ...} -- the base config plus every
    ``outputs.variants`` override, keyed the same way the figure filenames already are."""
    variants = {_variant_key(cfg): cfg}
    for overrides in cfg.outputs.get("variants", {}).values():
        variant_cfg = dataclasses.replace(cfg, **overrides)
        variants[_variant_key(variant_cfg)] = variant_cfg
    return variants


def pearson_r(cov) -> float:
    """Pearson correlation coefficient of (DT, DP) from their 2x2 covariance matrix."""
    sd_t = cov[0][0] ** 0.5
    sd_p = cov[1][1] ** 0.5
    if sd_t <= 0 or sd_p <= 0:
        return float("nan")
    return float(cov[0][1] / (sd_t * sd_p))


def compute(cfg) -> pd.DataFrame:
    """Per (variant, period): the Pearson r between DT and DP. Columns: variant, period, r."""
    rows = []
    for label, variant_cfg in _variant_configs(cfg).items():
        covs = LocaPipeline(variant_cfg)._compute_distribution()["covariances"]
        for period, cov in covs.items():
            rows.append({"variant": label, "period": period, "r": pearson_r(cov)})
    return pd.DataFrame(rows).sort_values(["variant", "period"]).reset_index(drop=True)


def plot_time_series(df: pd.DataFrame, *, title: str, out_path) -> Path:
    fig, ax = plt.subplots(figsize=(8, 5))
    for label, sub in sorted(df.groupby("variant")):
        pr, avg = label.split("-", 1)
        style = VARIANT_STYLES.get((pr, avg), {})
        sub = sub.sort_values("period")
        ax.plot(sub["period"], sub["r"], marker="o", markersize=3, linewidth=1.6, label=label, **style)

    ax.axhline(0.0, color="gray", linewidth=0.8, linestyle=":")
    ax.set_xlabel("Period")
    ax.set_ylabel("Pearson correlation coefficient (ΔT, ΔP)")
    ax.set_title(title, fontsize=11)
    ax.legend(loc="best", fontsize=9)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def run(cfg, out_dir=None, write: bool = True) -> dict:
    """Compute and (optionally) write the results CSV + the time-series figure."""
    df = compute(cfg)
    written: dict = {"csv": None, "figure": None}
    if write:
        processed = Path(out_dir) if out_dir else cfg.paths["processed_dir"]
        figures = Path(out_dir) if out_dir else cfg.paths["figures_dir"]
        processed.mkdir(parents=True, exist_ok=True)
        csv_path = processed / f"dt_dp_correlation_{cfg.name}.csv"
        df.to_csv(csv_path, index=False)
        written["csv"] = csv_path

        figures.mkdir(parents=True, exist_ok=True)
        fig_path = figures / f"dt_dp_correlation_{cfg.name}.svg"
        title = f"{cfg.name}: Pearson correlation between ΔT and ΔP"
        written["figure"] = plot_time_series(df, title=title, out_path=fig_path)
    return {"results": df, **written}
