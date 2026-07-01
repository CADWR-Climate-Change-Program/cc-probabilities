"""Minimum sliding-window N-yr average annual precipitation, by mvs, vs. the canonical historical
record.

A standalone analysis (not part of the bivariate-normal stress-test pipeline): for each
model_variant_ssp's raw annual precipitation timeseries, within each of three calendar periods
(the LOCA-2 historical run, and the 2021-2060 / 2061-2100 projection eras) and each window size in
``WINDOWS``, slide an overlapping N-yr window across that period's years and take the minimum
average across all window positions -- the driest N-yr stretch that model/scenario/era projects.
Overlapping (rather than non-overlapping/block) windows are used deliberately: a fixed-origin,
non-overlapping partition can split a real drought across a block boundary and understate it, and
its minimum is only guaranteed monotonic across window sizes that evenly divide one another (e.g.
a 15-yr and a 20-yr partition of the same record are NOT guaranteed comparable). A sliding window
searches every possible N-yr stretch instead, which is both the more standard drought-analysis
convention and far more robust to this artifact in practice. Windows are always confined to their
own calendar period (each period's years are filtered out *before* windowing), so a window can
never span across a period boundary. Plotted as one empirical-CDF figure per window size: a solid
line per period (black/blue/orange), each labeled with its sample size (the number of
model_variant_ssp series it draws from), plus a single fixed-color vertical line ("Historical
Livneh") marking the canonical Livneh historical record's own driest N-yr stretch
(``data/historical/``, produced by ``historical.py``), computed the same sliding-window way.

Two domains are supported: ``loca2-flow`` (reuses the existing annual-flow ``Ensemble``, already contiguous
1950-2100) and ``loca2-basin`` (aggregates the raw daily CSVs to annual here, since no pre-aggregated-annual
basin product exists -- only 30-yr-window averages, which the rest of the package uses instead).
"""

from __future__ import annotations

import glob
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: figures are written to files, never shown

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .config import REPO_ROOT
from .io import FilenameParser
from .projections import Ensemble

WINDOWS = (2, 5, 10, 15, 20)  # sliding-window sizes (years) for the minimum-average metric
HIST_END_YEAR = 2014    # inclusive; matches the LOCA-2 historical/projection file split

# period -> inclusive (start, end) year range; start=None means "from the beginning of the record"
PERIOD_RANGES = {
    "historical": (None, HIST_END_YEAR),
    "early_projection": (2021, 2060),
    "late_projection": (2061, 2100),
}
PERIOD_LABELS = {
    "historical": "Historical (LOCA-2)",
    "early_projection": "2021-2060",
    "late_projection": "2061-2100",
}
PERIOD_COLORS = {"historical": "black", "early_projection": "tab:blue", "late_projection": "tab:orange"}

DAILY_BASIN_GLOB = "data/loca2/loca2-projections-daily-basin/*_{basin_id}-*.csv"
HISTORICAL_DIR = REPO_ROOT / "data" / "historical"


def load_raw_annual(cfg) -> pd.DataFrame:
    """Raw annual precip per mvs, contiguous across historical + projection years.

    Columns: mvs, model, variant, ssp, y, pr.
    """
    if cfg.source_kind == "loca2-flow":
        return Ensemble.from_config(cfg).long_frame()[["mvs", "model", "variant", "ssp", "y", "pr"]]
    if cfg.source_kind == "loca2-basin":
        return _load_basin_daily_annual(cfg)
    raise ValueError(f"load_raw_annual does not support source_kind={cfg.source_kind!r}")


def _load_basin_daily_annual(cfg) -> pd.DataFrame:
    """Daily basin CSVs -> per-mvs annual precip, with the historical series replicated once per
    SSP the model+variant actually ran (mirrors Ensemble.from_config's loca2-flow rule -- the daily
    basin files split historical/projection by filename the same way the annual-flow files do)."""
    parser = FilenameParser(["model", "ssp", "variant"])
    pattern = str(REPO_ROOT / DAILY_BASIN_GLOB.format(basin_id=cfg.basin_id))
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No daily basin files match {pattern}")
    meta = [parser.parse(f) for f in files]
    meta_df = pd.DataFrame(meta)
    nonhist = meta_df[meta_df["ssp"] != "historical"]

    pieces = []
    for path, m in zip(files, meta):
        model, ssp, variant = m["model"], m["ssp"], m["variant"]
        annual = pd.read_csv(path).groupby("Year", as_index=False)["Pr (mm)"].sum()
        annual.columns = ["y", "pr"]
        ssps = (
            [ssp] if ssp != "historical"
            else nonhist[(nonhist.model == model) & (nonhist.variant == variant)]["ssp"].unique()
        )
        for s in ssps:
            c = annual.copy()
            c.insert(0, "mvs", FilenameParser.mvs(model, variant, s))
            c.insert(1, "model", model)
            c.insert(2, "variant", variant)
            c.insert(3, "ssp", s)
            pieces.append(c)

    frame = pd.concat(pieces, axis=0, ignore_index=True)
    return frame.sort_values(["mvs", "y"]).reset_index(drop=True)


def minimum_window_averages(frame: pd.DataFrame, windows=WINDOWS) -> pd.DataFrame:
    """Per mvs, per window: the minimum average annual precip across all overlapping (sliding)
    ``window``-yr windows in the series. A series shorter than ``window`` is excluded for that
    window size.
    """
    rows = []
    for mvs, data in frame.groupby("mvs"):
        pr = data.sort_values("y")["pr"].to_numpy()
        for window in windows:
            n_windows = len(pr) - window + 1
            if n_windows <= 0:
                continue
            cumsum = np.cumsum(np.insert(pr, 0, 0.0))
            window_avgs = (cumsum[window:] - cumsum[:-window]) / window
            rows.append({
                "mvs": mvs, "window": window, "n_windows": n_windows,
                "min_window_avg_pr": window_avgs.min(),
            })
    return pd.DataFrame(rows)


def compute(cfg, windows=WINDOWS) -> pd.DataFrame:
    """Minimum sliding-window average precip for every mvs, in each of the three periods and window
    sizes. Each period's years are filtered out first, so windows never cross a period boundary."""
    raw = load_raw_annual(cfg)
    pieces = []
    for period, (start, end) in PERIOD_RANGES.items():
        sub = raw[raw["y"] <= end] if start is None else raw[(raw["y"] >= start) & (raw["y"] <= end)]
        r = minimum_window_averages(sub, windows=windows)
        r.insert(0, "period", period)
        pieces.append(r)
    return pd.concat(pieces, axis=0, ignore_index=True)


def _historical_csv_path(cfg) -> Path:
    if cfg.source_kind == "loca2-flow":
        return HISTORICAL_DIR / "historical_daily_cv-flow-weighted.csv"
    if cfg.source_kind == "loca2-basin":
        return HISTORICAL_DIR / f"historical_daily_basin-{cfg.basin_id}.csv"
    raise ValueError(f"no canonical historical data for source_kind={cfg.source_kind!r}")


def wgen_min_window_average(cfg, window: int) -> float:
    """The canonical Livneh historical record's own minimum ``window``-yr sliding-window average."""
    daily = pd.read_csv(_historical_csv_path(cfg))
    annual = daily.groupby("Year", as_index=False)["Pr (mm)"].sum()
    annual.columns = ["y", "pr"]
    annual.insert(0, "mvs", "wgen_historical")
    result = minimum_window_averages(annual, windows=(window,))
    return float(result["min_window_avg_pr"].iloc[0])


WGEN_LINE_COLOR = "red"


def _assumptions_subtitle(window: int) -> str:
    """Short summary of this plot's analytical assumptions, for a small-text subtitle (one line per
    clause -- broken at the semicolon so each fits comfortably under the title)."""
    return (
        f"Minimum {window}-yr sliding window per model-variant-ssp\n"
        "dashed line = same metric on the Livneh historical record"
    )


def plot_cdf(results: pd.DataFrame, wgen_value: float, *, title: str, out_path) -> Path:
    """Empirical CDF of min_window_avg_pr for a single window size (``results`` pre-filtered to
    one ``window`` value): one solid line per period, each labeled with its sample size (the number
    of model_variant_ssp series drawn from), plus a single vertical reference line marking the
    canonical Livneh historical record's own value (fixed color, independent of period). A
    left-aligned, small-text subtitle below the title spells out the analytical assumptions.
    """
    window = int(results["window"].iloc[0])
    fig, ax = plt.subplots(figsize=(7, 5))

    for period in ("historical", "early_projection", "late_projection"):
        sub = results[results["period"] == period]
        vals = np.sort(sub["min_window_avg_pr"].to_numpy())
        if len(vals) == 0:
            continue
        color = PERIOD_COLORS[period]
        cdf = np.arange(1, len(vals) + 1) / len(vals)
        label = f"{PERIOD_LABELS[period]} (n={len(vals)})"
        ax.step(vals, cdf, where="post", color=color, linewidth=1.8, label=label)

    ax.axvline(wgen_value, color=WGEN_LINE_COLOR, linestyle="--", linewidth=1.5, label="Historical Livneh")

    ax.legend(loc="lower right")
    ax.set_xlabel(f"Minimum {window}-yr sliding-window average annual precipitation (mm)")
    ax.set_ylabel("Cumulative probability")
    ax.set_title(_assumptions_subtitle(window), loc="left", fontsize=8, style="italic", color="dimgray")
    ax.text(0.0, 1.12, title, transform=ax.transAxes, fontsize=12, ha="left", va="bottom")
    ax.set_ylim(0, 1)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def run(cfg, out_dir=None, write: bool = True, windows=WINDOWS) -> dict:
    """Compute and (optionally) write the results CSV (all window sizes) and one CDF figure per
    window size in ``windows``."""
    results = compute(cfg, windows=windows)
    wgen_values = {window: wgen_min_window_average(cfg, window=window) for window in windows}
    written: dict = {"csv": None, "figures": []}
    if write:
        processed = Path(out_dir) if out_dir else cfg.paths["processed_dir"]
        figures = Path(out_dir) if out_dir else cfg.paths["figures_dir"]
        processed.mkdir(parents=True, exist_ok=True)
        csv_path = processed / f"precip_block_minima_{cfg.name}.csv"
        results.to_csv(csv_path, index=False)
        written["csv"] = csv_path
        for window in windows:
            sub = results[results["window"] == window]
            title = f"{cfg.name}: {window}-yr precipitation minima"
            fig_path = figures / f"precip_block_minima_cdf_{window}yr_{cfg.name}.svg"
            written["figures"].append(plot_cdf(sub, wgen_values[window], title=title, out_path=fig_path))
    return {"results": results, "wgen_values": wgen_values, **written}
