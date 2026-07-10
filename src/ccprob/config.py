"""Domain configuration: replaces the ~25 hand-edited globals of the legacy scripts.

A ``DomainConfig`` fully describes one run (spatial domain + grid + filters + paths). Configs are
loaded from ``configs/<name>.yaml`` and all paths are resolved from the repository root, so there
is no ``setwd()`` or ``../`` working-directory assumption anywhere.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# src/ccprob/config.py -> parents[2] == repo root (works for an editable `pip install -e .`)
REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIGS_DIR = REPO_ROOT / "configs"


def r_seq(start: float, stop: float, by: float) -> np.ndarray:
    """Replicate R ``seq(start, stop, by=by)`` without float drift (endpoints inclusive)."""
    n = int(round((stop - start) / by))
    return start + by * np.arange(n + 1)


@dataclass(frozen=True)
class GridSpec:
    """The (Temperature x Precipitation) stress-test grid."""

    temp_min: float
    temp_max: float
    temp_increment: float
    precip_min: float
    precip_max: float
    precip_increment: float

    def temp_axis(self) -> np.ndarray:
        return r_seq(self.temp_min, self.temp_max, self.temp_increment)

    def precip_axis(self) -> np.ndarray:
        return r_seq(self.precip_min, self.precip_max, self.precip_increment)

    def grid(self) -> pd.DataFrame:
        """Grid points as R ``expand.grid(Temp, Precip)`` produces them: **Temp varies fastest**.

        This row ordering is load-bearing -- the biv_norm_vals output and the cumulative-area
        contour transform both depend on it.
        """
        temp = self.temp_axis()
        precip = self.precip_axis()
        return pd.DataFrame(
            {
                "T_lev": np.tile(temp, len(precip)),
                "P_lev": np.repeat(precip, len(temp)),
            }
        )


@dataclass(frozen=True)
class DomainConfig:
    """Everything needed to run one domain end-to-end."""

    name: str
    source_kind: str  # 'loca2-flow' | 'loca2-basin' | 'cmip5'
    base_center: int
    window: int
    yr_start: int
    filename_field_order: tuple[str, ...]
    time_column: str
    pr_type: str  # 'D_pr_lm' | 'D_pr'
    gcm_filter: tuple[str, ...]
    filter_gcms: bool
    filter_nmem_gcms: bool
    grid: GridSpec
    prob_levels: tuple[float, ...]
    prob_interval_count: int
    projection_change_years: tuple[int, ...]
    plot_periods: tuple[int, ...]
    animation_periods: tuple[int, int] | None
    ssp_colors: dict
    paths: dict  # resolved absolute Paths (+ any passthrough values)
    outputs: dict
    basin_id: str | None = None
    extras: dict = field(default_factory=dict)
    plot_grid: GridSpec | None = None  # finer grid for rendered figures only; None -> use `grid`

    # --- derived artifact paths (single source of truth for processed/ filenames) ---
    @property
    def lmfits_path(self) -> Path:
        return self.paths["processed_dir"] / f"loca2_lmfits_{self.yr_start}_{self.name}.csv"

    @property
    def warming_levels_path(self) -> Path:
        return self.paths["processed_dir"] / f"loca2_warming_levels_{self.name}.csv"

    @property
    def thirtyyr_avgs_path(self) -> Path:
        return self.paths["processed_dir"] / f"loca2_30yravgs_{self.name}.csv"

    @property
    def gcm_mean_path(self) -> Path:
        return self.paths["processed_dir"] / f"gcm_mean_loca2_varavg_lm_{self.name}.csv"

    @property
    def gcm_sigs_path(self) -> Path:
        return self.paths["processed_dir"] / f"gcm_sigs_loca2_varavg_lm_{self.name}.csv"

    @property
    def biv_norm_vals_path(self) -> Path:
        return self.paths["processed_dir"] / f"biv_norm_vals_dt-dp_loca2_{self.name}.csv"

    @property
    def gcm_sigs_decomposed_path(self) -> Path:
        return self.paths["processed_dir"] / f"gcm_sigs_loca2_varavg_lm_{self.name}_decomposed.csv"

    @property
    def gcm_points_path(self) -> Path:
        return self.paths["processed_dir"] / f"gcm_points_loca2_varavg_lm_{self.name}.csv"


def _parse_grid(g: dict) -> GridSpec:
    return GridSpec(
        temp_min=g["temp"]["min"],
        temp_max=g["temp"]["max"],
        temp_increment=g["temp"]["increment"],
        precip_min=g["precip"]["min"],
        precip_max=g["precip"]["max"],
        precip_increment=g["precip"]["increment"],
    )


def _resolve_path(value: str, repo_root: Path, basin_id: str | None) -> Path:
    if basin_id is not None:
        value = value.format(basin_id=basin_id)
    p = Path(value)
    return p if p.is_absolute() else (repo_root / p)


def load_domain(
    name: str,
    configs_dir: Path | None = None,
    repo_root: Path | None = None,
) -> DomainConfig:
    """Load and validate ``configs/<name>.yaml`` into a ``DomainConfig`` with resolved paths."""
    configs_dir = Path(configs_dir) if configs_dir else CONFIGS_DIR
    repo_root = Path(repo_root) if repo_root else REPO_ROOT

    cfg_path = configs_dir / f"{name}.yaml"
    if not cfg_path.exists():
        available = sorted(p.stem for p in configs_dir.glob("*.yaml"))
        raise FileNotFoundError(f"No config '{name}' in {configs_dir}. Available: {available}")

    with open(cfg_path, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)

    grid = _parse_grid(raw["grid"])
    plot_grid = _parse_grid(raw["plot_grid"]) if "plot_grid" in raw else None

    basin_id = raw.get("basin_id")
    paths = {
        k: (_resolve_path(v, repo_root, basin_id) if isinstance(v, str) else v)
        for k, v in raw.get("paths", {}).items()
    }

    anim = raw.get("animation_periods")
    animation_periods = (anim["start"], anim["stop"]) if anim else None

    known = {
        "name", "source_kind", "base_center", "window", "yr_start", "filename_field_order",
        "time_column", "pr_type", "gcm_filter", "filter_gcms", "filter_nmem_gcms", "grid",
        "plot_grid", "prob_levels", "prob_interval_count", "projection_change_years",
        "plot_periods", "animation_periods", "ssp_colors", "paths", "outputs", "basin_id",
    }
    extras = {k: v for k, v in raw.items() if k not in known}

    return DomainConfig(
        name=raw["name"],
        source_kind=raw["source_kind"],
        base_center=raw["base_center"],
        window=raw["window"],
        yr_start=raw["yr_start"],
        filename_field_order=tuple(raw["filename_field_order"]),
        time_column=raw["time_column"],
        pr_type=raw["pr_type"],
        gcm_filter=tuple(raw.get("gcm_filter", [])),
        filter_gcms=raw.get("filter_gcms", False),
        filter_nmem_gcms=raw.get("filter_nmem_gcms", False),
        grid=grid,
        prob_levels=tuple(raw.get("prob_levels", [0.68, 0.95])),
        prob_interval_count=raw.get("prob_interval_count", 100),
        projection_change_years=tuple(raw.get("projection_change_years", [])),
        plot_periods=tuple(raw.get("plot_periods", [])),
        animation_periods=animation_periods,
        ssp_colors=raw.get("ssp_colors", {}),
        paths=paths,
        outputs=raw.get("outputs", {}),
        basin_id=basin_id,
        extras=extras,
        plot_grid=plot_grid,
    )
