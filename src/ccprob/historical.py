"""Canonical historical daily temperature/precipitation, spatially aggregated to CalSim basins.

Standalone analysis (not part of the LOCA-2/CMIP5 bivariate-normal pipeline): assigns each cell of
the WGEN statewide daily gridded dataset to a CalSim basin region (see ``calsim_basins.py``) by
point-in-polygon membership of the cell's center, averages daily precip/tmax/tmin across a basin's
member cells, and combines the per-basin series into one Flow-Ratio-weighted valley-wide series --
the same weighting the existing ``cv-flow-weighted`` LOCA-2 product uses, renormalized over the
basins with a resolvable polygon (excludes Goose Lake, flow ratio 0 anyway, and the Delta, which
has no polygon in the geopackage).

The WGEN grid (~13,800 files, ~13 GB) lives outside this repo and is never copied in; only the much
smaller per-basin daily CSVs this module produces are written into ``data/historical/``. The basin
registry + point-in-polygon + flow-weighting steps (``_run_basin_aggregation``) are shared with the
sibling LOC95 grid ingestion (``loc95.py``) -- only the grid-cell lister and output filename prefix
differ between the two.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from shapely.geometry import Point
from shapely.geometry.base import BaseGeometry

from .calsim_basins import load_basin_registry

_CELL_RE = re.compile(r"data_(-?[0-9.]+)_(-?[0-9.]+)$")
WGEN_COLUMNS = ["year", "month", "day", "pr", "tmax", "tmin"]

OUTPUT_COLUMNS = {
    "year": "Year",
    "month": "Month",
    "day": "Day",
    "pr": "Pr (mm)",
    "tmax": "Tasmax (degC)",
    "tmin": "Tasmin (degC)",
}


def list_wgen_cells(
    wgen_dir: Path, pattern: str = "data_*", cell_re: re.Pattern = _CELL_RE
) -> list[tuple[float, float, Path]]:
    """Parse ``(lat, lon, path)`` for every WGEN-style grid-cell file matching ``pattern``, without
    opening any of them. ``pattern``/``cell_re`` are overridable so a sibling grid with a different
    filename prefix (e.g. LOC95's ``meteo_<lat>_<lon>``, see ``loc95.py``) can reuse this lister."""
    cells = []
    for path in Path(wgen_dir).glob(pattern):
        m = cell_re.match(path.name)
        if m:
            cells.append((float(m.group(1)), float(m.group(2)), path))
    return cells


def assign_cells_to_basins(
    cells: list[tuple[float, float, Path]], basin_polygons: dict[int, BaseGeometry]
) -> dict[int, list[Path]]:
    """Which WGEN cell files fall inside each basin's dissolved polygon (by cell-center point)."""
    assignment: dict[int, list[Path]] = {region_id: [] for region_id in basin_polygons}
    for lat, lon, path in cells:
        point = Point(lon, lat)  # shapely/WKB convention: x=lon, y=lat
        for region_id, polygon in basin_polygons.items():
            if polygon.contains(point):
                assignment[region_id].append(path)
    return assignment


def read_wgen_daily(path: Path) -> pd.DataFrame:
    """One WGEN grid-cell file: whitespace-delimited, no header."""
    return pd.read_csv(path, sep=r"\s+", header=None, names=WGEN_COLUMNS)


def aggregate_basin_daily(paths: list[Path]) -> pd.DataFrame:
    """Mean daily pr/tmax/tmin across a basin's member grid cells, one row per calendar day."""
    frames = [read_wgen_daily(p) for p in paths]
    combined = pd.concat(frames, axis=0, ignore_index=True)
    return combined.groupby(["year", "month", "day"], as_index=False)[["pr", "tmax", "tmin"]].mean()


def flow_weighted_aggregate(
    basin_frames: dict[int, pd.DataFrame], flow_ratios: dict[int, float]
) -> pd.DataFrame:
    """Combine per-basin daily frames into one Flow-Ratio-weighted series.

    Renormalizes ``flow_ratios`` over exactly the basins present in ``basin_frames`` (the caller
    only passes basins with a resolvable polygon), so a dropped region's weight is redistributed
    proportionally rather than silently lost.
    """
    total_ratio = sum(flow_ratios[region_id] for region_id in basin_frames)
    weights = {region_id: flow_ratios[region_id] / total_ratio for region_id in basin_frames}

    weighted = None
    for region_id, frame in basin_frames.items():
        contribution = frame.set_index(["year", "month", "day"])[["pr", "tmax", "tmin"]] * weights[region_id]
        weighted = contribution if weighted is None else weighted.add(contribution, fill_value=0)
    return weighted.reset_index()


def _write_daily_csv(frame: pd.DataFrame, out_path: Path) -> Path:
    out = frame.sort_values(["year", "month", "day"]).rename(columns=OUTPUT_COLUMNS)
    out = out[list(OUTPUT_COLUMNS.values())]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    return out_path


def _run_basin_aggregation(
    cells: list[tuple[float, float, Path]],
    gpkg_path: Path,
    xlsx_path: Path,
    out_dir: Path,
    file_prefix: str,
    basins: list[int] | None = None,
) -> dict:
    """Shared basin-registry + point-in-polygon + flow-weighting pipeline: assign ``cells`` to
    basins, aggregate, and write ``<file_prefix>_daily_basin-NN.csv`` / ``<file_prefix>_daily_cv-
    flow-weighted.csv``. Used by both ``run`` (WGEN historical) and ``loc95.run`` (LOC95) -- only
    the cell list and filename prefix differ."""
    registry = load_basin_registry(gpkg_path, xlsx_path)
    resolvable = {rid: r for rid, r in registry.items() if r.polygon is not None}
    if basins is not None:
        resolvable = {rid: r for rid, r in resolvable.items() if rid in basins}

    assignment = assign_cells_to_basins(cells, {rid: r.polygon for rid, r in resolvable.items()})

    out_dir = Path(out_dir)
    written: dict = {"basins": {}, "cv_flow_weighted": None}
    basin_frames: dict[int, pd.DataFrame] = {}
    for region_id in sorted(resolvable):
        paths = assignment[region_id]
        if not paths:
            continue
        frame = aggregate_basin_daily(paths)
        basin_frames[region_id] = frame
        out_path = out_dir / f"{file_prefix}_daily_basin-{region_id:02d}.csv"
        written["basins"][region_id] = _write_daily_csv(frame, out_path)

    if basins is None and basin_frames:
        flow_ratios = {rid: r.flow_ratio for rid, r in resolvable.items()}
        cv_frame = flow_weighted_aggregate(basin_frames, flow_ratios)
        written["cv_flow_weighted"] = _write_daily_csv(
            cv_frame, out_dir / f"{file_prefix}_daily_cv-flow-weighted.csv"
        )

    return written


def run(
    wgen_dir: Path,
    gpkg_path: Path,
    xlsx_path: Path,
    out_dir: Path,
    basins: list[int] | None = None,
) -> dict:
    """Build the basin registry, aggregate WGEN cells per basin, write the historical CSVs."""
    cells = list_wgen_cells(wgen_dir)
    return _run_basin_aggregation(cells, gpkg_path, xlsx_path, out_dir, "historical", basins=basins)
