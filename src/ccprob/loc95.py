"""LOC95 daily temperature/precipitation, spatially aggregated to CalSim basins.

Sibling to ``historical.py``: same point-in-polygon basin assignment and Flow-Ratio-weighted
valley-wide combination (``_run_basin_aggregation``), applied instead to the LOC95
stochastic-weather-generator grid -- ``meteo_<lat>_<lon>`` files, whitespace-delimited
year/month/day/pr/tmax/tmin with no header, the same six-column layout as the WGEN historical grid
(see ``historical.read_wgen_daily``), just a different filename prefix and an external, 1.8%-dried
precipitation scenario (see ``precip_extremes.LOC95_DRY_PCT``).

The LOC95 grid lives outside this repo and is never copied in; only the per-basin daily CSVs this
module produces are written into ``data/loc/``.
"""

from __future__ import annotations

import re
from pathlib import Path

from .historical import _run_basin_aggregation, list_wgen_cells

_LOC95_CELL_RE = re.compile(r"meteo_(-?[0-9.]+)_(-?[0-9.]+)$")


def list_loc95_cells(loc_dir: Path) -> list[tuple[float, float, Path]]:
    """Parse ``(lat, lon, path)`` for every LOC95 grid-cell file (``meteo_<lat>_<lon>``)."""
    return list_wgen_cells(loc_dir, pattern="meteo_*", cell_re=_LOC95_CELL_RE)


def run(
    loc_dir: Path,
    gpkg_path: Path,
    xlsx_path: Path,
    out_dir: Path,
    basins: list[int] | None = None,
) -> dict:
    """Build the basin registry, aggregate LOC95 cells per basin, write the loc95 CSVs."""
    cells = list_loc95_cells(loc_dir)
    return _run_basin_aggregation(cells, gpkg_path, xlsx_path, out_dir, "loc95", basins=basins)
