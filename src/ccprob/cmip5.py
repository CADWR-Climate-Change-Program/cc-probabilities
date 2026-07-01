"""The self-contained CMIP5 path (legacy ``createBivNormValues.R``).

CMIP5 differs from the LOCA-2 flow pipeline enough to warrant its own module rather than threading
flags through the LOCA-2 classes:

- No Python preprocessing / linear-model precip: ``D_pr`` is the raw 30-yr-window precip change.
- Two emission scenarios (rcp45 / rcp85), each with its OWN historical baseline workbook.
- Ensemble members are collapsed by **GCM family** (``data/gcm-families.csv``), not by realization.
- Per-year mean and 2x2 covariance are taken across the family points (both RCPs pooled).

It reuses ``BivariateNormalSurface`` (grid + bivariate normal) and the ``_rcompat`` writers; only the
data assembly and the per-period label (a ``FUT(YYYY-YYYY)`` window string) are CMIP5-specific.
"""

from __future__ import annotations

import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd

from ._rcompat import write_gcm_mean_csv, write_gcm_sigs_csv
from .distribution import BivariateNormalSurface
from .io import WorksheetReader
from .plotting import ContourPlotter

# CMIP5 period numbering counts windows from 1995 (window k -> period 1995 + k), independent of
# base_center; the first 30-yr window FUT(1982-2011) is k = 1 -> period 1996.
PERIOD_ORIGIN = 1995


class Cmip5Ensemble:
    """Load the CMIP5 climatology workbooks and assemble family-collapsed (ΔT, ΔP) points."""

    def __init__(self, cfg):
        self.cfg = cfg

    def _worksheets(self) -> list:
        return sorted(glob.glob(str(self.cfg.paths["input_dir"] / "*.xlsx")))

    def _families(self) -> pd.DataFrame:
        fam = pd.read_csv(self.cfg.paths["families_csv"])
        fam.columns = [c.strip() for c in fam.columns]
        fam = fam.rename(columns={"GCM Family": "GCM_Family"})
        fam["GCM"] = fam["GCM"].str.upper()
        return fam[["GCM", "GCM_Family"]]

    @staticmethod
    def _name_rcp(path) -> tuple[str, str]:
        # CaliforniaCentralValley_CMIP5_Climatology_Hist(..)_FUT(YYYY-YYYY)_rcpNN.xlsx
        parts = os.path.basename(str(path)).split("_")
        return parts[4], parts[5].replace(".xlsx", "")

    def build_deltas(self) -> pd.DataFrame:
        """Per-(GCM, RCP, window) raw deltas vs the RCP-matched 2006-centered baseline workbook."""
        cfg = self.cfg
        files = self._worksheets()
        # 1-based index of the rcp85 baseline file == the window centered on base_center.
        fi = (cfg.base_center - PERIOD_ORIGIN) * 2
        hist = {
            "rcp45": WorksheetReader.cmip5_totals(files[fi - 2]),
            "rcp85": WorksheetReader.cmip5_totals(files[fi - 1]),
        }

        pieces = []
        for i, path in enumerate(files, start=1):
            name, rcp = self._name_rcp(path)
            cur = WorksheetReader.cmip5_totals(path)
            base = hist[rcp]
            d = pd.DataFrame({"GCM": cur["GCM"]})
            d["D_pr"] = (cur["Pr_total"] - base["Pr_total"]) / base["Pr_total"] * 100.0
            d["D_tas"] = cur["Tas_mean"] - base["Tas_mean"]
            d["RCP"] = rcp
            d["Year"] = name
            d["period"] = PERIOD_ORIGIN + (i // 2 if i % 2 == 0 else (i - 1) // 2 + 1)
            pieces.append(d)

        gcm = pd.concat(pieces, axis=0, ignore_index=True)
        gcm = gcm.dropna(subset=["GCM", "D_pr", "D_tas"])  # R complete.cases
        gcm["GCM"] = gcm["GCM"].str.upper()
        return gcm.reset_index(drop=True)

    def collapse_families(self, gcm: pd.DataFrame) -> pd.DataFrame:
        """Average member GCMs within each (family, RCP, window) -> columns GCM_Family, RCP, Year, DT, DP.

        Left-joins the family table and keeps unmatched GCMs as an NaN-family group (``dropna=False``),
        matching R's ``left_join`` + ``group_by`` which does not drop NA families.
        """
        merged = gcm.merge(self._families(), on="GCM", how="left")
        return (
            merged.groupby(["GCM_Family", "RCP", "Year"], dropna=False, as_index=False)
            .agg(DT=("D_tas", "mean"), DP=("D_pr", "mean"))
        )


def period_means(family: pd.DataFrame) -> pd.DataFrame:
    """Per-window ensemble mean of (DT, DP) across all family x RCP points -> Year, DT_mean, DP_mean."""
    gm = family.groupby("Year", as_index=False)[["DT", "DP"]].mean()
    gm = gm.rename(columns={"DT": "DT_mean", "DP": "DP_mean"})
    return gm.sort_values("Year").reset_index(drop=True)


def period_covariances(family: pd.DataFrame) -> dict:
    """Per window, the 2x2 covariance of (DT, DP) (ddof=1), keyed by the Year window string."""
    out = {}
    for year, d in family.groupby("Year"):
        out[year] = np.cov(d["DT"].to_numpy(), d["DP"].to_numpy(), ddof=1)
    return out


class Cmip5Pipeline:
    """Orchestrate the CMIP5 path: deltas -> family collapse -> mean/cov -> surfaces / CSVs / figures."""

    def __init__(self, cfg):
        self.cfg = cfg

    def _compute(self) -> dict:
        ens = Cmip5Ensemble(self.cfg)
        deltas = ens.build_deltas()
        family = ens.collapse_families(deltas)
        gm = period_means(family)
        covs = period_covariances(family)
        means = gm.rename(columns={"Year": "period", "DT_mean": "DT", "DP_mean": "DP"})
        frames = BivariateNormalSurface(self.cfg).all_periods(means, covs)
        return {
            "deltas": deltas,
            "family": family,
            "gcm_mean": gm,
            "covariances": covs,
            "biv_norm_frames": frames,
        }

    def _gcm_mean_path(self) -> Path:
        return self.cfg.paths["processed_dir"] / "gcm_mean_cmip5_family.csv"

    def _gcm_sigs_path(self) -> Path:
        return self.cfg.paths["processed_dir"] / "gcm_sigs_cmip5_family.csv"

    def run_distribution(self, out_dir=None, write: bool = True) -> dict:
        """Compute and (optionally) write gcm_mean / gcm_sigs.

        ``biv_norm_vals`` for CMIP5 is a ~100 MB wide CSV on the fine 0.1 deg grid; it is computed
        in-memory (``biv_norm_frames``) but only written when ``outputs.biv_norm_vals`` is enabled.
        """
        cfg = self.cfg
        result = self._compute()
        if write:
            out = Path(out_dir) if out_dir else cfg.paths["processed_dir"]
            out.mkdir(parents=True, exist_ok=True)
            write_gcm_mean_csv(result["gcm_mean"], out / self._gcm_mean_path().name)
            write_gcm_sigs_csv(
                result["covariances"], out / self._gcm_sigs_path().name,
                row_labels=("DT", "DP"),
            )
            if cfg.outputs.get("biv_norm_vals"):
                from ._rcompat import write_biv_norm_vals_csv
                write_biv_norm_vals_csv(
                    result["biv_norm_frames"], out / "biv_norm_vals_dt-dp_cmip5.csv"
                )
        return result

    def run_plots(self, out_dir=None, write: bool = True) -> list:
        """Render the contour figure for each ``plot_periods`` offset (period = 1995 + offset)."""
        cfg = self.cfg
        result = self._compute()
        family, frames = result["family"], result["biv_norm_frames"]
        frames_by_year = {f["period"].iloc[0]: f for f in frames}
        period_to_year = dict(zip(result["deltas"]["period"], result["deltas"]["Year"]))

        plotter = ContourPlotter(cfg)
        out = Path(out_dir) if out_dir else cfg.paths["figures_dir"]
        written = []
        if write:
            out.mkdir(parents=True, exist_ok=True)
        for offset in cfg.plot_periods:
            year = period_to_year[PERIOD_ORIGIN + offset]
            frame = frames_by_year[year]
            points = family[family["Year"] == year]
            if write:
                name = f"biv-norm-prob-dt-dp-cmip5-{PERIOD_ORIGIN + offset}.svg"
                mean_row = result["gcm_mean"].loc[result["gcm_mean"]["Year"] == year]
                mean = (float(mean_row["DT_mean"].iloc[0]), float(mean_row["DP_mean"].iloc[0]))
                plotter.filled_contour(
                    frame, points, title=plotter.title_for(PERIOD_ORIGIN + offset),
                    out_path=out / name, mean=mean, cov=result["covariances"][year],
                )
                written.append(out / name)
        return written

    def run(self, out_dir=None, write: bool = True, plots: bool = True) -> dict:
        dist = self.run_distribution(out_dir=out_dir, write=write)
        figures = self.run_plots(out_dir=out_dir, write=write) if plots else []
        return {**dist, "figures": figures}
