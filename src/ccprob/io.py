"""Filename parsing and raw-projection readers.

One canonical ``FilenameParser`` (explicit per-source field order, configured by the domain) and
typed readers for the annual-flow and 30-yr-basin projection CSVs. Replaces the three near-identical
filename-parse + load loops in the legacy scripts. Worksheet (xlsx) and CSV writers are added in
later phases.
"""

from __future__ import annotations

import os

import openpyxl
import pandas as pd


class FilenameParser:
    """Parse model/variant/ssp out of a projection filename by configured field order.

    The flow and basin input directories order these fields differently, so the order is supplied
    by ``DomainConfig.filename_field_order`` rather than hard-coded. ``mvs`` is always the canonical
    ``model_variant_ssp`` regardless of on-disk order.
    """

    def __init__(self, field_order):
        self.field_order = list(field_order)

    def parse(self, path) -> dict:
        parts = os.path.basename(str(path)).split("_")
        return {field: parts[i] for i, field in enumerate(self.field_order)}

    @staticmethod
    def mvs(model: str, variant: str, ssp: str) -> str:
        return f"{model}_{variant}_{ssp}"


class ProjectionReader:
    """Read one raw projection file into a tidy per-year frame (columns: y, pr, tavg)."""

    def __init__(self, cfg):
        self.cfg = cfg

    def _aggregate(self, df: pd.DataFrame, time_column: str) -> pd.DataFrame:
        # monthly rows -> per-(window-)year totals: sum precip, mean temperature (matches legacy)
        agg = df.groupby(time_column, as_index=False).aggregate(
            {"Pr (mm)": "sum", "Tave (degC)": "mean"}
        )
        agg.columns = ["y", "pr", "tavg"]
        return agg

    def read_annual(self, path) -> pd.DataFrame:
        """Flow domain: monthly annual-flow CSV -> per-Year totals."""
        return self._aggregate(pd.read_csv(path), "Year")

    def read_30y(self, path) -> pd.DataFrame:
        """Basin domain: pre-averaged 30-yr CSV -> per-(30y start year) totals."""
        return self._aggregate(pd.read_csv(path), self.cfg.time_column)


class WorksheetReader:
    """Read the LOCA-2 summary worksheets (xlsx): per-period 30-yr monthly P/T climatologies.

    Each workbook has sheets pr_hist/tas_hist/pr_fut/tas_fut; columns 1-3 are GCM/Realization/SSP
    and columns 4-15 are the 12 monthly values. Precip total = row sum; temperature = row mean
    (matching the R rowSums/rowMeans).
    """

    @staticmethod
    def _sheet(path, sheet: str) -> pd.DataFrame:
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        try:
            rows = list(wb[sheet].iter_rows(values_only=True))
        finally:
            wb.close()
        return pd.DataFrame(rows[1:], columns=list(rows[0]))

    @classmethod
    def _totals(cls, path, pr_sheet: str, tas_sheet: str) -> pd.DataFrame:
        pr = cls._sheet(path, pr_sheet)
        tas = cls._sheet(path, tas_sheet)
        out = pr[["GCM", "Realization", "SSP"]].copy()
        out["Pr_total"] = pr[list(pr.columns[3:15])].sum(axis=1)
        out["Tas_mean"] = tas[list(tas.columns[3:15])].mean(axis=1)
        return out.reset_index(drop=True)

    @classmethod
    def hist_totals(cls, path) -> pd.DataFrame:
        return cls._totals(path, "pr_hist", "tas_hist")

    @classmethod
    def fut_totals(cls, path) -> pd.DataFrame:
        return cls._totals(path, "pr_fut", "tas_fut")

    @classmethod
    def cmip5_totals(cls, path) -> pd.DataFrame:
        """CMIP5 worksheets: GCM in column 0, the 12 monthly values in columns 1-12 of pr_fut/tas_fut.

        (No Realization/SSP columns -- the CMIP5 climatology layout differs from LOCA-2.) Precip
        total = row sum, temperature = row mean, matching the R ``rowSums(.[-1])`` / ``rowMeans``.
        """
        pr = cls._sheet(path, "pr_fut")
        tas = cls._sheet(path, "tas_fut")
        out = pd.DataFrame({"GCM": pr[pr.columns[0]]})
        out["Pr_total"] = pr[list(pr.columns[1:13])].sum(axis=1)
        out["Tas_mean"] = tas[list(tas.columns[1:13])].mean(axis=1)
        return out.reset_index(drop=True)
