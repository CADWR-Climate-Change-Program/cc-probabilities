"""Precipitation trends and warming levels (the former Stage-1 Python preprocessing).

One regression implementation (scipy ``linregress`` -- a superset of the legacy scipy/sklearn pair).
The flow-vs-basin differences are parameterized off ``source_kind`` rather than duplicated:

==================  =====================  ==========================
                    flow (loca2-flow)      basin (loca2-basin)
==================  =====================  ==========================
regression target   30-yr rolling mean pr  raw pr (already 30-yr avg)
post-load filter     y >= yr_start+window   y >= yr_start
projection x-offset  + window               + 0
warming-year offset  + 0                    + (window-1)
warming diff source  tavg_roll              tavg
emits r/p/se cols    yes                    no (legacy sklearn gave none)
==================  =====================  ==========================
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import linregress
from scipy.stats import t as student_t


@dataclass
class RegressionResult:
    slope: float
    intercept: float
    rvalue: float
    pvalue: float
    stderr: float
    intercept_stderr: float


def linear_regression(x, y) -> RegressionResult:
    """OLS fit of y on x via scipy.linregress (the single canonical implementation)."""
    r = linregress(x, y)
    return RegressionResult(
        slope=r.slope,
        intercept=r.intercept,
        rvalue=r.rvalue,
        pvalue=r.pvalue,
        stderr=r.stderr,
        intercept_stderr=r.intercept_stderr,
    )


def _tinv(p: float, df: int) -> float:
    """Two-sided inverse Student-t (legacy `tinv`)."""
    return abs(student_t.ppf(p / 2, df))


class TrendModel:
    def __init__(self, cfg):
        self.cfg = cfg
        self.rolling = cfg.source_kind == "loca2-flow"
        self.window = cfg.window

    # -- shared prep: rolling means (flow) + the post-load year filter --
    def prepare(self, ensemble) -> pd.DataFrame:
        df = ensemble.long_frame().copy()
        if self.rolling:
            df["pr_roll"] = df.groupby("mvs")["pr"].transform(
                lambda s: s.rolling(self.window).mean()
            )
            df["tavg_roll"] = df.groupby("mvs")["tavg"].transform(
                lambda s: s.rolling(self.window).mean()
            )
            df = df.dropna(axis=0, how="any")
            df = df[df["y"] >= self.cfg.yr_start + self.window]
        else:
            df = df[df["y"] >= self.cfg.yr_start]
        return df.reset_index(drop=True)

    # -- linear-model fits + projected % change at the target years --
    def lm_fits(self, prepared: pd.DataFrame) -> pd.DataFrame:
        cfg = self.cfg
        pr_col = "pr_roll" if self.rolling else "pr"
        offset = self.window if self.rolling else 0
        x_base = (cfg.base_center - 14) + offset

        rows = []
        for mvs, data in prepared.groupby("mvs", sort=True):
            reg = linear_regression(data["y"], data[pr_col])
            model, variant, ssp = mvs.split("_")
            base = reg.slope * x_base + reg.intercept

            row: dict = {"mvs": mvs, "m": model, "v": variant, "s": ssp,
                         "slope": reg.slope, "intercept": reg.intercept}
            if self.rolling:  # flow keeps the regression statistics; basin (legacy sklearn) did not
                ts = _tinv(0.05, len(data["y"]) - 2)
                row["rvalue"] = reg.rvalue      # legacy column was misnamed 'r2' (it holds rvalue, not r-squared)
                row["pvalue"] = reg.pvalue
                row["slope95"] = ts * reg.stderr
                row["intercept95"] = ts * reg.intercept_stderr
            row["base"] = base

            proj = {}
            for yr in cfg.projection_change_years:
                x = (yr - 15) + offset
                proj[yr] = reg.slope * x + reg.intercept
                row[str(yr)] = proj[yr]

            for yr in cfg.projection_change_years:
                row[f"pr_change_{yr}"] = round((proj[yr] - base) / base * 100, 1)
            rows.append(row)

        fits = pd.DataFrame(rows).sort_values(["m", "s"]).reset_index(drop=True)
        return fits

    # -- year each mvs first crosses +1..+5 C of warming --
    def warming_levels(self, prepared: pd.DataFrame, levels=(1, 2, 3, 4, 5)) -> pd.DataFrame:
        diff_col = "tavg_roll" if self.rolling else "tavg"
        year_offset = 0 if self.rolling else (self.window - 1)

        result = {}
        for mvs, data in prepared.groupby("mvs", sort=True):
            d = data.reset_index(drop=True)
            diffs = d[diff_col].to_numpy() - d[diff_col].to_numpy()[0]
            years = d["y"].to_numpy()
            result[mvs] = {
                t: int(years[np.searchsorted(diffs, t) - 1] + year_offset) for t in levels
            }

        wl = pd.DataFrame.from_dict(result, orient="index")
        wl.index.name = "mvs"
        wl = wl.reset_index()
        wl.columns = ["mvs", *levels]
        return wl
