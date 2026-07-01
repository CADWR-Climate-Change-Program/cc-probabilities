"""Climate deltas and the per-period bivariate-normal surface (the former Stage-2 R stats).

``ClimateDeltas`` turns the summary worksheets + lm fits into per-(mvs, period) (ΔT, ΔP) deltas,
collapses model realizations, and yields the per-period ensemble mean and 2x2 covariance.
``BivariateNormalSurface`` (added next) evaluates the normalized bivariate normal over the grid.
"""

from __future__ import annotations

import glob

import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal

from .io import WorksheetReader
from .trends import linear_regression


class ClimateDeltas:
    """Build per-period (ΔT, ΔP) deltas + lm fits into the per-(mvs, period) delta table.

    ``lm_fits`` is the lmfits table (columns m, v, s, slope, intercept, ...) -- passed in rather
    than recomputed, mirroring the R scripts which read the lmfits CSV. The source-specific load
    differs (flow reads the xlsx worksheets; basin reads the Python-written 30yravgs table, so
    ``thirtyyr_avgs`` is required for ``loca2-basin``) but both converge on the same column set
    (Model, Variant, SSP, Pr_total, Tas_mean, D_pr, D_pr_lm, D_tas, period) before
    ``collapse_variants`` / ``period_means`` / ``period_covariances``, which are source-agnostic.
    """

    def __init__(self, cfg, lm_fits: pd.DataFrame, thirtyyr_avgs: pd.DataFrame | None = None):
        self.cfg = cfg
        self.lm_fits = lm_fits.copy()
        self.thirtyyr_avgs = thirtyyr_avgs

    def _worksheets(self):
        return sorted(glob.glob(str(self.cfg.paths["worksheets_dir"] / "*.xlsx")))

    def build(self) -> pd.DataFrame:
        cfg = self.cfg
        gcm = self._build_basin() if cfg.source_kind == "loca2-basin" else self._build_flow()
        if cfg.filter_gcms:  # drop filtered models BEFORE uppercasing (matches R order)
            gcm = gcm[~gcm["Model"].isin(cfg.gcm_filter)]
        gcm["Model"] = gcm["Model"].str.upper()
        return gcm.reset_index(drop=True)

    def _build_flow(self) -> pd.DataFrame:
        cfg = self.cfg
        files = self._worksheets()
        hist = WorksheetReader.hist_totals(files[0])  # GCM/Realization/SSP/Pr_total/Tas_mean

        lm = self.lm_fits
        lm_pr_hist = lm["slope"] * (cfg.base_center - 14) + lm["intercept"]

        hist_rows = hist.rename(columns={"GCM": "Model", "Realization": "Variant"}).copy()
        hist_rows["D_pr"] = 0.0
        hist_rows["D_pr_lm"] = 0.0
        hist_rows["D_tas"] = 0.0
        hist_rows["period"] = cfg.base_center
        pieces = [hist_rows]

        for i, path in enumerate(files, start=1):
            fut = WorksheetReader.fut_totals(path)
            deltas = fut.rename(columns={"GCM": "Model", "Realization": "Variant"}).copy()
            # element-wise vs the positionally-aligned historical baseline (same row order in xlsx)
            deltas["D_pr"] = (fut["Pr_total"] - hist["Pr_total"]) / hist["Pr_total"] * 100.0
            deltas["D_tas"] = fut["Tas_mean"] - hist["Tas_mean"]
            deltas["period"] = cfg.base_center + i

            lm_pr_fut = lm["slope"] * (cfg.base_center + i - 14) + lm["intercept"]
            lm_dprlm = pd.DataFrame(
                {
                    "Model": lm["m"], "Variant": lm["v"], "SSP": lm["s"],
                    "D_pr_lm": (lm_pr_fut - lm_pr_hist) / lm_pr_hist * 100.0,
                }
            )
            deltas = deltas.merge(lm_dprlm, on=["Model", "Variant", "SSP"], how="inner")
            pieces.append(deltas)

        return pd.concat(pieces, axis=0, ignore_index=True)

    def _build_basin(self) -> pd.DataFrame:
        """Basin equivalent of ``_build_flow``, sourced from the Python-written 30yravgs table.

        The R script shifts the 30yravgs window-start year by +14 (the same window-center anchor
        ``lm_fits`` uses via ``base_center - 14``) so ``y`` lines up with ``period``, then takes
        positionally-aligned vectors per ``subset(y==...)``. We merge on ``mvs`` instead of relying
        on row order -- numerically identical here (one row per mvs per year, no gaps) and robust if
        that ever stops holding.
        """
        cfg = self.cfg
        avgs = self.thirtyyr_avgs[["mvs", "model", "variant", "ssp", "y", "pr", "tavg"]].copy()
        avgs["y"] = avgs["y"] + 14
        avgs = avgs[avgs["y"] >= cfg.base_center]
        avgs = avgs.rename(
            columns={
                "model": "Model", "variant": "Variant", "ssp": "SSP",
                "pr": "Pr_total", "tavg": "Tas_mean",
            }
        )

        hist = avgs[avgs["y"] == cfg.base_center]
        hist_join = hist[["mvs", "Pr_total", "Tas_mean"]]

        hist_rows = hist[["Model", "Variant", "SSP", "Pr_total", "Tas_mean"]].copy()
        hist_rows["D_pr"] = 0.0
        hist_rows["D_pr_lm"] = 0.0
        hist_rows["D_tas"] = 0.0
        hist_rows["period"] = cfg.base_center
        pieces = [hist_rows]

        lm = self.lm_fits
        lm_pr_hist = lm["slope"] * (cfg.base_center - 14) + lm["intercept"]

        for period in sorted(y for y in avgs["y"].unique() if y != cfg.base_center):
            fut = avgs[avgs["y"] == period].merge(hist_join, on="mvs", suffixes=("", "_hist"))
            fut["D_pr"] = (fut["Pr_total"] - fut["Pr_total_hist"]) / fut["Pr_total_hist"] * 100.0
            fut["D_tas"] = fut["Tas_mean"] - fut["Tas_mean_hist"]
            fut["period"] = period

            lm_pr_fut = lm["slope"] * (period - 14) + lm["intercept"]
            lm_dprlm = pd.DataFrame(
                {
                    "Model": lm["m"], "Variant": lm["v"], "SSP": lm["s"],
                    "D_pr_lm": (lm_pr_fut - lm_pr_hist) / lm_pr_hist * 100.0,
                }
            )
            fut = fut.merge(lm_dprlm, on=["Model", "Variant", "SSP"], how="inner")
            pieces.append(
                fut[["Model", "Variant", "SSP", "Pr_total", "Tas_mean", "D_pr", "D_pr_lm", "D_tas", "period"]]
            )

        return pd.concat(pieces, axis=0, ignore_index=True)

    def collapse_variants(self, gcm: pd.DataFrame) -> pd.DataFrame:
        """filter_nmem: keep models with >1 realization, average to per-(Model,SSP), then RECOMPUTE
        D_pr/D_tas from the collapsed totals vs the collapsed historical (not a mean of deltas)."""
        cfg = self.cfg
        gcm = gcm.copy()
        gcm["model_ssp"] = gcm["Model"] + "_" + gcm["SSP"]
        if not cfg.filter_nmem_gcms:
            return gcm

        base = gcm[gcm["period"] == cfg.base_center]
        counts = base.groupby(["Model", "SSP"]).size()
        nmem_keys = {f"{m}_{s}" for (m, s), n in counts.items() if n > 1}
        filt = gcm[gcm["model_ssp"].isin(nmem_keys)]

        num = ["Pr_total", "Tas_mean", "D_pr", "D_pr_lm", "D_tas"]
        hist_means = (
            filt[filt["period"] == cfg.base_center]
            .groupby(["Model", "SSP"], as_index=False)[num].mean()
        )
        out = [hist_means.assign(period=cfg.base_center)]
        for p in sorted(filt["period"].unique()):
            if p == cfg.base_center:
                continue
            fm = filt[filt["period"] == p].groupby(["Model", "SSP"], as_index=False)[num].mean()
            fm = fm.merge(
                hist_means[["Model", "SSP", "Pr_total", "Tas_mean"]],
                on=["Model", "SSP"], suffixes=("", "_hist"),
            )
            fm["D_pr"] = (fm["Pr_total"] - fm["Pr_total_hist"]) / fm["Pr_total_hist"] * 100.0
            fm["D_tas"] = fm["Tas_mean"] - fm["Tas_mean_hist"]
            fm["period"] = p
            out.append(fm[["Model", "SSP", *num, "period"]])
        return pd.concat(out, axis=0, ignore_index=True)

    def period_means(self, collapsed: pd.DataFrame) -> pd.DataFrame:
        """Per-period ensemble mean of (ΔT, ΔP) -> gcm_mean (columns period, DT, DP)."""
        gm = collapsed.groupby("period", as_index=False)[["D_tas", self.cfg.pr_type]].mean()
        gm.columns = ["period", "DT", "DP"]
        return gm

    def period_covariances(self, collapsed: pd.DataFrame) -> dict:
        """Per future period, the 2x2 covariance of (ΔT, ΔP) with ddof=1 (R cov default).

        Column/row order is (ΔT, ΔP) == (DTsig, D_pr_lm), matching the gcm_sigs contract.
        """
        cfg = self.cfg
        out = {}
        for p in sorted(collapsed["period"].unique()):
            if p == cfg.base_center:
                continue
            d = collapsed[collapsed["period"] == p]
            out[p] = np.cov(d["D_tas"].to_numpy(), d[cfg.pr_type].to_numpy(), ddof=1)
        return out


def decompose_precip_variance(prepared: pd.DataFrame, lm_fits: pd.DataFrame, collapsed: pd.DataFrame, cfg) -> pd.DataFrame:
    """var(ΔP) decomposed into within-model and between-model uncertainty, accumulated over the
    planning-horizon index k = period - base_center: var(ΔP) = A*k + B*k^2.

    D_pr_lm's own ensemble spread is *exactly* quadratic in k by construction (see
    ``ClimateDeltas.period_covariances``), with no linear component -- it only captures spread
    *between* models, not the year-to-year noise *within* each model's own projection. This combines
    both sources:

    - A (linear term): the year-to-year scatter of each model's 30-yr-rolling precip series around
      its OWN fitted trend line (the noise the regression doesn't explain), on the same % scale as
      D_pr_lm, averaged across models and divided by ``2*(window-1)`` -- an effective-sample-size
      correction for the autocorrelation of a rolling-window regression's residuals. This noise
      accumulates roughly linearly with lead time.
    - B (quadratic term): the variance, across the ensemble's distinct (model, ssp) points -- or
      (model, variant, ssp) points when ``cfg.filter_nmem_gcms`` is False, since ``collapsed`` then
      has one row per realization rather than one per (model, ssp) -- of the normalized lm slope,
      divided by ``2*n_unique_models`` -- the between-model spread term.
    """
    pr_col = "pr_roll" if cfg.source_kind == "loca2-flow" else "pr"
    lm_indexed = lm_fits.set_index("mvs")

    resid_pct_var = []
    for mvs, data in prepared.groupby("mvs"):
        y = data[pr_col].to_numpy()
        x = data["y"].to_numpy()
        reg = linear_regression(x, y)
        resid = y - (reg.slope * x + reg.intercept)
        lm_pr_hist = lm_indexed.loc[mvs, "slope"] * (cfg.base_center - 14) + lm_indexed.loc[mvs, "intercept"]
        resid_pct_var.append(np.var(resid / lm_pr_hist * 100.0, ddof=1))
    A = float(np.mean(resid_pct_var)) / (2 * (cfg.window - 1))

    group_cols = ["Model", "SSP"] if cfg.filter_nmem_gcms else ["Model", "Variant", "SSP"]
    beta = (
        collapsed[collapsed["period"] != cfg.base_center]
        .assign(k=lambda d: d["period"] - cfg.base_center)
        .assign(beta=lambda d: d["D_pr_lm"] / d["k"])
        .groupby(group_cols)["beta"].first()
    )
    n_unique_models = beta.index.get_level_values("Model").nunique()
    B = float(np.var(beta.to_numpy(), ddof=1)) / (2 * n_unique_models)

    periods = sorted(p for p in collapsed["period"].unique() if p != cfg.base_center)
    k = np.array([p - cfg.base_center for p in periods], dtype=float)
    return pd.DataFrame(
        {"period": periods, "var_dP": A * k + B * k**2, "A": A, "B": B}
    )


class BivariateNormalSurface:
    """Per-period normalized bivariate normal of (ΔT, ΔP) over the (Temp, Precip) stress grid.

    Replaces the R ``mvtnorm::dmvnorm(grid, mean, sigma)`` then ``probs / sum(probs)``. scipy's
    ``multivariate_normal.pdf`` equals ``dmvnorm`` exactly. Mean and covariance are paired by the
    SAME period (the legacy R loop paired mean[i] with sigma[i] offset by one -- treated as a bug,
    not reproduced, since we regenerate a consistent canonical set).
    """

    def __init__(self, cfg):
        self.cfg = cfg
        self._grid = cfg.grid.grid()  # columns T_lev, P_lev with Temp varying fastest (R parity)

    @property
    def grid(self) -> pd.DataFrame:
        return self._grid

    def evaluate(self, mean, cov) -> np.ndarray:
        """Normalized probability at each grid point for one period's (mean, cov)."""
        pts = self._grid[["T_lev", "P_lev"]].to_numpy()
        rv = multivariate_normal(
            mean=np.asarray(mean, float), cov=np.asarray(cov, float), allow_singular=True
        )
        dens = rv.pdf(pts)
        total = dens.sum()
        return dens / total if total > 0 else dens

    def period_surface(self, period, mean, cov) -> pd.DataFrame:
        out = self._grid.copy()
        out["Biv_Norm_Prob"] = self.evaluate(mean, cov)
        out["period"] = period
        return out

    def all_periods(self, means: pd.DataFrame, covs: dict) -> list:
        """One surface frame per future period, mean and cov aligned on the same period."""
        frames = []
        for period in sorted(covs):
            row = means.loc[means["period"] == period]
            mean = (float(row["DT"].iloc[0]), float(row["DP"].iloc[0]))
            frames.append(self.period_surface(period, mean, covs[period]))
        return frames

    @staticmethod
    def cumulative_area(period_df: pd.DataFrame) -> pd.DataFrame:
        """Sort by probability and take the cumulative sum (the contour 'area' transform)."""
        d = period_df.sort_values("Biv_Norm_Prob").copy()
        d["area"] = d["Biv_Norm_Prob"].cumsum()
        return d
