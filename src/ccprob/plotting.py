"""Decision-scaling contour figures (the former Stage-2 R plotting).

``ContourPlotter`` replaces the copy-pasted ``my.filled.contour`` R function. A per-period
bivariate-normal surface (``T_lev``, ``P_lev``, ``Biv_Norm_Prob``) is turned into the cumulative
"area" matrix R built (sort by probability, ``cumsum``, ``spread`` to a P x T grid), drawn as a
99-band grayscale filled contour, then overlaid with the per-SSP GCM points and the ``1 - area``
highest-density contour lines at 0.68 / 0.95.

The module forces the headless ``Agg`` backend -- this is batch figure generation, not interactive.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: figures are written to files, never shown

import matplotlib.patheffects as patheffects  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from scipy.stats import norm  # noqa: E402

# R: colorRampPalette(c("white", "#d0d0d0", "#373737")); the dark gray is also the contour-line color.
GRAY_RAMP = ["white", "#d0d0d0", "#373737"]
CONTOUR_COLOR = "#373737"

# two-sided 95% range of a normal distribution: mean +/- Z_95 * sd
_Z_95 = float(norm.ppf(0.975))


class ContourPlotter:
    """Render the (ΔP, ΔT) probability surface for one or more periods.

    Styling is config-driven so the LOCA-2 (gray ramp, SSP points) and CMIP5 (blue ramp, RCP points)
    figures share one code path. Defaults reproduce the LOCA-2 flow figure; CMIP5 overrides come from
    ``cfg.extras`` (``fill_ramp``, ``contour_color``, ``scenario_col``, ``point_x_col``,
    ``point_y_col``, ``point_marker``).
    """

    def __init__(self, cfg):
        self.cfg = cfg
        e = cfg.extras
        self.fill_ramp = e.get("fill_ramp", GRAY_RAMP)
        self.contour_color = e.get("contour_color", CONTOUR_COLOR)
        self.scenario_col = e.get("scenario_col", "SSP")
        self.point_x_col = e.get("point_x_col", cfg.pr_type)
        self.point_y_col = e.get("point_y_col", "D_tas")
        self.point_marker = e.get("point_marker", "o")

    # ------------------------------------------------------------------ transforms (no rendering)
    @staticmethod
    def area_pivot(period_df):
        """Surface frame -> ``(temp_axis, precip_axis, Z)`` with ``Z[temp, precip]`` = cumulative area.

        Mirrors the R pipeline exactly: order ascending by ``Biv_Norm_Prob`` (stable, to match R's
        ``order``), ``cumsum`` into ``area``, then ``spread(key=T_lev, value=area)``. R indexes the
        matrix rows by P_lev and cols by T_lev and plots it with x=Precip, y=Temp; matplotlib wants
        ``Z[row=y, col=x]`` so we return the transpose (``Z`` indexed temp-by-precip).
        """
        d = period_df.sort_values("Biv_Norm_Prob", kind="stable")
        d = d.assign(area=d["Biv_Norm_Prob"].cumsum())
        piv = (
            d.pivot(index="T_lev", columns="P_lev", values="area")
            .sort_index(axis=0)
            .sort_index(axis=1)
        )
        return piv.index.to_numpy(float), piv.columns.to_numpy(float), piv.to_numpy(float)

    def levels(self) -> np.ndarray:
        """Filled-contour breakpoints: R ``seq(0, 1, length.out = n) ** 5`` (clustered near 0)."""
        n = self.cfg.prob_interval_count
        return np.linspace(0.0, 1.0, n) ** 5

    def fill_colors(self) -> list:
        """``n - 1`` band colors from the configured ramp (R drops the last of ``n``)."""
        n = self.cfg.prob_interval_count
        cmap = LinearSegmentedColormap.from_list("ccprob_fill", self.fill_ramp, N=n)
        return [cmap(i) for i in range(n)][:-1]

    def title_for(self, period: int) -> str:
        """R's two-line title: projected year + the 30-yr baseline window around ``base_center``."""
        base0 = self.cfg.base_center - 14
        base1 = self.cfg.base_center + 15
        return (
            f"Projected range of likely climate changes by {period}\n"
            f"relative to the baseline 30-yr period {base0}-{base1}"
        )

    def _xlim(self, override=None):
        if override is not None:
            return tuple(override)
        x = self.cfg.extras.get("plot_xlim")
        return tuple(x) if x else (self.cfg.grid.precip_min, self.cfg.grid.precip_max)

    def _ylim(self, override=None):
        if override is not None:
            return tuple(override)
        y = self.cfg.extras.get("plot_ylim")
        return tuple(y) if y else (self.cfg.grid.temp_min, self.cfg.grid.temp_max)

    def _draw_crosshairs(self, ax, mean, cov):
        """Ensemble mean +/- 95% normal-fit range on each axis, as a crosshair with value labels.

        Drawn twice (a wide white halo under a narrow black line) rather than via ``path_effects``
        on ``errorbar`` -- the error-bar caps are ``LineCollection``s, not ``Line2D``s, and don't
        reliably take a ``path_effects`` kwarg the way a single ``Text``/``Line2D`` does; two passes
        works regardless of artist type and keeps the crosshair legible over both the dark contour
        fill and the light margins.
        """
        mean_t, mean_p = float(mean[0]), float(mean[1])
        half_t = _Z_95 * float(np.sqrt(cov[0][0]))
        half_p = _Z_95 * float(np.sqrt(cov[1][1]))
        xerr, yerr = [[half_p], [half_p]], [[half_t], [half_t]]

        for color, lw, ms, z in (("white", 4.5, 16, 6), ("black", 1.8, 11, 7)):
            ax.errorbar(
                [mean_p], [mean_t], xerr=xerr, yerr=yerr, fmt="+", color=color,
                markersize=ms, markeredgewidth=lw, elinewidth=lw, capsize=7, capthick=lw,
                zorder=z,
            )

        label = f"ΔT = {mean_t:.1f} ± {half_t:.1f}°C\nΔP = {mean_p:.1f} ± {half_p:.0f}%"
        txt = ax.annotate(
            label, xy=(mean_p, mean_t), xytext=(10, 10), textcoords="offset points",
            fontsize=9, fontweight="bold", color="black", zorder=8, ha="left", va="bottom",
        )
        txt.set_path_effects([patheffects.withStroke(linewidth=3, foreground="white")])

    # ------------------------------------------------------------------ drawing
    def _draw(
        self, ax, period_df, point_df=None, *, title=None, xlim=None, ylim=None, pr_col=None,
        mean=None, cov=None,
    ):
        """Draw one period onto ``ax`` (shared by ``filled_contour`` and ``animate``)."""
        cfg = self.cfg
        x_col = pr_col or self.point_x_col
        y_col = self.point_y_col
        temp, precip, Z = self.area_pivot(period_df)

        ax.contourf(precip, temp, Z, levels=self.levels(), colors=self.fill_colors())
        cs = ax.contour(
            precip, temp, 1.0 - Z, levels=list(cfg.prob_levels),
            colors=self.contour_color, linewidths=2.5,
        )
        if cs.allsegs and any(len(seg) for seg in cs.allsegs):
            ax.clabel(cs, fmt="%.2f", fontsize=9)

        if point_df is not None and len(point_df):
            for scenario, color in cfg.ssp_colors.items():
                sub = point_df[point_df[self.scenario_col] == scenario]
                if len(sub):
                    ax.scatter(
                        sub[x_col], sub[y_col], s=28, c=color, marker=self.point_marker,
                        edgecolors="none", zorder=5, label=scenario,
                    )
            ax.legend(loc="best", fontsize=8, framealpha=0.9)

        if mean is not None and cov is not None:
            self._draw_crosshairs(ax, mean, cov)

        ax.set_xlabel("Change in Precipitation (%)")
        ax.set_ylabel("Change in Temperature (C)")
        if title:
            # wrap=True matters most for animate(): a fixed-canvas GIF frame has no bbox_inches="tight"
            # to auto-expand around a long title (e.g. a variant's qualifier suffix), so an unwrapped
            # title would get clipped at the figure edge instead of just wrapping to a new line.
            ax.set_title(title, fontsize=10, wrap=True)
        ax.set_xlim(*self._xlim(xlim))
        ax.set_ylim(*self._ylim(ylim))
        return ax

    def filled_contour(
        self, period_df, point_df=None, *, period=None, title=None, xlim=None, ylim=None,
        pr_col=None, out_path=None, figsize=(5.5, 7.0), mean=None, cov=None,
    ):
        """Render one period. Writes ``out_path`` (svg/png by extension) or returns ``(fig, ax)``.

        ``mean``/``cov`` (the period's ensemble (ΔT, ΔP) mean and 2x2 covariance, same pairing as
        ``BivariateNormalSurface.all_periods``) are optional -- when given, a crosshair is drawn at
        the mean with normal-fit 95% range whiskers on each axis.
        """
        if title is None and period is not None:
            title = self.title_for(int(period))
        fig, ax = plt.subplots(figsize=figsize)
        self._draw(ax, period_df, point_df, title=title, xlim=xlim, ylim=ylim, pr_col=pr_col, mean=mean, cov=cov)
        if out_path is None:
            return fig, ax
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        return out_path

    def animate(
        self, frames, point_df=None, *, out_path, fps=4, figsize=(5.5, 7.0),
        xlim=None, ylim=None, pr_col=None, title_fn=None, dpi=150, means=None, covs=None,
    ):
        """One frame per surface in ``frames`` (period-ordered) -> animated GIF via PillowWriter.

        ``figsize`` matches ``filled_contour``'s default so the animation has the same aspect/axis
        appearance as the static per-period figures; ``xlim``/``ylim`` fall back to the same
        ``cfg.extras`` bounds via ``_draw`` when left ``None``, for the same reason. ``dpi`` defaults
        higher than matplotlib's usual 100 since a GIF is rasterized (unlike the vector SVG figures)
        and benefits from the extra resolution.

        ``title_fn(period) -> str`` overrides the per-frame title (defaults to ``self.title_for``)
        -- used to carry a variant's qualifier (e.g. "all variants, raw precipitation change") into
        its animation, matching that variant's static figures.

        ``means``/``covs`` (a ``gcm_mean``-shaped frame with ``period``/``DT``/``DP`` columns, and a
        ``{period: 2x2 cov}`` dict -- the same pair ``run_plots`` threads into ``filled_contour``)
        are optional; when both are given, each frame gets the same mean/95%-range crosshair as the
        matching static figure.
        """
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        ordered = sorted(frames, key=lambda f: f["period"].iloc[0])
        title_fn = title_fn or self.title_for

        fig, ax = plt.subplots(figsize=figsize)
        writer = PillowWriter(fps=fps)
        try:
            with writer.saving(fig, str(out_path), dpi=dpi):
                for fr in ordered:
                    ax.clear()
                    period = int(fr["period"].iloc[0])
                    pts = None
                    if point_df is not None:
                        pts = point_df[point_df["period"] == period]
                    mean = cov = None
                    if means is not None and covs is not None and period in covs:
                        row = means.loc[means["period"] == period]
                        if len(row):
                            mean = (float(row["DT"].iloc[0]), float(row["DP"].iloc[0]))
                            cov = covs[period]
                    self._draw(
                        ax, fr, pts, title=title_fn(period),
                        xlim=xlim, ylim=ylim, pr_col=pr_col, mean=mean, cov=cov,
                    )
                    writer.grab_frame()
        finally:
            plt.close(fig)
        return out_path
