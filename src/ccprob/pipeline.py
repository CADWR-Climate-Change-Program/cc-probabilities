"""Domain pipeline orchestration.

Wires the Stage-1 preprocessing (Ensemble -> TrendModel -> processed/ CSVs), the Stage-2 stats
(ClimateDeltas -> gcm_mean / gcm_sigs / biv_norm_vals), and the contour figures. ``run`` chains
preprocess -> distribution -> plots.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

from ._rcompat import (
    write_biv_norm_vals_csv,
    write_gcm_mean_csv,
    write_gcm_points_csv,
    write_gcm_sigs_csv,
)
from .config import DomainConfig
from .distribution import BivariateNormalSurface, ClimateDeltas, decompose_precip_variance
from .plotting import ContourPlotter
from .projections import Ensemble
from .trends import TrendModel


def _variant_tokens(variant_cfg: DomainConfig) -> tuple[str, str]:
    avg = "varavg" if variant_cfg.filter_nmem_gcms else "novaravg"
    pr = "lm" if variant_cfg.pr_type == "D_pr_lm" else "raw"
    return avg, pr


def _variant_csv_name(kind: str, variant_cfg: DomainConfig, suffix: str = "") -> str:
    avg, pr = _variant_tokens(variant_cfg)
    tail = f"_{suffix}" if suffix else ""
    return f"{kind}_loca2_{avg}_{pr}_{variant_cfg.name}{tail}.csv"


def _variant_figure_name(period: int, variant_cfg: DomainConfig) -> str:
    avg, pr = _variant_tokens(variant_cfg)
    return f"loca2-biv-norm-prob-dt-dp-{pr}-{avg}-{period}_{variant_cfg.name}.svg"


def _variant_animation_name(variant_cfg: DomainConfig) -> str:
    avg, pr = _variant_tokens(variant_cfg)
    return f"loca2-biv-norm-prob-dt-dp-{pr}-{avg}-animation_{variant_cfg.name}.gif"


def _mean_tuple(gcm_mean, period) -> tuple[float, float]:
    """(DT, DP) for one period from a ``gcm_mean``-shaped frame, for the crosshair overlay."""
    row = gcm_mean.loc[gcm_mean["period"] == period]
    return float(row["DT"].iloc[0]), float(row["DP"].iloc[0])


def _title_qualifier(variant_cfg: DomainConfig) -> str:
    """Third title line: which realizations were averaged + how the precip change was computed."""
    avg = "variant-averaged" if variant_cfg.filter_nmem_gcms else "all variants"
    pr = "trend-horizon" if variant_cfg.pr_type == "D_pr_lm" else "point-horizon"
    return f"{avg}, {pr} precipitation change"


def _variant_title(period: int, variant_cfg: DomainConfig) -> str:
    base0, base1 = variant_cfg.base_center - 14, variant_cfg.base_center + 15
    return (
        "Projected range of likely climate changes\n"
        f"by {period} relative to the baseline 30-yr period {base0}-{base1}\n"
        f"{_title_qualifier(variant_cfg)}"
    )


class LocaPipeline:
    def __init__(self, cfg: DomainConfig):
        self.cfg = cfg

    def run_preprocess(self, out_dir=None, write: bool = True) -> dict:
        """Compute and (optionally) write the Stage-1 artifacts for this domain.

        Writes never target the committed goldens unless ``out_dir`` resolves to ``processed/`` and
        ``write=True`` -- callers/tests pass an out_dir to stay clear of the baselines.
        """
        cfg = self.cfg
        ens = Ensemble.from_config(cfg)
        tm = TrendModel(cfg)
        prepared = tm.prepare(ens)

        artifacts = {
            "lmfits": tm.lm_fits(prepared),
            "warming_levels": tm.warming_levels(prepared),
            "prepared": prepared,
        }
        if cfg.source_kind == "loca2-basin":
            artifacts["thirtyyr_avgs"] = prepared[
                ["mvs", "model", "variant", "ssp", "y", "pr", "tavg"]
            ]

        if write:
            out = Path(out_dir) if out_dir else cfg.paths["processed_dir"]
            out.mkdir(parents=True, exist_ok=True)
            artifacts["lmfits"].to_csv(out / cfg.lmfits_path.name, index=False)
            artifacts["warming_levels"].to_csv(out / cfg.warming_levels_path.name, index=False)
            if "thirtyyr_avgs" in artifacts:
                # mirror the legacy stray-index column the basin R script tolerates
                artifacts["thirtyyr_avgs"].to_csv(out / cfg.thirtyyr_avgs_path.name, index=True)

        return artifacts

    def _compute_distribution(self) -> dict:
        """Shared Stage-2 compute (no I/O): collapsed delta points, per-period mean/cov, surfaces.

        Flow sources deltas from the xlsx worksheets; basin sources them from the (in-memory)
        30yravgs table. The covariance is always the empirical np.cov of the collapsed points.
        """
        cfg = self.cfg
        if cfg.source_kind not in ("loca2-flow", "loca2-basin"):
            raise NotImplementedError(
                f"distribution supports loca2-flow/loca2-basin; {cfg.source_kind!r} not yet ported"
            )
        artifacts = self.run_preprocess(write=False)
        deltas = ClimateDeltas(cfg, artifacts["lmfits"], thirtyyr_avgs=artifacts.get("thirtyyr_avgs"))
        collapsed = deltas.collapse_variants(deltas.build())
        gcm_mean = deltas.period_means(collapsed)
        covs = deltas.period_covariances(collapsed)
        frames = BivariateNormalSurface(cfg).all_periods(gcm_mean, covs)
        return {
            "collapsed": collapsed,
            "gcm_mean": gcm_mean,
            "covariances": covs,
            "biv_norm_frames": frames,
            "prepared": artifacts["prepared"],
            "lmfits": artifacts["lmfits"],
        }

    def run_distribution(self, out_dir=None, write: bool = True) -> dict:
        """Stage-2: per-period ensemble mean + covariance + normalized bivariate-normal surfaces.

        Computes lm fits in-memory so gcm_mean / gcm_sigs / biv_norm_vals are internally consistent.
        """
        cfg = self.cfg
        dist = self._compute_distribution()
        gcm_mean, covs, frames = dist["gcm_mean"], dist["covariances"], dist["biv_norm_frames"]

        if write:
            out = Path(out_dir) if out_dir else cfg.paths["processed_dir"]
            out.mkdir(parents=True, exist_ok=True)
            if cfg.outputs.get("gcm_mean"):
                write_gcm_mean_csv(gcm_mean, out / cfg.gcm_mean_path.name)
            if cfg.outputs.get("gcm_sigs"):
                write_gcm_sigs_csv(covs, out / cfg.gcm_sigs_path.name)
            if cfg.outputs.get("biv_norm_vals"):
                write_biv_norm_vals_csv(frames, out / cfg.biv_norm_vals_path.name)
            if cfg.outputs.get("gcm_points"):
                write_gcm_points_csv(dist["collapsed"], cfg.pr_type, out / cfg.gcm_points_path.name)
            self._maybe_write_decomposed(cfg, dist, out)
            for overrides in cfg.outputs.get("variants", {}).values():
                variant_cfg = dataclasses.replace(cfg, **overrides)
                vdist = LocaPipeline(variant_cfg)._compute_distribution()
                if cfg.outputs.get("gcm_mean"):
                    write_gcm_mean_csv(vdist["gcm_mean"], out / _variant_csv_name("gcm_mean", variant_cfg))
                if cfg.outputs.get("gcm_sigs"):
                    write_gcm_sigs_csv(vdist["covariances"], out / _variant_csv_name("gcm_sigs", variant_cfg))
                if cfg.outputs.get("gcm_points"):
                    write_gcm_points_csv(
                        vdist["collapsed"], variant_cfg.pr_type,
                        out / _variant_csv_name("gcm_points", variant_cfg),
                    )
                self._maybe_write_decomposed(cfg, vdist, out, variant_cfg)

        return {"gcm_mean": gcm_mean, "covariances": covs, "biv_norm_frames": frames}

    def _maybe_write_decomposed(
        self, cfg: DomainConfig, dist: dict, out: Path, variant_cfg: DomainConfig | None = None
    ) -> None:
        """Write a decomposed gcm_sigs (see ``decompose_precip_variance``) for ``variant_cfg`` (or
        the base config when omitted) -- only when it's lm-based; the decomposition's between-model
        term relies on D_pr_lm's exact proportionality to k, which raw D_pr does not have."""
        target_cfg = variant_cfg or cfg
        if not (cfg.outputs.get("gcm_sigs_decomposed") and target_cfg.pr_type == "D_pr_lm"):
            return
        decomp = decompose_precip_variance(dist["prepared"], dist["lmfits"], dist["collapsed"], target_cfg)
        var_dp_by_period = dict(zip(decomp["period"], decomp["var_dP"]))
        decomposed_covs = {p: m.copy() for p, m in dist["covariances"].items()}
        for p, m in decomposed_covs.items():
            m[1, 1] = var_dp_by_period[p]
        name = (
            cfg.gcm_sigs_decomposed_path.name
            if variant_cfg is None
            else _variant_csv_name("gcm_sigs", variant_cfg, suffix="decomposed")
        )
        write_gcm_sigs_csv(decomposed_covs, out / name)

    def run_plots(self, out_dir=None, write: bool = True, animate: bool = False) -> list:
        """Render the per-period contour figure(s) for ``cfg.plot_periods`` and (optionally) the GIF.

        ``plot_periods`` are offsets from ``base_center`` (64 -> 2070), matching the legacy
        ``climate_period`` index. Returns the list of written file paths.

        Figures are rendered on ``cfg.plot_grid`` when set (a finer grid than ``cfg.grid``, which
        stays coarse for the ``biv_norm_vals`` CSV contract) -- the cloud is recomputed at plot
        resolution rather than reusing ``_compute_distribution()``'s coarse-grid surfaces.
        """
        cfg = self.cfg
        dist = self._compute_distribution()
        collapsed = dist["collapsed"]
        plot_cfg = dataclasses.replace(cfg, grid=cfg.plot_grid) if cfg.plot_grid else cfg
        frames = BivariateNormalSurface(plot_cfg).all_periods(dist["gcm_mean"], dist["covariances"])
        frames_by_period = {int(f["period"].iloc[0]): f for f in frames}

        plotter = ContourPlotter(plot_cfg)
        out = Path(out_dir) if out_dir else cfg.paths["figures_dir"]
        written = []
        if write:
            out.mkdir(parents=True, exist_ok=True)

        for offset in cfg.plot_periods:
            period = cfg.base_center + offset
            frame = frames_by_period[period]
            points = collapsed[collapsed["period"] == period]
            if write:
                name = f"loca2-biv-norm-prob-dt-dp-lm-varavg-{period}_{cfg.name}.svg"
                plotter.filled_contour(
                    frame, points, title=_variant_title(period, cfg), out_path=out / name,
                    mean=_mean_tuple(dist["gcm_mean"], period), cov=dist["covariances"][period],
                )
                written.append(out / name)

        if animate and write and cfg.animation_periods:
            start, stop = cfg.animation_periods
            anim_frames = [
                frames_by_period[cfg.base_center + off]
                for off in range(start, stop + 1)
                if cfg.base_center + off in frames_by_period
            ]
            gif = out / f"loca2-biv-norm-prob-dt-dp-lm-varavg-animation_{cfg.name}.gif"
            plotter.animate(
                anim_frames, collapsed, out_path=gif,
                means=dist["gcm_mean"], covs=dist["covariances"],
                title_fn=lambda period, c=cfg: _variant_title(period, c),
            )
            written.append(gif)

        if write:
            animated_variants = set(cfg.outputs.get("animated_variants", []))
            for name, overrides in cfg.outputs.get("variants", {}).items():
                written += self._run_variant_plots(
                    plot_cfg, overrides, out, animate=(animate and name in animated_variants)
                )

        return written

    def _run_variant_plots(
        self, plot_cfg: DomainConfig, overrides: dict, out: Path, animate: bool = False
    ) -> list:
        """Figures from a distribution variant (e.g. no variant-averaging collapse, or raw
        unfitted D_pr instead of the lm-detrended precip change) -- same grid/styling as the
        default figures, distinguished by filename and title. Also renders that variant's
        ``animation_periods`` GIF when ``animate`` (config-gated via ``outputs.animated_variants``)."""
        variant_cfg = dataclasses.replace(plot_cfg, **overrides)
        variant_dist = LocaPipeline(variant_cfg)._compute_distribution()
        frames = BivariateNormalSurface(plot_cfg).all_periods(
            variant_dist["gcm_mean"], variant_dist["covariances"]
        )
        frames_by_period = {int(f["period"].iloc[0]): f for f in frames}
        plotter = ContourPlotter(variant_cfg)

        written = []
        for offset in self.cfg.plot_periods:
            period = self.cfg.base_center + offset
            frame = frames_by_period[period]
            points = variant_dist["collapsed"][variant_dist["collapsed"]["period"] == period]
            name = _variant_figure_name(period, variant_cfg)
            plotter.filled_contour(
                frame, points, title=_variant_title(period, variant_cfg), out_path=out / name,
                mean=_mean_tuple(variant_dist["gcm_mean"], period),
                cov=variant_dist["covariances"][period],
            )
            written.append(out / name)

        if animate and self.cfg.animation_periods:
            start, stop = self.cfg.animation_periods
            anim_frames = [
                frames_by_period[self.cfg.base_center + off]
                for off in range(start, stop + 1)
                if self.cfg.base_center + off in frames_by_period
            ]
            gif = out / _variant_animation_name(variant_cfg)
            plotter.animate(
                anim_frames, variant_dist["collapsed"], out_path=gif,
                title_fn=lambda period, vc=variant_cfg: _variant_title(period, vc),
                means=variant_dist["gcm_mean"], covs=variant_dist["covariances"],
            )
            written.append(gif)
        return written

    def run(self, out_dir=None, write: bool = True, plots: bool = True, animate: bool = False) -> dict:
        """End-to-end: preprocess -> distribution -> plots."""
        self.run_preprocess(out_dir=out_dir, write=write)
        dist = self.run_distribution(out_dir=out_dir, write=write)
        figures = self.run_plots(out_dir=out_dir, write=write, animate=animate) if plots else []
        return {**dist, "figures": figures}
