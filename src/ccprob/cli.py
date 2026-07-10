"""Command-line entry point: ``ccprob``.

Subcommands: ``list`` (domains), ``preprocess`` (Stage-1), ``distribution`` (Stage-2 CSVs),
``plot`` (contour figures), and ``run`` (the full preprocess -> distribution -> plots chain).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from . import correlation, historical, loc95, loc95_shift, precip_extremes
from .cmip5 import Cmip5Pipeline
from .config import CONFIGS_DIR, REPO_ROOT, load_domain
from .pipeline import LocaPipeline


def _pipeline(cfg):
    """Pick the pipeline for a domain (CMIP5 is self-contained; the rest run the LOCA-2 stages)."""
    return Cmip5Pipeline(cfg) if cfg.source_kind == "cmip5" else LocaPipeline(cfg)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="ccprob", description="Climate-change probability surfaces")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("list", help="list available domain configs")

    pp = sub.add_parser("preprocess", help="run Stage-1 preprocessing for a domain")
    pp.add_argument("domain")
    pp.add_argument("--out", default=None, help="output directory (default: processed/)")
    pp.add_argument("--no-write", action="store_true", help="compute only, do not write CSVs")

    dd = sub.add_parser("distribution", help="run Stage-2 (gcm_mean/gcm_sigs/biv_norm_vals)")
    dd.add_argument("domain")
    dd.add_argument("--out", default=None, help="output directory (default: processed/)")
    dd.add_argument("--no-write", action="store_true", help="compute only, do not write CSVs")

    pl = sub.add_parser("plot", help="render contour figure(s) for the configured plot periods")
    pl.add_argument("domain")
    pl.add_argument("--out", default=None, help="output directory (default: figures/)")
    pl.add_argument("--animate", action="store_true", help="also write the across-periods GIF")

    rn = sub.add_parser("run", help="full pipeline: preprocess -> distribution -> plots")
    rn.add_argument("domain")
    rn.add_argument("--out", default=None, help="output directory for CSVs (default: processed/)")
    rn.add_argument("--no-plots", action="store_true", help="skip figure rendering")
    rn.add_argument("--animate", action="store_true", help="also write the across-periods GIF")

    pe = sub.add_parser(
        "precip-extremes",
        help="sliding-window N-yr minimum average annual precip, by mvs, plus histogram figures",
    )
    pe.add_argument("domain")
    pe.add_argument("--out", default=None, help="output directory (default: processed/ + figures/)")
    pe.add_argument("--no-write", action="store_true", help="compute only, do not write CSV/figures")

    co = sub.add_parser(
        "correlation",
        help="time series of the Pearson correlation between DT and DP, across the lm/raw x "
        "varavg/novaravg distribution variants",
    )
    co.add_argument("domain")
    co.add_argument("--out", default=None, help="output directory (default: processed/ + figures/)")
    co.add_argument("--no-write", action="store_true", help="compute only, do not write CSV/figure")

    hi = sub.add_parser(
        "historical",
        help="aggregate the WGEN statewide daily grid to CalSim basin regions + flow-weighted CV average",
    )
    hi.add_argument(
        "--wgen-dir",
        default=r"C:\Users\warnold_la\Local\WGEN_NonDetrend_Unsplit_Statewide",
        help="directory of WGEN data_<lat>_<lon> grid-cell files (external to the repo)",
    )
    hi.add_argument(
        "--gpkg",
        default=r"C:\Users\warnold_la\Local\repos\SAC-SMA\data\calsim\gis\calsim3.gpkg",
        help="CalSim3 geopackage with basin watershed polygons",
    )
    hi.add_argument(
        "--xlsx",
        default=None,
        help="CalSim-19Basins-FlowContributionPercentage.xlsx (default: data/ in the repo)",
    )
    hi.add_argument("--out", default=None, help="output directory (default: data/historical/)")
    hi.add_argument(
        "--basins",
        default=None,
        help="comma-separated basin region ids to restrict to, e.g. 2,6 (default: all resolvable, "
        "plus the flow-weighted cv-flow-weighted aggregate)",
    )

    lo = sub.add_parser(
        "loc95",
        help="aggregate the LOC95 (1.8%%-dried) daily grid to CalSim basin regions + "
        "flow-weighted CV average",
    )
    lo.add_argument(
        "--loc-dir",
        default=r"C:\Users\warnold_la\Downloads\4\4",
        help="directory of LOC95 meteo_<lat>_<lon> grid-cell files (external to the repo)",
    )
    lo.add_argument(
        "--gpkg",
        default=r"C:\Users\warnold_la\Local\repos\SAC-SMA\data\calsim\gis\calsim3.gpkg",
        help="CalSim3 geopackage with basin watershed polygons",
    )
    lo.add_argument(
        "--xlsx",
        default=None,
        help="CalSim-19Basins-FlowContributionPercentage.xlsx (default: data/ in the repo)",
    )
    lo.add_argument("--out", default=None, help="output directory (default: data/loc/)")
    lo.add_argument(
        "--basins",
        default=None,
        help="comma-separated basin region ids to restrict to, e.g. 2,6 (default: all resolvable, "
        "plus the flow-weighted cv-flow-weighted aggregate)",
    )

    ls = sub.add_parser(
        "loc95-shift",
        help="overlay LOC95-vs-historical 30-yr moving-average (DT, DP) points on the "
        "raw+all-variants 2043 GCM contour figure (cv-flow-weighted only)",
    )
    ls.add_argument("--out", default=None, help="output figure path (default: figures/)")

    args = parser.parse_args(argv)

    if args.command == "list":
        for cfg_file in sorted(CONFIGS_DIR.glob("*.yaml")):
            print(cfg_file.stem)
        return 0

    if args.command == "preprocess":
        cfg = load_domain(args.domain)
        if cfg.source_kind == "cmip5":
            parser.error("cmip5 is self-contained (no preprocess stage); use distribution/plot/run")
        artifacts = LocaPipeline(cfg).run_preprocess(out_dir=args.out, write=not args.no_write)
        rows = {k: len(v) for k, v in artifacts.items()}
        print(f"preprocess complete for {args.domain}: {rows}")
        return 0

    if args.command == "distribution":
        cfg = load_domain(args.domain)
        result = _pipeline(cfg).run_distribution(out_dir=args.out, write=not args.no_write)
        print(
            f"distribution complete for {args.domain}: "
            f"{len(result['gcm_mean'])} periods, {len(result['covariances'])} covariances"
        )
        return 0

    if args.command == "plot":
        cfg = load_domain(args.domain)
        pipe = _pipeline(cfg)
        kwargs = {"animate": args.animate} if cfg.source_kind != "cmip5" else {}
        figures = pipe.run_plots(out_dir=args.out, write=True, **kwargs)
        print(f"plot complete for {args.domain}: {[str(p) for p in figures]}")
        return 0

    if args.command == "run":
        cfg = load_domain(args.domain)
        pipe = _pipeline(cfg)
        kwargs = {"animate": args.animate} if cfg.source_kind != "cmip5" else {}
        result = pipe.run(out_dir=args.out, write=True, plots=not args.no_plots, **kwargs)
        print(
            f"run complete for {args.domain}: {len(result['gcm_mean'])} periods, "
            f"{len(result['figures'])} figure(s)"
        )
        return 0

    if args.command == "precip-extremes":
        cfg = load_domain(args.domain)
        result = precip_extremes.run(cfg, out_dir=args.out, write=not args.no_write)
        print(f"precip-extremes complete for {args.domain}: {len(result['results'])} rows, csv={result['csv']}")
        return 0

    if args.command == "correlation":
        cfg = load_domain(args.domain)
        result = correlation.run(cfg, out_dir=args.out, write=not args.no_write)
        print(
            f"correlation complete for {args.domain}: {len(result['results'])} rows, "
            f"csv={result['csv']}, figure={result['figure']}"
        )
        return 0

    if args.command == "historical":
        xlsx_path = Path(args.xlsx) if args.xlsx else REPO_ROOT / "data" / "CalSim-19Basins-FlowContributionPercentage.xlsx"
        out_dir = Path(args.out) if args.out else REPO_ROOT / "data" / "historical"
        basins = [int(b) for b in args.basins.split(",")] if args.basins else None
        result = historical.run(args.wgen_dir, args.gpkg, xlsx_path, out_dir, basins=basins)
        print(
            f"historical complete: {len(result['basins'])} basin file(s), "
            f"cv_flow_weighted={result['cv_flow_weighted']}"
        )
        return 0

    if args.command == "loc95":
        xlsx_path = Path(args.xlsx) if args.xlsx else REPO_ROOT / "data" / "CalSim-19Basins-FlowContributionPercentage.xlsx"
        out_dir = Path(args.out) if args.out else REPO_ROOT / "data" / "loc"
        basins = [int(b) for b in args.basins.split(",")] if args.basins else None
        result = loc95.run(args.loc_dir, args.gpkg, xlsx_path, out_dir, basins=basins)
        print(
            f"loc95 complete: {len(result['basins'])} basin file(s), "
            f"cv_flow_weighted={result['cv_flow_weighted']}"
        )
        return 0

    if args.command == "loc95-shift":
        result = loc95_shift.plot_shift_overlay(out_path=args.out)
        print(f"loc95-shift complete: {len(result['shift'])} points, figure={result['figure']}")
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
