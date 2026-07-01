# Legacy scripts (frozen historical reference)

These are the original six research scripts this repo ran on before the `ccprob` Python package
replaced them (see the top-level `README.md` / `CLAUDE.md`). They are kept, unmodified except for
the move itself, as the historical oracle the package was validated against — every `processed/`
golden CSV these scripts produced is what `ccprob`'s test suite regresses against.

They are **not maintained going forward** and are not part of the installable package. Use
`ccprob` for anything new. This directory is slated for removal after one release once `ccprob` has
been the primary path long enough to build confidence.

| Script | Replaced by |
|---|---|
| `loca2-pr-trends-flow-weighted.py`, `loca2-pr-trends-basin.py` | `ccprob preprocess <domain>` (`ccprob.trends`, `ccprob.projections`) |
| `createBivNormValues - LOCA-2.R`, `createBivNormValues - LOCA-2-basin.R` | `ccprob distribution <domain>` / `ccprob plot <domain>` (`ccprob.distribution`, `ccprob.plotting`) |
| `createBivNormValues.R` (CMIP5) | `ccprob distribution cmip5` / `ccprob plot cmip5` (`ccprob.cmip5`) |
| `loca2-pr-percentiles.py` | not ported (exploratory daily-precip analysis, not part of the bivariate-normal pipeline; had no golden output to validate against) |
