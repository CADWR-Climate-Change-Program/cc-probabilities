# cc-probabilities

Generates probability surfaces — joint distributions of temperature change (ΔT, °C) and
precipitation change (ΔP, %) — describing the range of plausible future climate for the California
Central Valley, for use in water-resources vulnerability and decision-scaling studies. For each
future period it fits a bivariate normal to a GCM ensemble and renders a filled-contour "stress
test" figure (ΔP on the x-axis, ΔT on the y-axis) with the 68% and 95% probability contours and
the individual models overlaid. Two ensembles are supported: legacy **CMIP5** and current
**CMIP6 / LOCA-2** downscaled projections. The fitted parameters (`gcm_mean_*` + `gcm_sigs_*`) are
also the input to the companion **`risk-informed-scenarios`** repo, which converts them into
water-system "Levels of Concern."

The pipeline is the installable Python package **`ccprob`** (`src/ccprob/`), driven by a CLI and
per-domain YAML config (`configs/*.yaml`) — no hardcoded paths or globals to hand-edit between runs.
The six original Python/R research scripts this package replaced are kept under `legacy/` as a
frozen historical reference (and the oracle the package's golden tests were validated against); see
`legacy/README.md`.

## Pipeline

Every run targets a spatial **domain** — `cv-flow-weighted` (Central-Valley-wide, flow-weighted),
a basin such as `basin-06`, or `cmip5` — selected by name and resolved entirely from
`configs/<domain>.yaml`. ΔT and ΔP are always changes relative to a 30-year baseline window centered
on 2006 (1992–2021). The LOCA-2 domains run two stages in-process (preprocess, then distribution);
CMIP5 is self-contained (no preprocessing stage).

```
LOCA-2 (CMIP6) — ccprob preprocess <domain>, then ccprob distribution/plot <domain>
────────────────────────────────────────────────────────────────────────────────────
  data/loca2/loca2-projections-*/        ccprob.trends / .projections     processed/
  raw projection CSVs                     · 30-yr rolling P/T trend       loca2_lmfits_1981_<domain>.csv
  daily · annual · 30-yr,       ───────▶   · +1…+5 °C crossing years ──▶  loca2_warming_levels_<domain>.csv
  one file per model_variant_ssp                                          loca2_30yravgs_basin-<NN>.csv
                                                                                      │
  data/loca2/...summary-worksheets/*.xlsx                                            │
  per-period 30-yr P & T climatologies ──────────────────┐                           │
                                                          ▼                           ▼
                                              ccprob.distribution (reads both inputs)
                                               · ensemble mean + covariance of (ΔT, ΔP) per period
                                               · evaluate a bivariate normal over a ΔT × ΔP grid
                                               · normalize to a probability surface
                                                          │
                       ┌──────────────────────────────────┼─────────────────────────────────┐
                       ▼                                   ▼                                  ▼
               processed/                            figures/                       handoff to downstream
               biv_norm_vals_dt-dp_*.csv             *.svg contour                  gcm_mean_loca2_varavg_lm_<domain>.csv
               gcm_mean_loca2_varavg_lm_<domain>.csv "stress test" plots            gcm_sigs_loca2_varavg_lm_<domain>.csv
               gcm_sigs_loca2_varavg_lm_<domain>.csv (+ animation .gif, flow only)          │
                                                                                            ▼
                                                                          ../risk-informed-scenarios

CMIP5 — self-contained, no preprocessing stage
───────────────────────────────────────────────
  data/cmip5/*.xlsx + data/gcm-families.csv ──▶ ccprob.cmip5 ──▶ processed/gcm_mean_cmip5_family.csv
                                                                  figures/biv-norm-prob-dt-dp-cmip5-*.svg
```

## Running

Requires Python ≥3.10 with `numpy pandas scipy openpyxl matplotlib pyyaml` (`pip install -e .` from
the repo root installs these and the `ccprob` console script; `pip install -e ".[dev]"` adds
`pytest`).

```
ccprob list                                  # available domains (configs/*.yaml)
ccprob preprocess cv-flow-weighted           # Stage-1: lmfits / warming-levels / 30yr-avgs
ccprob distribution cv-flow-weighted         # Stage-2: gcm_mean / gcm_sigs / biv_norm_vals
ccprob plot cv-flow-weighted [--animate]     # contour figure(s) (+ across-periods GIF)
ccprob run cv-flow-weighted                  # the full chain; --out DIR to avoid processed/figures/
```

## Notes

- **`model_variant_ssp` ("mvs")**, e.g. `ACCESS-CM2_r1i1p1f1_ssp245`, is the unit of one projection
  time series and the thing most loops iterate over.
- The two LOCA-2 input directories order filename fields differently — `annual-flow` is
  `model_ssp_variant`, `30y-basin` is `model_variant_ssp`. This is config-driven
  (`filename_field_order` per domain in `configs/*.yaml`), not duplicated parsing code.
- `D_pr_lm` (linear-model-detrended precip change) vs `D_pr` (raw 30-yr change) is the `pr_type`
  setting per domain; `filter_nmem_gcms` averages a model's multiple realizations into one point
  before fitting.
- **Regenerating the downstream handoff:** `ccprob distribution cv-flow-weighted` (or `ccprob run`)
  writes `gcm_mean_loca2_varavg_lm_cv-flow-weighted.csv` and
  `gcm_sigs_loca2_varavg_lm_cv-flow-weighted.csv` into `processed/` — copy both into
  `risk-informed-scenarios/data/`. Preserve their on-disk layout exactly (the downstream code reads
  `gcm_sigs` as a transposed, row-label-indexed table); see `_rcompat.write_gcm_sigs_csv` for the
  exact layout contract.
