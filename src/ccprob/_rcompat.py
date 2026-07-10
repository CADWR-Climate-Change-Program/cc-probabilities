"""R-parity helpers for the wide CSV layouts the legacy ``write.csv`` produced.

These layouts are a contract: ``gcm_sigs`` in particular is read by the downstream
risk-informed-scenarios repo via ``pd.read_csv(...).T`` + a two-rows-per-period slice. The grid
ordering (Temp varies fastest) lives in ``config.GridSpec.grid``; here we provide R-style
``make.unique`` headers and the wide writers.
"""

from __future__ import annotations

import pandas as pd


def make_unique(names) -> list[str]:
    """Replicate R ``make.unique``: first occurrence bare, repeats get .1, .2, ... suffixes."""
    seen: dict[str, int] = {}
    out = []
    for n in names:
        if n in seen:
            seen[n] += 1
            out.append(f"{n}.{seen[n]}")
        else:
            seen[n] = 0
            out.append(n)
    return out


def write_gcm_mean_csv(gcm_mean: pd.DataFrame, path) -> None:
    """Header ``period,DT,DP``; no index (matches the committed gcm_mean)."""
    gcm_mean.to_csv(path, index=False)


def write_gcm_sigs_csv(covs: dict, path, row_labels=("DTsig", "D_pr_lm"), col_labels=None) -> None:
    """Wide gcm_sigs: 2 rows x 2N cols, R make.unique headers, leading label column.

    Per period p the 2x2 [[var(DT), cov], [cov, var(DP)]] occupies two columns:
      col 2k (1st label)   = [var(DT), cov]   (row 0, row 1)
      col 2k+1 (2nd label) = [cov, var(DP)]
    The downstream consumer of the LOCA-2 flow file depends on exactly this layout -- do not
    "tidy" it. ``row_labels``/``col_labels`` differ by source: LOCA-2 uses (DTsig, D_pr_lm);
    CMIP5 uses (DT, DP).
    """
    col_labels = tuple(col_labels) if col_labels is not None else tuple(row_labels)
    dt_row, dp_row, names = [], [], []
    for p in sorted(covs):
        m = covs[p]
        dt_row.extend([m[0, 0], m[0, 1]])  # row 0 across the period's two columns
        dp_row.extend([m[1, 0], m[1, 1]])  # row 1 across the period's two columns
        names.extend(col_labels)
    df = pd.DataFrame([dt_row, dp_row], index=list(row_labels), columns=make_unique(names))
    df.to_csv(path, index=True)


def write_gcm_points_csv(collapsed: pd.DataFrame, pr_col: str, path) -> None:
    """Per-GCM (period, Model, SSP, DT, DP) points -- the un-aggregated inputs to gcm_mean/gcm_sigs.

    ``pr_col`` selects which precip-change column (``D_pr`` or ``D_pr_lm``) becomes ``DP``, so the
    same writer covers both the point-horizon and trend-horizon precipitation variants.
    """
    out = collapsed[["period", "Model", "SSP", "D_tas", pr_col]].rename(
        columns={"D_tas": "DT", pr_col: "DP"}
    )
    out.to_csv(path, index=False)


def write_biv_norm_vals_csv(period_frames, path) -> None:
    """Wide biv_norm_vals: per period a (T_lev, P_lev, Biv_Norm_Prob, period) block, concatenated
    side by side with make.unique headers and a 1-based row index (matches the committed layout)."""
    blocks = [
        f[["T_lev", "P_lev", "Biv_Norm_Prob", "period"]].reset_index(drop=True)
        for f in period_frames
    ]
    wide = pd.concat(blocks, axis=1)
    wide.columns = make_unique(list(wide.columns))
    wide.index = range(1, len(wide) + 1)
    wide.to_csv(path, index=True)
