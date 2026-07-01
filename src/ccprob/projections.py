"""Ensemble loading: combine the per-file projection series into one tidy long frame.

Centralizes the per-file load loop (duplicated three times in the legacy scripts) and the
historical-replication rule (flow only) into a single place.
"""

from __future__ import annotations

import glob

import pandas as pd

from .io import FilenameParser, ProjectionReader


class Ensemble:
    """A combined long frame of every model_variant_ssp (mvs) series for a domain.

    Columns: ``mvs, model, variant, ssp, y, pr, tavg`` (one row per mvs-year).
    """

    def __init__(self, frame: pd.DataFrame, cfg):
        self.frame = frame
        self.cfg = cfg

    @classmethod
    def from_config(cls, cfg, reader: ProjectionReader | None = None) -> "Ensemble":
        reader = reader or ProjectionReader(cfg)
        parser = FilenameParser(cfg.filename_field_order)
        files = sorted(glob.glob(str(cfg.paths["input_glob"])))
        if not files:
            raise FileNotFoundError(f"No projection files match {cfg.paths['input_glob']}")
        meta = [parser.parse(f) for f in files]

        pieces: list[pd.DataFrame] = []
        if cfg.source_kind == "loca2-flow":
            meta_df = pd.DataFrame(meta)
            nonhist = meta_df[meta_df["ssp"] != "historical"]
            for path, m in zip(files, meta):
                model, variant, ssp = m["model"], m["variant"], m["ssp"]
                series = reader.read_annual(path)
                series.insert(0, "model", model)
                series.insert(1, "variant", variant)
                if ssp == "historical":
                    # replicate the historical series once per ssp this model+variant actually ran
                    ssps = nonhist[(nonhist.model == model) & (nonhist.variant == variant)][
                        "ssp"
                    ].unique()
                    for s in ssps:
                        c = series.copy()
                        c.insert(2, "ssp", s)
                        c.insert(0, "mvs", FilenameParser.mvs(model, variant, s))
                        pieces.append(c)
                else:
                    series.insert(2, "ssp", ssp)
                    series.insert(0, "mvs", FilenameParser.mvs(model, variant, ssp))
                    pieces.append(series)
        elif cfg.source_kind == "loca2-basin":
            for path, m in zip(files, meta):
                model, variant, ssp = m["model"], m["variant"], m["ssp"]
                series = reader.read_30y(path)
                series.insert(0, "model", model)
                series.insert(1, "variant", variant)
                series.insert(2, "ssp", ssp)
                series.insert(0, "mvs", FilenameParser.mvs(model, variant, ssp))
                pieces.append(series)
        else:
            raise ValueError(f"Ensemble does not support source_kind={cfg.source_kind!r}")

        frame = pd.concat(pieces, axis=0, ignore_index=True)
        # guarantee per-mvs year-ascending order so the rolling window is well-defined
        frame = frame.sort_values(["mvs", "y"]).reset_index(drop=True)
        return cls(frame, cfg)

    def long_frame(self) -> pd.DataFrame:
        return self.frame

    @property
    def mvs_list(self):
        return self.frame["mvs"].unique()
