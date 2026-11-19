"""Retrieve and strictly audit public EltonTraits 1.0 bird nectar diet percentages.

This is an independent bird dietary descriptor only, *not* a range map,
island guild observation, or evidence of effective pollination. No eBird.
Original ecological data: Wilman et al. 2014 Ecology, DOI 10.1890/13-1917.1.
Official figshare collection 3306933, BirdFuncDat file 5631081 (CC BY 4.0).
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

SOURCE_DOI = "10.6084/m9.figshare.c.3306933.v1"
PAPER_DOI = "10.1890/13-1917.1"
OFFICIAL = (
    "https://ndownloader.figshare.com/files/5631081",
    "https://figshare.com/ndownloader/files/5631081",
)
MAX_BYTES = 10_000_000
MIN_SPECIES = 8000
MAX_SPECIES = 12000


def get_original_bytes(urls: tuple[str, ...] = OFFICIAL) -> tuple[bytes, str]:
    errors = []
    for url in urls:
        if url not in OFFICIAL:
            raise ValueError("Unapproved or non-original EltonTraits URI")
        try:
            req = Request(url, headers={"User-Agent": "fd-public-eltontraits-bird-v1",
                                        "Accept": "text/plain,text/tab-separated-values,*/*"})
            with urlopen(req, timeout=40) as f:
                blob = f.read(MAX_BYTES + 1)
            if not blob or len(blob) > MAX_BYTES:
                raise ValueError("Empty or excessively large original EltonTraits file")
            if b"<html" in blob[:300].lower():
                raise ValueError("Webpage substituted for source data")
            return blob, url
        except Exception as exc:
            errors.append(type(exc).__name__ + ": " + str(exc)[:150])
    raise RuntimeError("Official EltonTraits figshare data unavailable: " + "; ".join(errors))


def audit_bird_diet(blob: bytes, *, enforce_species_range: bool = True
                    ) -> tuple[dict, pd.DataFrame]:
    """Require raw original BirdFuncDat schema; never infer missing nectar as zero."""
    digest = hashlib.sha256(blob).hexdigest()
    errors = []
    frame = None
    for encoding in ("utf-8-sig", "latin-1"):
        try:
            frame = pd.read_csv(io.BytesIO(blob), sep="\t", dtype=str,
                                encoding=encoding, keep_default_na=False)
            break
        except (UnicodeError, pd.errors.ParserError) as exc:
            errors.append(type(exc).__name__)
    if frame is None:
        raise ValueError(f"Unable to parse original tab-delimited text: {errors}")
    if frame.columns.duplicated().any():
        raise ValueError("Duplicate original EltonTraits column names")
    required = {"Scientific", "Diet.Nect"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Original BirdFuncDat header missing {sorted(required-set(frame.columns))}: {list(frame.columns)[:12]}")
    if enforce_species_range and not MIN_SPECIES <= len(frame) <= MAX_SPECIES:
        raise ValueError(f"Unexpected EltonTraits source bird population: {len(frame)}")
    names = frame["Scientific"].str.strip().str.replace("_", " ", regex=False)
    if names.eq("").any() or names.duplicated().any():
        raise ValueError("Missing or repeated scientific species name; source manual audit required")
    raw = frame["Diet.Nect"].astype(str).str.strip()
    numeric = pd.to_numeric(raw.replace("", np.nan), errors="coerce")
    if not np.isfinite(numeric).all() or numeric.lt(0).any() or numeric.gt(100).any():
        raise ValueError("Missing/invalid/non-percentage Diet.Nect; do not impute zeros")
    weights = pd.DataFrame({
        "original_scientific_name": names,
        "diet_nectar_fraction": numeric.to_numpy(float) / 100.0,
    })
    receipt = {
        "contract": "fd_original_eltontraits_bird_nectar_v1",
        "source_collection_doi": SOURCE_DOI,
        "source_paper_doi": PAPER_DOI,
        "source_filename": "BirdFuncDat.txt",
        "source_file_id": 5631081,
        "original_sha256": digest,
        "n_original_bird_records": int(len(frame)),
        "n_recorded_nectar_positive": int(numeric.gt(0).sum()),
        "n_recorded_nectar_zero": int(numeric.eq(0).sum()),
        "sum_of_nectar_fraction_by_species": float((numeric / 100).sum()),
        "mean_nectar_fraction_by_species": float((numeric / 100).mean()),
        "source_license": "CC BY 4.0; attribute Wilman et al. (2014)",
        "has_eBird_taxonomy_or_occurrences": False,
        "observational_guild_presence": False,
        "licensed_spatial_range_overlays_available": False,
        "island_source_animal_join_completed": False,
        "original_flower_effect_estimated": False,
    }
    return receipt, weights


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-file", type=Path, help="Previously downloaded original BirdFuncDat.txt")
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    if args.input_file:
        blob = args.input_file.read_bytes()
        origin = "local_original_figshare_file"
    else:
        blob, origin = get_original_bytes()
    receipt, weights = audit_bird_diet(blob)
    receipt["source_delivery"] = origin
    args.output_dir.mkdir(parents=True, exist_ok=True)
    weights.to_csv(args.output_dir / "eltontraits_bird_nectar_weights.csv", index=False)
    (args.output_dir / "source_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "PUBLIC_DIET_ONLY_NOT_GUILD_MAP",
                      "n_source_birds": receipt["n_original_bird_records"],
                      "n_nectar_positive": receipt["n_recorded_nectar_positive"],
                      "bird_range_or_floral_effect": False}, sort_keys=True))


if __name__ == "__main__":
    main()
