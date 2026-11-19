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
    if "Scientific" not in frame.columns:
        raise ValueError("Original BirdFuncDat header missing Scientific")
    # Figshare's original archived BirdFuncDat.txt uses a literal hyphen;
    # older secondary packages may recode it to a dot when importing to R.
    # Never accept BOTH interpretations silently if future schema drifts.
    nectar_matches = [c for c in ("Diet-Nect", "Diet.Nect") if c in frame.columns]
    if len(nectar_matches) != 1:
        raise ValueError(f"Original BirdFuncDat header missing/ambiguous nectar: {list(frame.columns)[:20]}")
    nectar_col = nectar_matches[0]
    if enforce_species_range and not MIN_SPECIES <= len(frame) <= MAX_SPECIES:
        raise ValueError(f"Unexpected EltonTraits source bird population: {len(frame)}")
    names = frame["Scientific"].astype(str).str.strip().str.replace(
        "_", " ", regex=False)
    raw = frame[nectar_col].astype(str).str.strip()
    numeric = pd.to_numeric(raw.where(raw.ne(""), np.nan), errors="coerce")
    if numeric.lt(0).any() or numeric.gt(100).any() or np.isinf(numeric).any():
        raise ValueError("Invalid/non-percentage Diet.Nect; original value outside [0,100]")
    # Never interpret an undocumented missing value as zero. Duplicated
    # accepted names sometimes arise in original taxonomic archives; identical
    # values may be collapsed, but discordant duplicate name values are
    # quarantined (not averaged or chosen by row order).
    candidate = pd.DataFrame({"name": names, "percent": numeric})
    blank = candidate["name"].eq("")
    unknown = candidate["percent"].isna()
    valid = candidate.loc[~blank & ~unknown].copy()
    duplicated = valid["name"].duplicated(keep=False)
    conflict_names = set(
        valid.loc[duplicated].groupby("name")["percent"].nunique().loc[
            lambda v: v.gt(1)].index)
    clean = valid.loc[~valid["name"].isin(conflict_names)].drop_duplicates(
        subset=["name"], keep="first")
    if len(clean) < MIN_SPECIES and enforce_species_range:
        raise ValueError(
            f"Insufficient uniquely identified source diet species: {len(clean)}"
            f" (raw={len(frame)}, blank={int(blank.sum())},"
            f" missing_diet={int(unknown.sum())},"
            f" duplicated_names={int(duplicated.sum())},"
            f" discordant_names={len(conflict_names)})")
    weights = clean.rename(columns={
        "name": "original_scientific_name",
        "percent": "diet_nectar_fraction"})
    weights["diet_nectar_fraction"] = weights["diet_nectar_fraction"] / 100.0
    receipt = {
        "contract": "fd_original_eltontraits_bird_nectar_v1",
        "source_collection_doi": SOURCE_DOI,
        "source_paper_doi": PAPER_DOI,
        "source_filename": "BirdFuncDat.txt",
        "source_file_id": 5631081,
        "source_nectar_column": nectar_col,
        "original_sha256": digest,
        "n_original_bird_records": int(len(frame)),
        "n_name_unique_valid_diet_species": int(len(clean)),
        "n_blank_original_names": int(blank.sum()),
        "n_missing_diet_original_rows": int(unknown.sum()),
        "n_original_rows_with_nonunique_names": int(duplicated.sum()),
        "n_original_names_with_conflicting_values_quarantined": int(len(conflict_names)),
        "n_original_records_quarantined_or_deduplicated": int(len(frame)-len(clean)),
        "n_recorded_nectar_positive": int(weights.diet_nectar_fraction.gt(0).sum()),
        "n_recorded_nectar_zero": int(weights.diet_nectar_fraction.eq(0).sum()),
        "sum_of_nectar_fraction_by_species": float(weights.diet_nectar_fraction.sum()),
        "mean_nectar_fraction_by_species": float(weights.diet_nectar_fraction.mean()),
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
