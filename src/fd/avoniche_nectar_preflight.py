"""Independent public AVONICHE BirdLife-2020 nectar-use source preflight.

This does NOT create island availability, pollinator visits, or guild absence.
No eBird-derived AVONICHE file is used.
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

DATASET_DOI = "10.5061/dryad.1zcrjdg56"
BIRDLIFE_FILE = "Avoniche1_Birdlife2020.csv"
DIRECT_URL = "https://datadryad.org/downloads/file_stream/4513904"
API_URL = "https://datadryad.org/api/v2/files/4513904/download"
DIET = ("Pa", "Pt", "Ne", "Se", "Fr", "In", "Ap", "Vt", "Ca")
EXPECTED_ROWS = 10999
MAX_BYTES = 8_000_000


def download_public_csv(*, urls: tuple[str, ...] = (DIRECT_URL, API_URL)) -> tuple[bytes, str]:
    """Use only the exact BirdLife2020 public original file, not another taxonomy."""
    errors = []
    for url in urls:
        if not url.startswith("https://datadryad.org/"):
            raise ValueError("Only the official public Dryad source is allowed")
        try:
            req = Request(url, headers={"User-Agent": "fd-avoniche-research-preflight/1.0",
                                        "Accept": "text/csv,application/octet-stream,*/*"})
            with urlopen(req, timeout=40) as handle:
                blob = handle.read(MAX_BYTES + 1)
            if len(blob) > MAX_BYTES:
                raise ValueError("Download exceeds public CSV size cap")
            if not blob or b"<html" in blob[:200].lower():
                raise ValueError("Downloaded HTML/empty content, not a CSV")
            return blob, url
        except Exception as exc:  # one alternate official URL, then fail closed
            errors.append(type(exc).__name__ + ": " + str(exc)[:150])
    raise RuntimeError("AVONICHE official public download unavailable: " + "; ".join(errors))


def _choose_column(columns: list[str], *, explicit: str | None,
                   guesses: tuple[str, ...], label: str) -> str:
    if explicit is not None:
        if explicit not in columns:
            raise ValueError(f"Requested {label} column {explicit!r} is missing")
        return explicit
    norm = {str(c).strip().casefold().replace(" ", "_"): c for c in columns}
    for option in guesses:
        if option.casefold() in norm:
            return norm[option.casefold()]
    raise ValueError(f"Unverified {label} schema. Columns: {columns!r}")


def inspect_public_birdlife_csv(blob: bytes, *, expected_rows: int = EXPECTED_ROWS,
                                species_col: str | None = None,
                                nectar_col: str | None = None) -> tuple[dict, pd.DataFrame]:
    """Return auditable source receipt and bird nectar weights, NOT geographic guild maps.

    Percentage scale is accepted only when all nine dietary niches independently
    certify a sum of ~100; otherwise a manual source review is required.
    """
    sha = hashlib.sha256(blob).hexdigest()
    try:
        frame = pd.read_csv(io.BytesIO(blob), dtype=str, encoding="utf-8-sig",
                            keep_default_na=False)
    except (UnicodeError, pd.errors.ParserError) as exc:
        raise ValueError("AVONICHE BirdLife CSV cannot be parsed") from exc
    if int(len(frame)) != expected_rows:
        raise ValueError(f"AVONICHE BirdLife expected {expected_rows} species, got {len(frame)}")
    if frame.columns.duplicated().any():
        raise ValueError("Ambiguous duplicate AVONICHE columns")
    species = _choose_column(list(frame), explicit=species_col,
                             guesses=("Scientific_Name", "ScientificName", "Scientific name",
                                      "Species_Name", "species", "binomial", "Sp1", "Taxon"),
                             label="species")
    nectar = _choose_column(list(frame), explicit=nectar_col,
                            guesses=("Ne", "Diet_Nectar", "Diet.Nectar", "Nectar"),
                            label="nectar diet")
    names = frame[species].astype(str).str.strip()
    if names.eq("").any() or names.duplicated().any():
        raise ValueError("Empty or duplicated accepted bird names: manual taxonomic review required")
    values = pd.to_numeric(frame[nectar].replace("", np.nan), errors="coerce")
    if values.isna().any() or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Nonfinite, missing or negative AVONICHE nectar diet values")

    cols = {short: next((col for col in frame if str(col).casefold() == short.casefold()), None)
            for short in DIET}
    diet_all = all(col is not None for col in cols.values())
    if diet_all:
        diet_frame = frame[[cols[short] for short in DIET]].apply(pd.to_numeric, errors="coerce")
        if diet_frame.isna().any().any():
            raise ValueError("Nine diet components contain missing/non-numeric entries")
        sums = diet_frame.sum(axis=1)
        if (sums.sub(1.0).abs() <= .06).mean() >= .99:
            scale = 1
        elif (sums.sub(100.0).abs() <= 6.0).mean() >= .99:
            scale = 100
        else:
            raise ValueError("Nine AVONICHE diets fail 1/100 scale certification")
    elif float(values.max()) <= 1:
        scale = 1
    else:
        raise ValueError("Nectar score >1 without independent nine-diet scale certification")
    fraction = values / scale
    if (fraction > 1 + 1e-9).any():
        raise ValueError("Nectar fractions exceed 1 after unit validation")
    weights = pd.DataFrame({"birdlife2020_scientific_name": names,
                            "nectar_fraction": fraction.astype(float)})
    receipt = {
        "contract": "fd_public_avoniche_birdlife2020_nectar_v1",
        "role": "external_animal_diet_not_geographic_guild_presence",
        "doi": DATASET_DOI,
        "file": BIRDLIFE_FILE,
        "sha256": sha,
        "bytes": len(blob),
        "n_birdlife2020_species": int(len(weights)),
        "n_recorded_nectar_use_positive": int(fraction.gt(0).sum()),
        "n_recorded_nectar_use_zero": int(fraction.eq(0).sum()),
        "nectar_fraction_sum_across_taxa": float(fraction.sum()),
        "input_species_column": species,
        "input_nectar_column": nectar,
        "diet_scale": scale,
        "all_nine_diet_components_verified": bool(diet_all),
        "uses_ebird_taxonomy_or_occurrence": False,
        "includes_bird_distribution_or_pollinator_effectiveness": False,
        "ready_for_birdlife_range_polygon_join": True,
        "range_polygon_join_completed": False,
        "sphingidae_map_join_completed": False,
        "floral_guild_test_executed": False,
    }
    return receipt, weights


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path,
                        help="Already acquired original Dryad BirdLife2020 CSV")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--species-column")
    parser.add_argument("--nectar-column")
    args = parser.parse_args()
    if args.input_csv:
        blob = args.input_csv.read_bytes()
        origin = "verified_local_public_csv"
    else:
        blob, origin = download_public_csv()
    receipt, weights = inspect_public_birdlife_csv(
        blob, species_col=args.species_column, nectar_col=args.nectar_column)
    receipt["public_source_delivery"] = origin
    args.output_dir.mkdir(parents=True, exist_ok=True)
    weights.to_csv(args.output_dir / "birdlife2020_nectar_weights.csv", index=False)
    (args.output_dir / "public_source_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    print(json.dumps({"status": "PUBLIC_BIRD_DIET_ONLY_VERIFIED",
                      "n_birds": receipt["n_birdlife2020_species"],
                      "n_nectar_positive": receipt["n_recorded_nectar_use_positive"],
                      "map_join_complete": False}, sort_keys=True))


if __name__ == "__main__":
    main()
