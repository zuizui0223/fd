"""Pinned-public EuPPollNet flower-availability denominator audit.

This is SOURCE QUALITY/FEASIBILITY, never inferred visitor preference.
The original Flower_counts and Interaction_data join on an explicit
published Flower_data_merger key. Repeated key matches are NOT treated
as independent flower observations or automatically summed/merged.
No original GIFT island data/accepted names are loaded here.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
from io import BytesIO
import numpy as np
import pandas as pd
from .euppollnet_public_preflight import (
    SOURCE_COMMIT,SOURCE_FILES,qualify_original_blob,exact_git_blob_sha,
)

FLOWER_SOURCE_FILE="Flower_counts.csv"
FLOWER_SOURCE_SIZE=20315455
FLOWER_GIT_BLOB="8543e446a53f3b7cffece76e1699489e7169f7fd"
JOIN="Flower_data_merger"


def source_original_flower_receipt(data:bytes)->dict:
    if len(data)!=FLOWER_SOURCE_SIZE or exact_git_blob_sha(data)!=FLOWER_GIT_BLOB:
        raise ValueError("Original public EuPPollNet Flower_counts source changed")
    return dict(original_source_file=FLOWER_SOURCE_FILE,
        original_git_commit=SOURCE_COMMIT,
        original_file_bytes=len(data),
        original_git_blob_sha1=FLOWER_GIT_BLOB,
        original_file_sha256=sha256(data).hexdigest())


def _nonmissing(q:pd.Series)->pd.Series:
    v=q.dropna().astype(str).str.strip()
    return v[v.ne("")&~v.str.lower().isin(
        ("na","nan","null","none","unknown"))]


def audit_public_flower_source(flower_blob:bytes, interaction_blob:bytes,
                              *, allow_synthetic_source:bool=False)->dict:
    """Source-original joint-key presence, no abundance-model invention."""
    if not allow_synthetic_source:
        receipt=source_original_flower_receipt(flower_blob)
        visitor=qualify_original_blob("Interaction_data.csv.gz",interaction_blob)
    else:
        receipt=dict(original_source_file="SYNTHETIC",
                     original_file_bytes=len(flower_blob))
        visitor=dict(original_public_source="SYNTHETIC")
    try:
        f=pd.read_csv(BytesIO(flower_blob),dtype=str,low_memory=False,
            encoding="utf-8",encoding_errors="replace")
        i=pd.read_csv(BytesIO(interaction_blob),compression="gzip",
            dtype=str,low_memory=False,encoding="utf-8",
            encoding_errors="replace",
            usecols=lambda c:c==JOIN)
    except (ValueError,UnicodeError,EOFError) as exc:
        raise ValueError("Original public source CSV impossible to parse") from exc
    if JOIN not in f or JOIN not in i:
        raise ValueError("Published Flower_data_merger key missing")
    fkey=_nonmissing(f[JOIN]);ikey=_nonmissing(i[JOIN])
    set_f,set_i=set(fkey),set(ikey)
    common=set_f&set_i
    freq_f=fkey.value_counts()
    freq_i=ikey.value_counts()
    # Any source join key with multiple flower-count rows requires survey
    # grain validation BEFORE use as botanical exposure/availability.
    ambiguous=set(freq_f.index[freq_f.gt(1)])
    potential={x for x in common if x not in ambiguous}
    # Headers only; never publish original per-plant/site/taxon join key values.
    numeric=[c for c in f.columns if c!=JOIN and any(x in c.lower()
             for x in ("flower_count","flowers_count","flower_n",
                       "n_flowers","flower_abundance","flower_number"))]
    availability={}
    for col in numeric:
        q=pd.to_numeric(f[col],errors="coerce")
        good=np.isfinite(q)
        availability[col]=dict(n_nonnegative_numeric=int((q[good]>=0).sum()),
            n_negative_numeric=int((q[good]<0).sum()),
            n_missing_or_nonnumeric=int((~good).sum()))
    return dict(
        status="PUBLIC_EUPPOLLNET_ORIGINAL_FLOWER_EXPOSURE_KEY_FEASIBILITY",
        source_original_flower_receipt=receipt,
        source_original_visitor_receipt=visitor,
        n_original_flowers_raw_rows=int(len(f)),
        n_original_interaction_raw_rows=int(len(i)),
        n_nonempty_flower_join_rows=int(len(fkey)),
        n_nonempty_interaction_join_rows=int(len(ikey)),
        n_distinct_flower_join_keys=int(len(set_f)),
        n_distinct_interaction_join_keys=int(len(set_i)),
        n_original_distinct_join_keys_shared=int(len(common)),
        n_distinct_flower_join_keys_with_multiple_source_rows=int(len(ambiguous)),
        n_original_join_keys_shared_and_unique_flower_count_rows=int(len(potential)),
        n_original_interaction_rows_on_any_shared_join_key=int(
            ikey.isin(common).sum()),
        n_original_interaction_rows_with_unique_flower_key=int(
            ikey.isin(potential).sum()),
        flower_source_column_names=list(f.columns),
        original_quantitative_flower_abundance_column_candidates=availability,
        conclusions=dict(
            source_species_taxonomy_and_original_GIFT_flora_not_yet_joined=True,
            original_plant_flower_abundance_not_assumed_available_for_every_visit=True,
            flower_source_duplicate_join_keys_not_auto_aggregated=True,
            original_source_codebook_needed_to_interpret_count_units=True,
            valid_source_flower_count_is_not_measured_pollen_deposition=True,
            original_European_study_not_tropical_island_replication=True,
            no_effect_coefficient_or_pvalue_estimated=True,
        ))
