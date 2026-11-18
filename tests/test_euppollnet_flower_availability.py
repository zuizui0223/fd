"""Source-pinned original European flower count availability: synthetic control."""
import gzip
import pandas as pd
import pytest
from fd.euppollnet_flower_availability import (
    audit_public_flower_source,source_original_flower_receipt,
)


def _fixture():
    f=pd.DataFrame([
        dict(Flower_data_merger="X",Flower_count="10",
             Plant_accepted_name="A alba"),
        dict(Flower_data_merger="Y",Flower_count="5",
             Plant_accepted_name="B alba"),
        dict(Flower_data_merger="Y",Flower_count="7",
             Plant_accepted_name="C alba"),
        dict(Flower_data_merger="Z",Flower_count=None,
             Plant_accepted_name="D alba"),
    ]).to_csv(index=False).encode()
    i=gzip.compress(pd.DataFrame([
        dict(Flower_data_merger="X",Study_id="One"),
        dict(Flower_data_merger="X",Study_id="One"),
        dict(Flower_data_merger="Y",Study_id="Two"),
        dict(Flower_data_merger="Q",Study_id="Three"),
    ]).to_csv(index=False).encode())
    return f,i


def test_ambiguous_public_flower_count_keys_not_counted_as_verified_availability():
    f,i=_fixture()
    q=audit_public_flower_source(f,i,allow_synthetic_source=True)
    assert q["n_original_flowers_raw_rows"]==4
    assert q["n_original_interaction_raw_rows"]==4
    assert q["n_original_distinct_join_keys_shared"]==2
    assert q["n_distinct_flower_join_keys_with_multiple_source_rows"]==1
    assert q["n_original_join_keys_shared_and_unique_flower_count_rows"]==1
    assert q["n_original_interaction_rows_on_any_shared_join_key"]==3
    assert q["n_original_interaction_rows_with_unique_flower_key"]==2
    assert q["original_quantitative_flower_abundance_column_candidates"][
        "Flower_count"]["n_nonnegative_numeric"]==3
    assert q["conclusions"]["flower_source_duplicate_join_keys_not_auto_aggregated"]


def test_missing_published_flower_merger_key_fails_closed():
    f,i=_fixture()
    wrong=f.replace(b"Flower_data_merger",b"unknown_field")
    with pytest.raises(ValueError,match="Flower_data_merger"):
        audit_public_flower_source(wrong,i,allow_synthetic_source=True)


def test_unpinned_flower_count_bytes_not_accepted():
    f,_=_fixture()
    with pytest.raises(ValueError,match="source changed"):
        source_original_flower_receipt(f)
