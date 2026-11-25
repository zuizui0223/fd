"""Synthetic original-source validity gates; no real species/locations committed."""
import pandas as pd
import pytest
from fd.elton_globi_guild_bridge import audit


def make():
    w = pd.DataFrame({
        "original_scientific_name": ["Fictus alpha","Fictus beta","Fictus gamma"],
        "diet_nectar_fraction": [.8,0,.5],
    })
    cat = pd.DataFrame({
        "channel_id": ["flower_visiting_birds"]*3 + ["lepidoptera"]*2,
        "pollinator_species": ["Fictus alpha","Fictus beta","Fictus unmatched",
                               "Imago moth1","Imago moth2"],
        "pollinator_family": ["Fictidae"]*3 + ["Sphingidae","NotSphingidae"],
        "catalog_tier": ["confirmatory"]*5,
        "n_pollination_rows": [3,0,1,4,0],
        "n_flower_visit_rows": [6,2,1,2,3],
    })
    return w,cat


def test_actual_same_species_diet_and_guild_audit_not_site_association():
    w,c=make()
    r=audit(w,c,expect_source_counts=False)
    assert r["n_exact_name_matched_bird_species"] == 2
    assert r["n_bird_catalog_unmatched_diet_species"] == 1
    assert r["n_matched_bird_diet_nectar_positive"] == 1
    assert r["n_matched_bird_diet_nectar_zero"] == 1
    assert r["n_pollination_claim_positive_diet"] == 1
    assert r["n_no_pollination_claim_zero_diet"] == 1
    assert r["n_globi_confirmatory_sphingidae"] == 1
    assert r["n_sphingidae_with_pollination_claim"] == 1
    assert r["independent_spatial_guild_maps_available"] is False
    assert r["no_eBird_data_used"] is True
    assert "Fictus" not in str(r)


def test_unknown_names_are_not_imputed_from_genus():
    w,c=make()
    c.loc[0,"pollinator_species"]="Fictus unknown_same_genus"
    r=audit(w,c,expect_source_counts=False)
    assert r["n_exact_name_matched_bird_species"] == 1
    assert r["n_bird_catalog_unmatched_diet_species"] == 2


def test_ambiguous_guild_name_and_unknown_diet_fail_closed():
    w,c=make()
    c.loc[1,"pollinator_species"]="Fictus alpha"
    with pytest.raises(ValueError,match="confirmatory species catalog drifted"):
        audit(w,c,expect_source_counts=False)
    w,c=make()
    w.loc[0,"diet_nectar_fraction"]=float("nan")
    with pytest.raises(ValueError,match="Nonfinite"):
        audit(w,c,expect_source_counts=False)


def test_no_same_accepted_species_is_insufficient():
    w,c=make()
    c.loc[:2,"pollinator_species"]=["Other alpha","Other beta","Other gamma"]
    with pytest.raises(ValueError,match="No exact"):
        audit(w,c,expect_source_counts=False)


def test_original_counts_are_pinned_and_quarantined():
    w,c=make()
    with pytest.raises(ValueError,match="denominator"):
        audit(w,c,expect_source_counts=True)
    c.loc[0,"catalog_tier"]="sensitivity"
    with pytest.raises(ValueError,match="confirmatory"):
        audit(w,c,expect_source_counts=False)
