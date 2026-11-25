"""Synthetic strict original EltonTraits tests; no botanical or animal IDs."""
import pandas as pd
import pytest
from fd.eltontraits_bird_nectar import audit_bird_diet, get_original_bytes


def source(rows):
    return pd.DataFrame(rows).to_csv(index=False, sep="\t").encode("utf-8")


def test_percentage_scale_without_eBird():
    summary, weight = audit_bird_diet(
        source([{"Scientific": "Bird_a", "Diet.Nect": "75.0"},
                {"Scientific": "Bird_b", "Diet.Nect": "0"}]),
        enforce_species_range=False)
    assert weight["original_scientific_name"].tolist() == ["Bird a", "Bird b"]
    assert weight["diet_nectar_fraction"].tolist() == [.75, 0.]
    assert summary["n_recorded_nectar_positive"] == 1
    assert summary["has_eBird_taxonomy_or_occurrences"] is False
    assert summary["island_source_animal_join_completed"] is False


@pytest.mark.parametrize("unknown", ["", "NaN", "0.5x"])
def test_missing_nectar_never_becomes_zero(unknown):
    summary, weights = audit_bird_diet(
        source([{"Scientific": "Bird_a", "Diet.Nect": unknown}]),
        enforce_species_range=False)
    assert summary["n_missing_diet_original_rows"] == 1
    assert len(weights) == 0
    assert summary["n_recorded_nectar_zero"] == 0


@pytest.mark.parametrize("invalid", ["101", "-2"])
def test_out_of_bounds_diet_fails_closed(invalid):
    with pytest.raises(ValueError, match="Invalid/non-percentage"):
        audit_bird_diet(
            source([{"Scientific": "Bird_a", "Diet.Nect": invalid}]),
            enforce_species_range=False)


def test_real_original_archive_hyphenated_header_is_accepted():
    summary, weights = audit_bird_diet(
        source([{"Scientific": "Bird_1", "Diet-Nect": "25"},
                {"Scientific": "Bird_2", "Diet-Nect": "0"}]),
        enforce_species_range=False,
    )
    assert summary["source_nectar_column"] == "Diet-Nect"
    assert weights.diet_nectar_fraction.tolist() == [.25, 0.]


def test_ambiguous_dual_nectar_column_is_rejected():
    with pytest.raises(ValueError, match="missing/ambiguous"):
        audit_bird_diet(
            source([{"Scientific":"Bird_1","Diet-Nect":25,"Diet.Nect":25}]),
            enforce_species_range=False)


def test_wrong_source_header_and_size_are_rejected():
    with pytest.raises(ValueError, match="header"):
        audit_bird_diet(source([{"Species": "Bird_a", "nectar": 20}]),
                        enforce_species_range=False)
    with pytest.raises(ValueError, match="population"):
        audit_bird_diet(source([{"Scientific": "Bird_a", "Diet.Nect": 20}]))


def test_exact_duplicate_is_collapsed_and_conflicting_values_quarantined():
    rows = [{"Scientific": "Bird_a", "Diet.Nect": 20},
            {"Scientific": "Bird_a", "Diet.Nect": 20}]
    summary, weights = audit_bird_diet(source(rows), enforce_species_range=False)
    assert summary["n_original_rows_with_nonunique_names"] == 2
    assert summary["n_original_names_with_conflicting_values_quarantined"] == 0
    assert weights["diet_nectar_fraction"].tolist() == [.2]

    conflicting = [{"Scientific": "Bird_a", "Diet.Nect": 20},
                   {"Scientific": "Bird_a", "Diet.Nect": 40}]
    summary, weights = audit_bird_diet(source(conflicting), enforce_species_range=False)
    assert summary["n_original_names_with_conflicting_values_quarantined"] == 1
    assert len(weights) == 0


def test_unapproved_download_url_is_rejected_before_network():
    with pytest.raises(ValueError, match="Unapproved"):
        get_original_bytes(("https://example.com/BirdFuncDat.txt",))
