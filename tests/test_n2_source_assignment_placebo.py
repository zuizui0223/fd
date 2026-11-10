"""Synthetic and privacy safeguards for native mainland source-placebo matching.

No original island, taxon, genus or GIFT source identities are loaded.
"""
import numpy as np
import pytest

from fd.n2_source_assignment_placebo import (
    matching_candidates, summarize_placebo, POLICIES,
)


def _catalog():
    return {
        "focal": {
            "region": "tropical", "spatial_block": "block0",
            "source_entities": frozenset(["near1", "near2"]),
            "genera": {"GenA", "GenB"}, "n_donor_names": 100,
        },
        "wrong_regime": {
            "region": "north", "spatial_block": "block4",
            "source_entities": frozenset(["other1"]),
            "genera": {"GenA", "GenB"}, "n_donor_names": 100,
        },
        "same_block": {
            "region": "tropical", "spatial_block": "block0",
            "source_entities": frozenset(["other2"]),
            "genera": {"GenA", "GenB"}, "n_donor_names": 100,
        },
        "same_sources": {
            "region": "tropical", "spatial_block": "block2",
            "source_entities": frozenset(["near1", "near2"]),
            "genera": {"GenA", "GenB"}, "n_donor_names": 100,
        },
        "too_sparse": {
            "region": "tropical", "spatial_block": "block3",
            "source_entities": frozenset(["new1"]),
            "genera": {"GenA"}, "n_donor_names": 120,
        },
        "comparable1": {
            "region": "tropical", "spatial_block": "block1",
            "source_entities": frozenset(["new2", "new3"]),
            "genera": {"GenA", "GenB"}, "n_donor_names": 106,
        },
        "comparable2": {
            "region": "tropical", "spatial_block": "block2",
            "source_entities": frozenset(["new3", "new4"]),
            "genera": {"GenA", "GenB"}, "n_donor_names": 98,
        },
        "donor_size_outlier": {
            "region": "tropical", "spatial_block": "block5",
            "source_entities": frozenset(["distant"]),
            "genera": {"GenA", "GenB"}, "n_donor_names": 400,
        },
    }


def test_two_outcome_blind_coverage_matched_source_controls():
    candidates = _catalog()
    counts = {"GenA": 49, "GenB": 51}
    loose, n_loose = matching_candidates(
        "focal", "tropical", "block0", counts, candidates, POLICIES[0])
    strict, n_strict = matching_candidates(
        "focal", "tropical", "block0", counts, candidates, POLICIES[1])
    assert n_loose == 3 and n_strict == 2
    assert len(loose) == len(strict) == 2
    assert {v["source_island"] for v in strict} == {
        "comparable1", "comparable2"}
    assert all(abs(v["coverage_diff"]) <= .02 for v in strict)
    assert all(.5 <= v["donor_ratio"] <= 2 for v in strict)
    assert not any(v["source_island"] in (
        "focal", "same_block", "same_sources", "wrong_regime", "too_sparse")
        for v in loose)
    assert strict == matching_candidates(
        "focal", "tropical", "block0", counts, candidates, POLICIES[1])[0]
    with pytest.raises(ValueError):
        matching_candidates("focal", "north", "block0", counts,
                            candidates, POLICIES[0])
    with pytest.raises(ValueError):
        matching_candidates("focal", "tropical", "block0", {},
                            candidates, POLICIES[0])


def test_no_controls_if_fewer_than_two_eligible():
    candidates = _catalog()
    candidates.pop("comparable2")
    candidates.pop("donor_size_outlier")
    result, n = matching_candidates(
        "focal", "tropical", "block0",
        {"GenA": 49, "GenB": 51}, candidates, POLICIES[1])
    assert not result and n == 1


def test_summary_clusters_geography_and_hides_private_source_identifiers():
    rows = []
    for policy in POLICIES:
        for i in range(7):
            rows.append({
                "scope": "all", "policy": policy, "fold": i%5,
                "analysis_regime": "tropical",
                "spatial_block": f"secretblock{i}",
                "source_island": f"secretisland{i}",
                "accepted_species": f"secretgenus species{i}",
                "matched_kernel_gain": .002,
                "matched_brier_gain": .001,
                "coverage_difference": .0004,
                "control_source_donor_count_ratio": 1.02,
                "mean_overlap_source_entity_jaccard": .12,
                "n_eligible_control_assignments": 12,
            })
    summary, folds = summarize_placebo(rows, draws=199)
    assert len(summary) == 4  # regional+all for both policies
    for r in summary:
        assert r["mean_kernel_gain"] == pytest.approx(.002)
        assert r["kernel_block_q025"] == pytest.approx(.002)
        assert r["kernel_block_q975"] == pytest.approx(.002)
        assert r["leave_one_block_out_kernel_gain_min"] == pytest.approx(.002)
        assert r["mean_coverage_gain"] == pytest.approx(.0004)
        assert not {"source_island", "accepted_species", "spatial_block"} & r.keys()
    assert all(not {"source_island", "accepted_species", "spatial_block"} & r.keys()
               for r in folds)


def test_small_spatial_cell_suppresses_results():
    rows = [
        dict(scope="direct", policy=POLICIES[1], fold=0,
             analysis_regime="tropical", spatial_block="same",
             matched_kernel_gain=.4, matched_brier_gain=.4,
             coverage_difference=0.,
             control_source_donor_count_ratio=1.,
             mean_overlap_source_entity_jaccard=0.,
             n_eligible_control_assignments=5)
        for _ in range(6)
    ]
    summary, _ = summarize_placebo(rows, draws=99)
    assert len(summary) == 2
    assert all(r["status"] == "insufficient_geographic_support" for r in summary)
    assert all("mean_kernel_gain" not in r for r in summary)
