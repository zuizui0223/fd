"""Synthetic-only N2 five-geography-fold prediction and cluster uncertainty."""
import numpy as np
import pandas as pd
import pytest

from fd.n2_crossfit_macro import (
    METRICS, block_bootstrap_mean, summarize,
)


def test_block_bootstrap_clusters_not_individual_islands_and_is_deterministic():
    # 10 dependent islands in one block + one island in each other block.
    groups = ["A"] * 10 + ["B", "C", "D", "E", "F"]
    values = [2.] * 10 + [-1.] * 5
    data = pd.DataFrame({"spatial_block": groups, "delta": values})
    mu, low, high = block_bootstrap_mean(
        data, response="delta", draws=999, seed=42)
    assert mu == pytest.approx(1.)
    assert low < mu < high
    assert (mu, low, high) == block_bootstrap_mean(
        data, response="delta", draws=999, seed=42)
    with pytest.raises(ValueError):
        block_bootstrap_mean(
            data.assign(spatial_block=None), response="delta", draws=999)
    with pytest.raises(ValueError):
        block_bootstrap_mean(
            data.loc[data.spatial_block.isin(["A", "B"])],
            response="delta", draws=999)


def test_grouped_summary_keeps_the_five_fold_source_contrast_sign():
    rows = []
    for scope in ("all", "direct"):
        for block in range(6):
            rows.append({
                "scope": scope, "analysis_regime": "tropical",
                "spatial_block": f"B{block}", "island_id": f"secret{block}",
                "global_coverage": .8, "local_coverage": .3,
                "prior_logloss": 6., "global_logloss": 5.,
                "local_logloss": 5.1,
                "genus_gain": 1., "source_gain": -.1,
                "source_distance_gain": .02,
                "genus_brier_gain": .3, "source_brier_gain": .004,
                "genus_kernel_gain": .5, "source_kernel_gain": .01,
            })
    output = summarize(rows, draws=299, min_cell=5)
    assert len(output) == 4  # one global and one tropical row per scope
    assert all(r["mean_source_gain"] == pytest.approx(-.1) for r in output)
    assert all(r["source_gain_block_q975"] < 0 for r in output)
    assert all(r["genus_gain_block_q025"] > 0 for r in output)
    assert all(r["mean_source_distance_gain"] == pytest.approx(.02)
               for r in output)
    assert all(r["mean_source_brier_gain"] == pytest.approx(.004) for r in output)
    assert all(r["mean_source_kernel_gain"] == pytest.approx(.01) for r in output)
    assert all(r["source_kernel_gain_block_q025"] > 0 for r in output)
    assert not any("island_id" in row or "spatial_block" in row for row in output)


def test_scant_spatial_replication_suppresses_predictive_claim():
    rows = [{
        "scope": "all", "analysis_regime": "southern_extratropical",
        "spatial_block": "ONE", "global_coverage": .1,
        "local_coverage": .05, "prior_logloss": 5.0,
        "global_logloss": 4.0, "local_logloss": 4.5,
        "genus_gain": 1.0, "source_gain": -.5,
        "source_distance_gain": .01,
    } for _ in range(8)]
    out = summarize(rows, draws=99)
    assert len(out) == 2
    assert all(r["status"] == "insufficient_independent_geographic_support"
               for r in out)
    assert all("mean_genus_gain" not in r for r in out)
    assert METRICS == (
        "genus_gain", "source_gain", "source_distance_gain",
        "genus_brier_gain", "source_brier_gain",
        "genus_kernel_gain", "source_kernel_gain",
    )
