"""Synthetic checks that probability scores are proper for categorical flowers."""
import numpy as np
import pandas as pd
import pytest

from fd.proper_scoring import inspect_kernel, heldout_probability_scores


class TinyModel:
    n_profiles = 2

    def probabilities(self, island, genus):
        assert island == "island" and genus == "Gen"
        return np.array([
            [.5, .5],        # no-genus prior
            [.5, .5],        # unused local-mainland prior
            [.25, .75],      # global genus
            [.1, .9],        # local source genus
        ])


def test_reference_similarity_is_psd_and_rejects_indefinite_kernel():
    k = inspect_kernel(np.array([[0., 1.], [1., 0.]]))
    assert k["rank_positive_1e8"] == 2
    assert k["minimum_eigenvalue"] == pytest.approx(1.)
    assert np.allclose(k["kernel"], np.eye(2))
    # A deliberately invalid similarity has negative eigenvalues.
    d = np.array([[0., 0., 1.], [0., 0., 0.], [1., 0., 0.]])
    with pytest.raises(ValueError, match="positive semidefinite"):
        inspect_kernel(d)
    with pytest.raises(ValueError):
        inspect_kernel(np.array([[0., 2.], [2., 0.]]))


def test_brier_equals_kernel_score_for_discrete_identity_similarity():
    part = pd.DataFrame({
        "accepted_species": ["Gen a", "Gen b"],
        "genus": ["Gen", "Gen"], "profile_id": [1, 1],
    })
    result = heldout_probability_scores(part, TinyModel(), "island", np.eye(2))
    assert result["global_prior_brier"] == pytest.approx(.5)
    assert result["global_genus_brier"] == pytest.approx(.125)
    assert result["local_genus_brier"] == pytest.approx(.02)
    for model in ("global_prior", "global_genus", "local_genus"):
        assert result[model+"_kernel"] == pytest.approx(result[model+"_brier"])
        assert result[model+"_logloss"] == pytest.approx(
            -np.log({"global_prior": .5, "global_genus": .75, "local_genus": .9}[model]))
    assert result["local_genus_kernel"] < result["global_genus_kernel"]
    with pytest.raises(ValueError):
        heldout_probability_scores(
            pd.concat([part, part.iloc[[0]]], ignore_index=True),
            TinyModel(), "island", np.eye(2))


def test_positive_semidefinite_nonidentity_kernel_proper_distance_score():
    part = pd.DataFrame({
        "accepted_species": ["Gen first", "Gen second"],
        "genus": ["Gen", "Gen"], "profile_id": [0, 1],
    })
    k = np.array([[1., .5], [.5, 1.]])
    result = heldout_probability_scores(part, TinyModel(), "island", k)
    for label, p in (
        ("global_prior", np.array([.5, .5])),
        ("global_genus", np.array([.25, .75])),
        ("local_genus", np.array([.1, .9])),
    ):
        expected_kernel = 1 + p@k@p - 2*np.mean([(k@p)[0], (k@p)[1]])
        assert result[label+"_kernel"] == pytest.approx(expected_kernel)
        assert result[label+"_trait_distance_loss"] == pytest.approx(
            np.mean([1-(k@p)[0], 1-(k@p)[1]]))
