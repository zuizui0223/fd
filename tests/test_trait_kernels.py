"""Synthetic mathematical checks: additive floral kernel misses true pairings."""
import numpy as np
import pytest

from fd.core import distance
from fd.trait_kernels import (
    kernel_quality, phenotype_kernels, quadratic_kernel_score,
)


def _four_combinations():
    return {
        0: (frozenset({"red"}), frozenset({"tube"})),
        1: (frozenset({"red"}), frozenset({"open"})),
        2: (frozenset({"blue"}), frozenset({"tube"})),
        3: (frozenset({"blue"}), frozenset({"open"})),
    }


def test_additive_kernel_is_exact_original_colour_form_d():
    profiles = _four_combinations()
    kernels, audit = kernel_quality(profiles)
    for i in profiles:
        for j in profiles:
            assert kernels["additive"][i, j] == pytest.approx(
                1-distance(profiles[i], profiles[j]))
    assert np.allclose(kernels["joint_product"], np.eye(4))
    assert all(audit[k]["minimum_eigenvalue"] >= -1e-8
               for k in kernels)
    assert audit["joint_product"]["rank_positive_1e8"] == 4


def test_identical_colour_and_form_marginals_can_hide_opposite_pairings():
    k = phenotype_kernels(_four_combinations())
    p = np.array([.5, 0., 0., .5])  # red/tube + blue/open
    q = np.array([0., .5, .5, 0.])  # red/open + blue/tube
    delta = p-q
    for key in ("colour", "form", "additive"):
        assert delta @ k[key] @ delta == pytest.approx(0., abs=1e-12)
    # The joint product detects cross-trait COUPLING missed by additive D.
    assert delta @ k["joint_product"] @ delta == pytest.approx(1.)
    expected = sum(q[y]*(
        quadratic_kernel_score(p, y, k["joint_product"]) -
        quadratic_kernel_score(q, y, k["joint_product"]))
        for y in range(4))
    assert expected == pytest.approx(1.)


def test_multistate_profiles_remain_intact_recorded_sets():
    profiles = {
        0: (frozenset({"white", "blue"}), frozenset({"tubular"})),
        1: (frozenset({"white"}), frozenset({"tubular"})),
        2: (frozenset({"blue"}), frozenset({"open"})),
    }
    k = phenotype_kernels(profiles)
    assert k["colour"][0, 1] == pytest.approx(.5)
    assert k["form"][0, 1] == pytest.approx(1.)
    assert k["joint_product"][0, 1] == pytest.approx(.5)
    assert len(k["additive"]) == len(profiles) == 3
    assert k["joint_product"][0, 2] == pytest.approx(0.)


def test_probability_score_checks_shapes_and_names():
    k = np.eye(3)
    p = np.array([.2,.3,.5])
    assert quadratic_kernel_score(p, 2, k) == pytest.approx(
        1+sum(p**2)-2*p[2])
    with pytest.raises(ValueError):
        quadratic_kernel_score(np.array([.2,.2,.5]), 2, k)
    with pytest.raises(ValueError):
        quadratic_kernel_score(p, 3, k)
    with pytest.raises(ValueError):
        quadratic_kernel_score(p, 1, np.eye(4))
