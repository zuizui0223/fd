"""Synthetic-only nearest geographic floral support/reference safeguards."""
import numpy as np
import pytest

from fd.local_geography_ses import (
    _support_summary, great_circle_distance_km,
    local_pool, nearest_regional_donors,
)
from fd.distance_ses import sample_null


def test_great_circle_haversine_and_coordinate_validation():
    dist = great_circle_distance_km(
        0, 0, np.array([0., 0., 30.]), np.array([0., 90., 0.]))
    assert dist[0] == pytest.approx(0.)
    assert dist[1] == pytest.approx(10007.557, rel=1e-4)
    assert dist[2] == pytest.approx(3335.852, rel=1e-4)
    with pytest.raises(ValueError):
        great_circle_distance_km(95, 0, np.array([0.]), np.array([0.]))
    with pytest.raises(ValueError):
        great_circle_distance_km(0, 0, np.array([91.]), np.array([0.]))


def test_nearest_donor_excludes_focal_and_other_regime_even_when_closer():
    m = {
        key: {"analysis_regime": regime}
        for key, regime in (
            ("focal", "tropical"), ("a", "tropical"),
            ("b", "tropical"), ("c", "tropical"),
            ("d", "northern_midlatitude"),
        )
    }
    geo = {
        "focal": {"island_latitude": "0", "island_longitude": "0"},
        "a": {"island_latitude": "0", "island_longitude": "1"},
        "b": {"island_latitude": "0", "island_longitude": "-1"},
        "c": {"island_latitude": "0", "island_longitude": "3"},
        "d": {"island_latitude": "0", "island_longitude": "0"},
    }
    out, radius = nearest_regional_donors("focal", "tropical", m, geo, k=2)
    assert out == ["a", "b"]  # equal distance, tie resolved by ID
    assert radius == pytest.approx(111.195, abs=.01)
    assert "focal" not in out and "d" not in out
    assert nearest_regional_donors("focal","tropical",m,geo,k=4)[0] == []


def test_local_pool_weights_donor_island_incidence_not_carrier_number():
    donors = ["a", "b", "c"]
    carriers = {
        "a": {0: 100, 1: 1}, "b": {1: 50, 2: 1},
        "c": {0: 1, 2: 1},
    }
    pool, weight = local_pool(donors, carriers, {0: 0, 1: 1, 2: 2})
    assert pool.tolist() == [0, 1, 2]
    assert weight.tolist() == [2, 2, 2]
    # A weighted null cannot manufacture new profiles outside the local pool.
    d=np.array([[0, .2, .8], [.2,0,.5], [.8,.5,0]])
    scores=sample_null(d,pool,3,draws=40,seed=11,weights=weight)
    assert np.unique(scores[:,0]).size == 1
    with pytest.raises(ValueError):
        local_pool(["a", "missing"],carriers,{0:0,1:1,2:2})


def test_small_geographic_support_cell_does_not_release_radius():
    common = dict(scope="all", geological_origin="oceanic",
                  has_neighbors=True, has_full=False, has_equal30=True,
                  n_donor_profiles=70, radius_km=150.)
    rows = [{**common, "island_id": f"private_{i}"} for i in range(4)]
    small = _support_summary(rows)
    assert len(small) == 1 and small[0]["status"] == "small_cell_suppressed"
    assert "median_10th_neighbor_distance_km" not in small[0]
    rows.append({**common, "island_id": "private_4"})
    aggregate = _support_summary(rows)[0]
    assert aggregate["median_10th_neighbor_distance_km"] == pytest.approx(150.)
    assert aggregate["n_n1c_equal30_estimable"] == 5
    assert "island_id" not in aggregate
