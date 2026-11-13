"""True original branch PhyloSor beta and floral support Jaccard are separate metrics."""
import numpy as np
import pandas as pd
import pytest

from fd.phylo_beta_functional_beta_macro import (
    OriginalTreeSubtree,functional_reference_beta,_weighted_slope
)


def original_tiny_tree():
    edge=np.array([[7,5],[7,6],[5,1],[5,2],[6,3],[6,4]],
                  dtype=np.int64)
    lengths=np.array([2.,3.,1.,1.5,2.,1.],dtype=float)
    return OriginalTreeSubtree(edge,lengths,np.array([1,2,3,4]))


def test_original_rooted_phylo_beta_uses_shared_branch_length_not_kde():
    p=original_tiny_tree()
    # [1,3] and [2,4] use both major backbone branches, but
    # exact tips differ. Shared original branch length 2+3=5.
    beta=p.phylosor_beta([1,3],[2,4])
    assert beta==pytest.approx(1-2*5/15.5)
    assert p.phylosor_beta([1,2],[3,4])==pytest.approx(1.)
    assert p.phylosor_beta([1,3],[1,3])==pytest.approx(0.)
    with pytest.raises(ValueError,match="Non-original"):
        p.phylosor_beta([1,5],[2,3])


def test_source_bitset_jaccard_functional_beta_is_not_profile_count():
    assert functional_reference_beta(0b1011,0b1101)==pytest.approx(0.5)
    assert functional_reference_beta(0b1001,0b1001)==pytest.approx(0.)
    assert functional_reference_beta(0b0011,0b1100)==pytest.approx(1.)
    with pytest.raises(ValueError,match="No actually observed"):
        functional_reference_beta(0,0)


def synthetic_original_region():
    rng=np.random.default_rng(13)
    n=45
    x=np.linspace(.13,.78,n)+rng.normal(0,.018,n)
    g=rng.uniform(.1,.8,n)
    f=rng.uniform(.15,.72,n)
    source=rng.uniform(3,8,n)
    area=rng.uniform(2,9,n)
    dist=rng.uniform(1,7,n)
    data=pd.DataFrame(dict(
        phylosor_beta=x,
        functional_support_beta=.38*x+.1*g+.025*f+rng.normal(0,.025,n),
        original_genus_jaccard_turnover=g,
        original_colxr_candidate_family_jaccard_turnover=f,
        log_source_exact_tip_candidates=source,
        log_island_area_km2=area,
        log_distance_to_continent_km=dist,
        spatial_block=[f"block{i%15}" for i in range(n)],
        source_signature=[f"source{i%12}" for i in range(n)],
    ))
    return data


def test_phylo_functional_correlation_cluster_weighted_not_repetition_pseudoreplication():
    z=synthetic_original_region()
    for cluster in ("spatial_block","source_signature"):
        raw=_weighted_slope(z,"phylosor_beta","functional_support_beta",
                            cluster,adjusted=False)
        adjusted=_weighted_slope(z,"phylosor_beta","functional_support_beta",
                                 cluster,adjusted=True)
        assert raw>0
        assert adjusted>0
        assert np.isfinite(adjusted)


def test_not_enough_independent_blocks_is_not_a_zero_slope():
    z=synthetic_original_region()
    z.spatial_block="one_block"
    with pytest.raises(ValueError,match="Insufficient"):
        _weighted_slope(z,"phylosor_beta","functional_support_beta",
                        "spatial_block",adjusted=False)
