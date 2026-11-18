import numpy as np
import pytest

from fd.dopi_degree_fixed_visitor_null import (
    bipartite_degree_swap,site_degree_preserving_null,
    audit_original_degree_fixed_networks,
    _degrees,
)


def source_network():
    edges={
        ("Alpha pos","Bombus terrestris"),
        ("Alpha pos","Episyrphus balteatus"),
        ("Alpha neg","Bombus lapidarius"),
        ("Alpha neg","Syrphus ribesii"),
        ("Beta pos","Bombus lapidarius"),
        ("Beta pos","Syrphus ribesii"),
        ("Beta neg","Bombus terrestris"),
        ("Beta neg","Episyrphus balteatus"),
    }
    guild={
        "Bombus terrestris":"Bombus",
        "Bombus lapidarius":"Bombus",
        "Episyrphus balteatus":"Syrphidae_hoverfly",
        "Syrphus ribesii":"Syrphidae_hoverfly",
    }
    return edges,guild


def motifs():
    known={"Alpha pos","Alpha neg","Beta pos","Beta neg"}
    return {
        "open_radial":dict(known=known,joint={"Alpha pos","Beta pos"},axis_count=1),
        "tubular_bell":dict(known=known,joint={"Alpha pos","Beta pos"},axis_count=1),
        "zygomorphic":dict(known=known,joint={"Alpha pos","Beta pos"},axis_count=1),
        "deep_tube":dict(known=known,joint={"Alpha pos","Beta pos"},axis_count=1),
    }


def test_degree_swaps_exactly_preserve_real_plant_and_real_insect_species_degrees():
    original,guild=source_network()
    before=_degrees(original)
    rng=np.random.default_rng(20261118)
    randomized,n=bipartite_degree_swap(original,rng,attempts=1000)
    assert n>0
    assert len(randomized)==len(original)
    assert _degrees(randomized)==before
    assert randomized!=original or n>1
    assert set(v for _,v in randomized)==set(guild)
    assert len(original)==8


def test_two_chain_null_uses_actual_species_guild_labels_and_limits():
    edges,guild=source_network()
    out=site_degree_preserving_null(edges,guild,motifs(),seed="synthetic",draws=60)
    assert out["status"]=="ORIGINAL_SOURCE_DEGREE_FIXED_SITE_PROFILES"
    x=out["motif_profiles"]["open_radial"]
    assert x["original_observed_guild_richness_difference"]==pytest.approx(0)
    assert x["status"]=="DESCRIPTIVE_FIXED_BOTH_SIDES_DEGREES_NETWORK_NULL"
    assert min(x["n_accepted_degree_swaps_each_chain"])>=10
    assert np.isfinite(x["excess_over_network_degree_null"])


def test_full_bipartite_network_cannot_randomize_and_does_not_mean_zero_effect():
    plants=["Alpha pos","Beta pos","Alpha neg","Beta neg"]
    visitors=["Bombus terrestris","Episyrphus balteatus","Syrphus ribesii"]
    guild={"Bombus terrestris":"Bombus",
           "Episyrphus balteatus":"Syrphidae_hoverfly",
           "Syrphus ribesii":"Syrphidae_hoverfly"}
    edges={(p,v) for p in plants for v in visitors}
    z=site_degree_preserving_null(edges,guild,motifs(),seed="full",draws=40)
    assert z["status"]=="ORIGINAL_SOURCE_DEGREE_FIXED_SITE_PROFILES"
    x=z["motif_profiles"]["open_radial"]
    assert x["status"]=="NOT_IDENTIFIABLE_SOURCE_NETWORK_NO_SWAPPABLE_OPPORTUNITY"
    assert "excess_over_network_degree_null" not in x


def test_real_network_null_rejects_non_real_taxonomy_or_undefined_specified_guild():
    edges,guild=source_network()
    with pytest.raises(ValueError,match="guild identity"):
        site_degree_preserving_null(edges,{"Bombus terrestris":"Bombus"},
                                    motifs(),seed="missing",draws=40)
    with pytest.raises(ValueError,match="snapshots"):
        site_degree_preserving_null(edges,guild,motifs(),seed="small",draws=1)


def test_twelve_independent_source_studies_and_no_island_visitation_claim():
    edges,guild=source_network()
    idx=dict(
        private_positive_edges={},
        documented_guild_absence_not_inferred=True,
        n_original_source_studies_qualified=12,
        source_row_counts=dict(original_DoPI_v2_source_rows=676900),
    )
    for i in range(12):
        idx["private_positive_edges"][(f"study{i}",("published_site","53.123","-1.456"),"2025","6")]={
            (p,v,guild[v]) for p,v in edges}
    report=audit_original_degree_fixed_networks(
        idx,motifs=motifs(),draws=199,site_draws=40)
    assert report["original_GIFT_botanical_island_population"]==4505
    assert report["n_original_published_site_months"]==12
    assert len(report["panels"])==4
    for p in report["panels"]:
        assert p["n_source_publications_with_degree_fixed_null"]==12
        assert p["status"]=="DESCRIPTIVE_SOURCE_VISITOR_DEGREE_FIXED_FLORAL_GUILD_ASSOCIATION"
        assert np.isfinite(p["residual_over_degree_fixed_network_null"])
    assert report["caveats"]["each_insect_species_positive_host_degree_preserved"]


def test_no_true_morphological_contrast_is_not_a_guild_null():
    edges,guild=source_network()
    missing=dict(motifs())
    missing["open_radial"]=dict(known=set(motifs()["open_radial"]["known"]),
                                joint=set(),axis_count=1)
    z=site_degree_preserving_null(edges,guild,missing,seed="missing",draws=40)
    assert "open_radial" not in z["motif_profiles"]
