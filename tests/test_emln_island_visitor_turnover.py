"""Synthetic contract tests for EMLN source and site-matched visit turnover.

Never add fake pollinator taxa to original empirical results.
"""
import pandas as pd
import pytest

from fd.emln_island_visitor_turnover import (
    _pivot, _jaccard, original_emln_source,
    canary_same_plant_partner_turnover, mauritius_native_partner_rewiring
)


def _small_source():
    node=pd.DataFrame([
      ("1","node_name","NativePlant"),
      ("1","type","plant"),
      ("2","node_name","BeeA"),
      ("2","type","pollinator"),
    ],columns=["node_id","attribute","value"])
    edge=pd.DataFrame([
      ("1","node_from","NativePlant"),("1","node_to","BeeA"),
      ("1","layer_from","loc1"),("1","layer_to","loc1"),
      ("1","weight","4"),("1","type","pollination"),
    ],columns=["interaction_id","attribute","value"])
    layer=pd.DataFrame([
      ("1","name","loc1"),("1","type","time"),
      ("1","type","perturbation"),
    ],columns=["layer","attribute","value"])
    return node,edge,layer


def test_original_emln_long_to_wide_repeated_layer_type_is_legitimate():
    n,e,l=_small_source()
    nodes,edges,layers=original_emln_source(n,e,l,expected=1)
    assert len(nodes)==2 and len(edges)==1 and len(layers)==1
    assert edges.iloc[0].plant=="NativePlant"
    assert edges.iloc[0].visitor=="BeeA"
    assert layers.iloc[0].type=="time|perturbation"
    assert edges.iloc[0].weight==4


def test_not_treating_duplicate_source_taxon_as_new_population():
    n,e,l=_small_source()
    e=pd.concat([e,e.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="Conflicting"):
        original_emln_source(n,e,l,expected=1)


def test_jaccard_respects_observed_partner_identity():
    assert _jaccard({"A","B"},{"B","C"})==pytest.approx(1/3)
    assert _jaccard({"A"},{"X"})==0
    assert _jaccard({"B"},{"B"})==1


def test_one_network_never_promoted_to_5_independent_islands():
    n,e,l=_small_source()
    nodes,edges,layers=original_emln_source(n,e,l,expected=1)
    with pytest.raises(ValueError,match="source support"):
        canary_same_plant_partner_turnover(nodes,edges,layers)
    with pytest.raises(ValueError,match="source panel"):
        mauritius_native_partner_rewiring(nodes,edges,layers)
