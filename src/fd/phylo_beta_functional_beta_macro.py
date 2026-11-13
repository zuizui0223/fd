"""Phylogenetic beta vs functional beta of *the same* original island/source floras.

This is a cross-sectional macroecological source-conditional design, not a
test that pollinator loss caused island floral evolution. Exact published
GBOTB tree tips ONLY (no genus/family grafts), frozen original 928-island
high/medium Direct 3-trait florets, and original five-nearest GIFT candidate
mainland native source floras (all focal island names excluded).

Primary phylogenetic beta: true branch-length PhyloSor beta = 1 -
  2 * shared ROOT-TO-TIP SUBTREE branch length /
  (island subtree branch length + source subtree branch length).
Primary functional beta: Jaccard dissimilarity of the original globally
observed 100-profile reference that lies within ONE original category
mismatch from each sampled flora (finite support, NOT 17-D KDE).
Secondary endpoint: rare-minus-common original profile coverage difference
(island-source) corrected for exact number k of recorded profile types.

Both communities are thinned to the SAME N=15 actual published tree-tip
plant species on 20 repetitions per island. Compute ONE mean predictor and
outcome per physical island; repetitions are not biological independent data.
Within-region source-supply, geographic, genus-/family-composition guards,
whole original 10-degree blocks and first-nearest GIFT-source cluster
bootstrap. Correlation != evolutionary process, causality or real pollinators.

Only regional anonymous summaries are returned. Original species, donor,
genus, family, island name, coordinate and tree tip IDs stay in protected CI.
"""
from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
import hashlib
import numpy as np
import pandas as pd

from .discrete_functional_volume import TRAITS, _resample_group
from .rare_profile_decomposition import RarityReference, _seedstr, genus
from .true_backbone_phylo_source_match import _strict_source_taxon_groups

N_RAREFY=15
N_REPEATS=20
N_PERM=99
N_BS=999
REGION_MIN=20
BLOCK_MIN=10
SOURCE_MIN=10


class OriginalTreeSubtree:
    """Path-compressed *original* published rooted tree; no grafts."""
    def __init__(self,edge:np.ndarray,branch_lengths:np.ndarray,
                 original_tips:np.ndarray):
        edge=np.asarray(edge,dtype=np.int64)
        lengths=np.asarray(branch_lengths,dtype=float)
        tips=np.asarray(original_tips,dtype=np.int64)
        if edge.ndim!=2 or edge.shape[1]!=2 or len(lengths)!=len(edge):
            raise ValueError("Original source rooted tree and branches incompatible")
        if np.any(lengths<=0) or not np.isfinite(lengths).all():
            raise ValueError("Published original tree branch lengths missing/invalid")
        nodes=int(edge.max())+1
        if edge.min()<1 or nodes>400000 or len(edge)!=nodes-2:
            raise ValueError("Original GBOTB edge list not single rooted tree")
        parent=np.zeros(nodes,dtype=np.int32)
        weight=np.zeros(nodes,dtype=float)
        parent[edge[:,1]]=edge[:,0]
        weight[edge[:,1]]=lengths
        if np.any(np.bincount(edge[:,1],minlength=nodes)[1:]>1):
            raise ValueError("Source rooted tree has nodes with multiple parents")
        root=np.setdiff1d(np.unique(edge[:,0]),np.unique(edge[:,1]))
        if len(root)!=1:
            raise ValueError("Source root not unique")
        self.root=int(root[0])
        if len(set(tips))!=len(tips) or np.any(tips<1):
            raise ValueError("Original phylogenetic tip species mapping not unique")
        self.original_tip_paths={}
        for tip in tips:
            path=[]
            current=int(tip)
            seen=set()
            while current!=self.root:
                if current in seen or current<1 or current>=nodes or parent[current]==0:
                    raise ValueError("Original plant tree branch path invalid")
                seen.add(current)
                path.append(current)  # each node denotes its incoming edge
                current=int(parent[current])
            self.original_tip_paths[int(tip)]=frozenset(path)
        self.edge_lengths=weight

    def branch_set(self,tips:list[int]|np.ndarray)->set[int]:
        nodes=set()
        for tip in tips:
            if int(tip) not in self.original_tip_paths:
                raise ValueError("Non-original (grafted or absent) phylogenetic tip used")
            nodes.update(self.original_tip_paths[int(tip)])
        return nodes

    def phylosor_beta(self,tips_a:list[int]|np.ndarray,
                      tips_b:list[int]|np.ndarray)->float:
        a=self.branch_set(tips_a)
        b=self.branch_set(tips_b)
        la=float(sum(self.edge_lengths[list(a)]))
        lb=float(sum(self.edge_lengths[list(b)]))
        if la<=0 or lb<=0:
            raise ValueError("Published source flora has zero evolutionary branch length")
        shared=float(sum(self.edge_lengths[list(a&b)]))
        beta=1-2*shared/(la+lb)
        if beta<(-1e-9) or beta>1+1e-9:
            raise ValueError("PhyloSor original dated branch normalized beta invalid")
        return float(np.clip(beta,0,1))


def functional_reference_beta(mask_a:int,mask_b:int)->float:
    union=mask_a|mask_b
    if union==0:
        raise ValueError("No actually observed floral reference profile support")
    return 1-float((mask_a&mask_b).bit_count()/union.bit_count())


def _phylo_functional_draw(
    reference:RarityReference,tree:OriginalTreeSubtree,local:list[str],
    donor:list[str],tip_map:dict[str,int],taxon_profiles:dict[str,tuple],
    original_family:dict[str,str],
):
    if len(local)!=N_RAREFY or len(donor)!=N_RAREFY:
        raise ValueError("Both source floras must have exactly 15 original taxa")
    if set(local)&set(donor):
        raise ValueError("Focal accepted island species not disjoint from GIFT source")
    if any(x not in tip_map or x not in taxon_profiles for x in local+donor):
        raise ValueError("One source floral species not an ORIGINAL published dated tree tip")
    pA=[taxon_profiles[x][0] for x in local]
    pB=[taxon_profiles[x][0] for x in donor]
    a={z for z in pA}
    b={z for z in pB}
    ma=0
    mb=0
    for t in a:ma|=reference.volume.masks[t][1]
    for t in b:mb|=reference.volume.masks[t][1]
    functional=functional_reference_beta(ma,mb)
    direct_profile_jaccard=1-len(a&b)/len(a|b)
    ea=np.array([(ma&reference.bins[c]).bit_count()/reference.n_bins[c]
                for c in ("rare","middle","common")],dtype=float)
    eb=np.array([(mb&reference.bins[c]).bit_count()/reference.n_bins[c]
                for c in ("rare","middle","common")],dtype=float)
    null_a=reference.expected_uniform[len(a)]
    null_b=reference.expected_uniform[len(b)]
    signed=float((ea[0]-ea[2])-(eb[0]-eb[2]))
    k_signed=signed-float(
        (null_a[0]-null_a[2])-(null_b[0]-null_b[2]))
    phy=tree.phylosor_beta([tip_map[x] for x in local],
                           [tip_map[x] for x in donor])
    gensA={genus(x) for x in local}
    gensB={genus(x) for x in donor}
    famA={original_family[x] for x in local if original_family.get(x)}
    famB={original_family[x] for x in donor if original_family.get(x)}
    return dict(
        phylosor_beta=phy,
        functional_support_beta=functional,
        exact_profile_jaccard_beta=float(direct_profile_jaccard),
        signed_rare_minus_common_source_gap=signed,
        k_adjusted_rare_minus_common_source_gap=k_signed,
        original_genus_jaccard_turnover=1-len(gensA&gensB)/len(gensA|gensB),
        original_colxr_candidate_family_jaccard_turnover=(
           1-len(famA&famB)/len(famA|famB) if famA|famB else np.nan),
        n_original_shared_genera=len(gensA&gensB),
        n_original_shared_families=len(famA&famB),
    )


def _weighted_slope(rows:pd.DataFrame,xcol:str,ycol:str,cluster:str,
                    *,adjusted:bool)->float:
    if len(rows)<20 or rows[cluster].nunique()<10:
        raise ValueError("Insufficient independent original clusters for regression")
    d=rows.copy()
    x=d[xcol].to_numpy(float)
    y=d[ycol].to_numpy(float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Nonfinite original phylogeny/functional response")
    sd=float(np.std(x))
    if sd<1e-12:
        raise ValueError("No physical-island variation in phylogenetic beta")
    z=(x-np.mean(x))/sd
    predictors=[z]
    if adjusted:
        for key in (
          "original_genus_jaccard_turnover",
          "original_colxr_candidate_family_jaccard_turnover",
          "log_source_exact_tip_candidates",
          "log_island_area_km2",
          "log_distance_to_continent_km",
        ):
            values=d[key].to_numpy(float)
            ss=np.std(values)
            if ss<1e-12:
                raise ValueError("No original regional confounder variation")
            predictors.append((values-values.mean())/ss)
    X=np.column_stack([np.ones(len(d)),*predictors])
    # Strict original source regional/regression covariance rank gate. Do not
    # silently drop a confounder to rescue a directional finding.
    if np.linalg.matrix_rank(X)<X.shape[1]:
        raise ValueError("Original phylogenetic/richness confounders collinear")
    weights=1/d.groupby(cluster)[cluster].transform("size").to_numpy(float)
    weights=weights/weights.sum()
    root=np.sqrt(weights)
    return float(np.linalg.lstsq(X*root[:,None],y*root,rcond=None)[0][1])


def _boot_beta(rows:pd.DataFrame,ycol:str,*,cluster:str,adjusted:bool,
               draws:int,seed:int):
    if rows[cluster].nunique()<10:
        return None
    try:
        observed=_weighted_slope(rows,"phylosor_beta",ycol,cluster,adjusted=adjusted)
    except ValueError:
        return None
    groups={str(k):v for k,v in rows.groupby(cluster)}
    keys=sorted(groups)
    rng=np.random.default_rng(seed)
    effects=[]
    for _ in range(draws):
        sampled=rng.choice(keys,size=len(keys),replace=True)
        # resampled original independent clusters: maintain distinct
        # bootstrap cluster IDs, not recompute weights against duplicate
        # exact strings after a given source was sampled twice.
        blocks=[]
        for j,k in enumerate(sampled):
            g=groups[str(k)].copy()
            g["_bootstrap_cluster_id"]=j
            blocks.append(g)
        data=pd.concat(blocks,ignore_index=True)
        try:
            effects.append(_weighted_slope(
                data,"phylosor_beta",ycol,"_bootstrap_cluster_id",
                adjusted=adjusted))
        except ValueError:
            continue
    if len(effects)<0.8*draws:
        return dict(status="NONESTIMABLE_BOOTSTRAP_PREDICTOR_RANK",
                    n_valid=len(effects),
                    observed_std_phylo_beta_slope=observed)
    lo,hi=np.quantile(effects,[.025,.975])
    return dict(status="EXPLORATORY_CLUSTER_CONDITIONAL_ASSOCIATION",
                observed_std_phylo_beta_slope=observed,
                cluster_025=float(lo),cluster_975=float(hi),
                n_valid_block_resamples=len(effects))


def macro_original_phylo_functional_beta(
    triple:pd.DataFrame,floras:pd.DataFrame,geo:pd.DataFrame,
    original_source_assignment:pd.DataFrame,
    original_GIFT_native_mainland:pd.DataFrame,
    tree:OriginalTreeSubtree,exact_original_tip_names:list[str],
    exact_tip_node_map:dict[str,int], original_family:dict[str,str],
    *,tree_version:str,repeats:int=20,reference_null_draws:int=199,
    bootstrap:int=999,seed:int=20261114,
)->dict:
    if tree_version not in ("TPL","LCVP","WP"):
        raise ValueError("No original dated GBOTB source version")
    if repeats<20 or reference_null_draws<199 or bootstrap<999:
        raise ValueError("Original phylo/functional beta protocol support weakened")
    if len(triple)!=2234 or len(exact_original_tip_names)<1100:
        raise ValueError("Pinned original Direct/source phylogenetic species cohort changed")
    if len(set(exact_original_tip_names))!=len(exact_original_tip_names):
        raise ValueError("Source phylogenetic exact tip names duplicated")
    original,places,gift_src,gift,universe=_strict_source_taxon_groups(
        triple,floras,geo,original_source_assignment,original_GIFT_native_mainland)
    if len(universe)!=928:
        raise ValueError("Original 928 island Direct macro population changed")
    if not {"island_id","log_island_area_km2",
            "log_distance_to_continent_km"}.issubset(geo):
        raise ValueError("Original island area and geographical isolation needed")
    geographic=geo.set_index("island_id")
    ref=RarityReference.from_original(triple,null_draws=reference_null_draws,seed=seed)
    phen={
        str(name):tuple(row) for name,row in zip(
            triple.index,triple[list(TRAITS)].to_numpy(dtype=str))}
    true=set(exact_original_tip_names)
    if not true.issubset(phen) or not true.issubset(original_family):
        raise ValueError("Only original source verified family/trait true tips may enter")
    if true!=set(exact_tip_node_map) or (
        set(exact_tip_node_map.values())!=set(tree.original_tip_paths)):
        raise ValueError("Dated source original tree tip-node mapping not strictly equivalent")
    per_island={i:v for i,v in original.groupby("island_id")}
    store=[]
    for island in sorted(universe):
        original_flora=per_island[island]
        full=set(original_flora.accepted_species)
        pool=(set().union(*(gift.get(x,set()) for x in gift_src[island]))
              & true)-full
        if len(pool)<N_RAREFY:
            continue
        for scope,names in (
            ("all_recorded",full),
            ("source_native",set(original_flora.loc[
                original_flora.origin_status.eq("native"),"accepted_species"]))):
            if len(names)<30:
                continue
            local=sorted(names&true)
            donor=sorted(pool)
            if len(local)<N_RAREFY:
                continue
            random=np.random.default_rng(_seedstr(seed,tree_version,island,scope))
            reps=[]
            for _ in range(repeats):
                A=random.choice(local,size=N_RAREFY,replace=False).tolist()
                B=random.choice(donor,size=N_RAREFY,replace=False).tolist()
                reps.append(_phylo_functional_draw(
                    ref,tree,A,B,exact_tip_node_map,phen,original_family))
            summarized=pd.DataFrame(reps).mean(numeric_only=True).to_dict()
            store.append(dict(
                island_id=island,scope=scope,
                regime=str(places.at[island,"analysis_regime"]),
                spatial_block=str(places.at[island,"spatial_block"]),
                source_signature=gift_src[island][0],
                log_source_exact_tip_candidates=float(np.log1p(len(pool))),
                log_island_area_km2=float(
                    geographic.at[island,"log_island_area_km2"]),
                log_distance_to_continent_km=float(
                    geographic.at[island,"log_distance_to_continent_km"]),
                **summarized))
    df=pd.DataFrame(store)
    reports=[]
    responses=("functional_support_beta",
       "exact_profile_jaccard_beta",
       "k_adjusted_rare_minus_common_source_gap")
    for (scope,regime),g in df.groupby(["scope","regime"]):
        n=len(g)
        blocks=g.spatial_block.nunique()
        sources=g.source_signature.nunique()
        pass_gate=n>=REGION_MIN and blocks>=BLOCK_MIN and sources>=SOURCE_MIN
        report=dict(scope=scope,regime=regime,
            n_original_eligible_islands=n,n_physical_spatial_blocks=blocks,
            n_first_nearest_GIFT_source_signatures=sources,
            n_true_tip_focal_species_per_side=N_RAREFY,
            original_flora_repeats_per_island=repeats,
            passed_20islands_10blocks_10source_support=bool(pass_gate),
            median_original_phylo_beta=float(g.phylosor_beta.median()),
            median_original_functional_beta=float(g.functional_support_beta.median()),
            median_original_genus_turnover=float(
                g.original_genus_jaccard_turnover.median()),
            median_original_family_turnover=float(
                g.original_colxr_candidate_family_jaccard_turnover.median()),
        )
        if pass_gate:
            for i,y in enumerate(responses):
                for adj in (False,True):
                    for j,cluster in enumerate(("spatial_block","source_signature")):
                        result=_boot_beta(g,ycol=y,cluster=cluster,
                           adjusted=adj,draws=bootstrap,seed=seed+97*i+31*int(adj)+j)
                        report[f"{y}_{'adjusted' if adj else 'raw'}_"+
                            ("block" if j==0 else "source")]=result
        reports.append(report)
    return dict(
        status="ORIGINAL_TRUE_DATED_PHYLO_BETA_VS_CATEGORICAL_FUNCTIONAL_BETA",
        n_original_928_direct_island_population=928,
        n_original_global_direct_flower_species=2234,
        n_original_true_phylogenetic_tip_species=len(true),
        original_tree_backbone_version=tree_version,
        region_support_and_categorical_beta_relationships=reports,
        n_original_private_island_scope_matched_rows=len(df),
        original_design={
            "phylogenetic_beta":"one_minus_rooted_original_branch_length_PhyloSor",
            "functional_beta":"one_minus_Jaccard_original_reference_one_trait_neighborhood",
            "secondary":"signed_rare_minus_common_original_k_adjusted_coverage",
            "same_15_true_dated_plant_tips_on_island_and_GIFT_source":True,
            "n_per_island_original_source_random_rarefactions":repeats,
            "n_reference_observed_global_flower_profiles":100,
            "no_taxonomically_imputed_phylogenetic_plant_tips":True,
            "remove_all_same_island_species_from_source_GIFT_floras":True,
            "family_genus_source_area_isolation_adjusted_regressions":True,
            "equal_physical_block_and_GIFT_source_bootstraps":bootstrap,
            "one_original_physical_island_one_regression_row":True,
            "conditional_estimates_not_source_independent_ancestries":True,
            "original_plant_trait_assignments_are_species_fixed_not_evolving":True,
            "phylogenetic_source_distance_does_not_prove_ancestor":True,
        },
        no_actual_pollinator_guild_decline_or_evolutionary_trait_change_measured=True,
        original_source_tip_genus_family_island_or_coordinate_ids_not_exported=True,
    )
