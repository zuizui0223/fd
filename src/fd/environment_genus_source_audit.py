"""Source-only audit of original 4505 GIFT island environmental genus profiles.

Input must be the protected original anonymous genotype/genus forecast result.
This is a checklist GENUS COMPOSITION scoring audit, not a causal genus
establishment or observed mainland / pollinator-guild transition test.
Both multiclass Brier and proper log-loss, and BOTH original independent
physical-block/first source-entity group equal bootstrap intervals,
must support any directional source predictor versus regional baseline.
"""
from __future__ import annotations

import math

REGIONS=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude")
SCOPES=("all_recorded","original_native")
MODELS=("geo_nearest","climate_nearest","environment_nearest",
        "environment_source_nearest")
GROUPS=("physical_10degree_block","first_candidate_source_entity_group")
QUALIFIED="ORIGINAL_GIFT_TAXONOMIC_GENUS_SOURCE_COMPOSITION_GROUP_PREDICTION"
FULL_STATUS="ORIGINAL_GIFT_GENUS_OCCUPANCY_COMPOSITION_ENVIRONMENTAL_TEST_BOTH_SOURCE_GROUPS"


def _finite(z):
    return isinstance(z,(int,float)) and type(z)!=bool and math.isfinite(z)


def _group(r:dict,key:str):
    if r.get("status")!=QUALIFIED or r.get("n_source_groups",0)<10:
        return None
    val=r.get(key+"_group_equal_mean")
    interval=r.get(key+"_group_resample_95")
    if (not _finite(val) or not isinstance(interval,list)
        or len(interval)!=2 or not all(_finite(t) for t in interval)
        or interval[0]>interval[1]):
        return None
    return dict(value=float(val),ci_95=[float(x) for x in interval],
                n_source_groups=int(r["n_source_groups"]),
                n_scored_floras=int(r["n_original_source_scored_floras"]))


def _direction(source_row:dict,model:str):
    comparisons={}
    for field in ("Brier","logloss"):
        for g in GROUPS:
            val=source_row.get(g)
            comparisons[field+"__"+g]=(_group(val,
                model+"_minus_region_equal_"+field) if isinstance(val,dict)
                else None)
    if not all(comparisons.values()):
        decision="not_identifiable"
    elif all(x["ci_95"][0]>0 for x in comparisons.values()):
        decision="predictor_better_than_region_equal"
    elif all(x["ci_95"][1]<0 for x in comparisons.values()):
        decision="predictor_worse_than_region_equal"
    else:
        decision="ambiguous_or_disagreeing_losses_or_source_groups"
    return dict(status=decision,source_grouped_proper_loss_evidence=comparisons)


def audit_geo_source_genus_prediction(original:dict)->dict:
    if (original.get("status")!=
        "ORIGINAL_GIFT_4505_NO_FOCAL_GENUS_SOURCE_ISLAND_ENVIRONMENTAL_GENUS_ASSEMBLY"
        or original.get("original_GIFT_parent_island_floras")!=4505
        or original.get("source_genus_models")!=list(MODELS)):
        raise ValueError("Original 4505 global GIFT genus source analysis not authenticated")
    source=original.get("controls",{})
    required=(
        "original_4505_island_botanical_flora_population_retained",
        "original_focal_true_genus_community_not_used_to_predict",
        "original_true_accepted_species_genus_used_ONLY_FOR_SCORING",
        "original_training_GIFT_islands_in_focal_10degree_block_EXCLUDED",
        "original_training_GIFT_islands_with_same_first_GIFT_candidate_mainland_source_EXCLUDED",
        "GIFT_source_donor_island_taxonomy_NOT_independent_mainland_ancestor",
        "source_species_normalized_by_island_and_UNSEEN_GENUS_held_in_common_unknown_category",
        "source_taxon_genus_nonlisting_is_NOT_confirmed_ecological_absence",
        "no_true_island_genus_colonization_probability_OR_independent_trait_adaptation",
        "no_pollinator_guild_selection_or_seedset_causality_claim",
        "source_group_bootstrap_exploratory_not_multiplicity_adjusted",
    )
    for label in required:
        if source.get(label) is not True:
            raise ValueError("Original source scientific provenance ceiling missing: "+label)
    cells=original.get("original_no_oracle_original_island_genus_composition_source_predictions")
    if not isinstance(cells,list) or len(cells)!=8:
        raise ValueError("Original two nested status × four broad regions incomplete")
    seen=set()
    result=[]
    for row in cells:
        key=(row.get("scope"),row.get("region"))
        if key in seen or key[0] not in SCOPES or key[1] not in REGIONS:
            raise ValueError("Unexpected original GIFT source geography/status panel")
        seen.add(key)
        parent=row.get("n_original_region_botanical_floras")
        eligible=row.get("n_original_source_islands_with_full_geo_and_genus_baseline",0)
        scored=row.get("n_original_floras_scored_for_genus_community_environmental_prediction",0)
        missing=row.get("n_original_focal_islands_no_complete_geographic_or_full_genus_record",0)
        dropout=row.get("n_unscored_donor_coverage_by_reason",{})
        if (not all(type(x)==int and x>=0 for x in (parent,eligible,scored,missing))
                or scored>eligible or eligible+missing!=parent
                or not isinstance(dropout,dict) or
                sum(dropout.values())+scored!=eligible):
            raise ValueError("Original 4505 source genus source/geo eligibility waterfall wrong")
        outcome={m:_direction(row,m) for m in MODELS}
        if row.get("status")==FULL_STATUS and (
            not all(outcome[m]["status"]!="not_identifiable" for m in MODELS)):
            raise ValueError("Fully spatially qualified original output lacks proper score metrics")
        result.append(dict(
            scope=key[0],region=key[1],n_original_region_islands=parent,
            n_original_geo_and_flora_eligible_islands=eligible,
            n_original_islands_scored_without_local_genus_oracle=scored,
            n_missing_source_geography_or_minimum_checklist=missing,
            n_source_donor_insufficient_islands=sum(dropout.values()),
            source_nonidentifiability_counts=dropout,
            original_source_models=outcome,
            first_source_candidate_may_be_geography_proxy_not_measured_immigration=True,
        ))
    if seen!={(scope,r) for scope in SCOPES for r in REGIONS}:
        raise ValueError("Original 4-regime global flora study missing status panels")
    summary=[]
    for scope in SCOPES:
        for model in MODELS:
            cells_scope=[r for r in result if r["scope"]==scope]
            summary.append(dict(
                botanical_scope=scope,original_genus_predictor=model,
                n_independently_qualified_macroregions=sum(
                    r["original_source_models"][model]["status"]!=
                    "not_identifiable" for r in cells_scope),
                n_regions_outperform_equal_baseline_on_both_proper_losses_and_groups=sum(
                    r["original_source_models"][model]["status"]==
                    "predictor_better_than_region_equal" for r in cells_scope),
                n_regions_below_equal_baseline_on_both_proper_losses_and_groups=sum(
                    r["original_source_models"][model]["status"]==
                    "predictor_worse_than_region_equal" for r in cells_scope),
            ))
    return dict(
        status="ORIGINAL_GIFT_4505_GEOGRAPHIC_AND_ENVIRONMENTAL_GENUS_ASSEMBLY_SOURCE_AUDIT",
        n_original_global_botanical_floras=4505,
        original_full_geographic_composition_cells=result,
        original_macroregion_proper_loss_source_count_summary=summary,
        limitations=dict(
            original_focal_botanical_genus_input_never_used=True,
            genus_community_richness_does_not_measure_genus_colonization_probability=True,
            whole_species_checklist_not_seen_is_NOT_biological_genus_absence=True,
            original_source_islands_not_independent_mainland_ancestors=True,
            original_candidate_mainland_entity_similarity_NOT_dispersion_or_selection=True,
            paired_environment_and_environment_plus_source_effect_NOT_tested_for_direct_increment=True,
            original_genus_community_source_not_yet_flora_trait_prediction=True,
            original_direct_flower_selection_pollination_service_NOT_estimated=True,
            sample_composition_annotation_not_phylogenetic_independent_replicates=True,
            group_equal_source_intervals_exploratory_not_multiplicity_corrected=True,
        ),
    )
