"""Source-pinned original EMLN island flower-visitor network contrasts.

Scientific estimands:
EMLN60 Canary 5 geographic islands, six paired island locations versus the
single Western Sahara mainland location (2 replicate plots each): same plant
species' documented VISITOR-SPECIES turnover. Match plant identity and equal
one-site versus one-site pairings; also require the plant to occur in ALL
four site lists for the strict site-replicate negative control.
EMLN71 Mauritius ONE invaded vs ONE actively managed site across 12
consecutive 2-week network periods: native plant taxon richness, documented
visitor portfolio and matched same-plant rewiring.
No historical pollinator decline, bird/bee replacement, effective pollen
delivery, seed set, plant morphological trait changes or causal site effect.
"""
from __future__ import annotations
from collections import Counter
import numpy as np
import pandas as pd

NODES=("node_id","attribute","value")
EDGES=("interaction_id","attribute","value")
LAYERS=("layer","attribute","value")
CANARY_SITE_NAMES={
    "El_Hierro","Fuerteventura","Gran_Canaria","La_Gomera",
    "Tenerife_South","Tenerife_Teno","Western_Sahara"}
EXPECTED_SITE_LAYER_COUNT=14
EXPECTED_MAURITIUS_LAYER_COUNT=24


def _pivot(raw:pd.DataFrame,keys:tuple[str,...])->pd.DataFrame:
    if not set(keys+("attribute","value")).issubset(raw.columns):
        raise ValueError("Original EMLN format missing source keys or attributes")
    if raw[list(keys)].isna().any().any():
        raise ValueError("Original EMLN layer/interaction/taxon IDs missing")
    base=raw[list(keys)+["attribute","value"]].copy().fillna("")
    # Original EMLN71 has two legitimate layer.type labels: time and
    # perturbation. All other key×attribute duplicates are forbidden.
    duplicates=base.duplicated(list(keys)+["attribute"],keep=False)
    if duplicates.any():
        d=base.loc[duplicates]
        if not (keys==("layer",) and d.attribute.eq("type").all()):
            raise ValueError("Conflicting original EMLN node/edge/source fields")
    dedup=base.drop_duplicates(list(keys)+["attribute","value"])
    wide=(dedup.groupby(list(keys)+["attribute"],sort=False).value
          .agg(lambda v:"|".join(dict.fromkeys(map(str,v))))
          .unstack("attribute").reset_index().fillna(""))
    return wide


def original_emln_source(
    node_long:pd.DataFrame,edge_long:pd.DataFrame,layer_long:pd.DataFrame,
    *, expected:int,
)->tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    nodes=_pivot(node_long,("node_id",))
    edges=_pivot(edge_long,("interaction_id",))
    layers=_pivot(layer_long,("layer",))
    if len(layers)!=expected:
        raise ValueError("Original EMLN source-layer count changed")
    if not {"node_name","type"}.issubset(nodes):
        raise ValueError("Original EMLN plant/pollinator node roles absent")
    if not {"node_from","node_to","layer_from","layer_to","weight"}.issubset(edges):
        raise ValueError("Original EMLN visitor link format changed")
    if not {"name"}.issubset(layers):
        raise ValueError("Original EMLN network names absent")
    if nodes.node_name.duplicated().any() or layers.name.duplicated().any():
        raise ValueError("EMLN node or layer names not unique")
    t=nodes.set_index("node_name").type
    a=edges.node_from.map(t)
    b=edges.node_to.map(t)
    if not ((a.eq("plant")&b.eq("pollinator")) |
            (a.eq("pollinator")&b.eq("plant"))).all():
        raise ValueError("Source interaction endpoints not plant-pollinator")
    if not (edges.layer_from.eq(edges.layer_to)).all():
        raise ValueError("Cross-layer visitor links do not represent local visits")
    edges["plant"]=np.where(a.eq("plant"),edges.node_from,edges.node_to)
    edges["visitor"]=np.where(a.eq("pollinator"),edges.node_from,edges.node_to)
    edges["weight"]=pd.to_numeric(edges.weight,errors="coerce")
    if edges.weight.isna().any() or edges.weight.le(0).any():
        raise ValueError("Original flower visit weights must be positive")
    return nodes,edges,layers


def _jaccard(a:set[str],b:set[str])->float:
    return len(a & b)/len(a | b) if (a|b) else 0.


def canary_same_plant_partner_turnover(
    nodes:pd.DataFrame, edges:pd.DataFrame, layers:pd.DataFrame,
)->dict:
    """Equal individual-site comparison and STRICT 4-site same-plant control."""
    if len(layers)!=EXPECTED_SITE_LAYER_COUNT or len(edges)!=969:
        raise ValueError("Original pinned Canary network source support drifted")
    if "latitude" not in layers or "longitude" not in layers:
        raise ValueError("Original canary study site geography missing")
    if not all(layers.name.str.contains(r"_site_[12]$",regex=True)):
        raise ValueError("Canary two-site source contract changed")
    ll=layers.copy()
    ll["place"]=ll.name.str.replace(r"_site_[12]$","",regex=True)
    if set(ll.place)!=CANARY_SITE_NAMES or not ll.groupby("place").size().eq(2).all():
        raise ValueError("Original 7 locations each with two sites changed")
    v=edges.merge(ll[["name","place"]],left_on="layer_from",
                  right_on="name",validate="many_to_one")
    per_site={
        site:sub.groupby("plant").visitor.agg(set).to_dict()
        for site,sub in v.groupby("name")
    }
    site_names={place:sorted(group.name) for place,group in ll.groupby("place")}
    mainland=site_names["Western_Sahara"]
    records=[]
    for place in sorted(CANARY_SITE_NAMES-{"Western_Sahara"}):
        ia,ib=site_names[place]
        ma,mb=mainland
        island_common=set(per_site[ia])&set(per_site[ib])
        site_within=np.array([
            _jaccard(per_site[ia][p],per_site[ib][p])
            for p in island_common],float)
        between=[]
        between_plants=set()
        for i in (ia,ib):
            for m in (ma,mb):
                common=set(per_site[i])&set(per_site[m])
                between_plants|=common
                between.append(
                    float(np.mean([_jaccard(per_site[i][plant],
                                            per_site[m][plant])
                                 for plant in common])) if common else np.nan)
        strict=(set(per_site[ia])&set(per_site[ib])&
                set(per_site[ma])&set(per_site[mb]))
        within_strict=np.array([
            _jaccard(per_site[ia][p],per_site[ib][p])
            for p in sorted(strict)],float)
        mainland_strict=np.array([
            _jaccard(per_site[ma][p],per_site[mb][p])
            for p in sorted(strict)],float)
        cross_strict=np.array([
            np.mean([_jaccard(per_site[i][p],per_site[m][p])
                for i in (ia,ib) for m in (ma,mb)])
            for p in sorted(strict)],float)
        records.append(dict(
            place=place,
            geographical_island="Tenerife" if place.startswith("Tenerife") else place,
            n_shared_plants_between_two_island_sites=len(island_common),
            n_distinct_same_plants_against_mainland_sites=len(between_plants),
            n_strict_same_plants_present_all_four_sites=len(strict),
            mean_within_island_site_partner_jaccard=float(site_within.mean()),
            mean_between_single_sites_partner_jaccard=float(np.nanmean(between)),
            strictly_matched_same_plant_within_site_jaccard=float(within_strict.mean()),
            strictly_matched_same_plant_mainland_site_jaccard=float(mainland_strict.mean()),
            strictly_matched_same_plant_island_mainland_jaccard=float(cross_strict.mean()),
            strict_within_minus_mainland_partner_overlap=(
                float((within_strict-cross_strict).mean())),
        ))
    x=pd.DataFrame(records)
    by_island=x.groupby("geographical_island")[[
        "strict_within_minus_mainland_partner_overlap",
        "strictly_matched_same_plant_within_site_jaccard",
        "strictly_matched_same_plant_island_mainland_jaccard",
    ]].mean()
    # Geographically nearest alternative island: plant names must occur in
    # all THREE source locations. Control island is picked on distance alone.
    # This is a DESCRIPTIVE 6-location sensitivity, with only 2–5 plants per
    # three-way match and one continental locality; not an evolutionary
    # phylogenetic or independent source-region contrast.
    centre=ll.groupby("place")[["latitude","longitude"]].first().astype(float)
    aggregate={
        place:g.groupby("plant").visitor.agg(set).to_dict()
        for place,g in v.groupby("place")
    }
    def greatcircle(a,b):
        ax=np.deg2rad(centre.loc[a,["latitude","longitude"]].to_numpy(float))
        bx=np.deg2rad(centre.loc[b,["latitude","longitude"]].to_numpy(float))
        h=(np.sin((ax[0]-bx[0])/2.)**2+
           np.cos(ax[0])*np.cos(bx[0])*np.sin((ax[1]-bx[1])/2.)**2)
        return float(6371.0088*2*np.arcsin(np.sqrt(np.clip(h,0,1))))
    geo_sensitivity=[]
    for place in sorted(CANARY_SITE_NAMES-{"Western_Sahara"}):
        distance_mainland=greatcircle(place,"Western_Sahara")
        alternatives=sorted(
            (name for name in CANARY_SITE_NAMES
             if name not in (place,"Western_Sahara")
             and not (place.startswith("Tenerife_") and
                      name.startswith("Tenerife_"))))
        alternate=min(alternatives,key=lambda name:
                      abs(greatcircle(place,name)-distance_mainland))
        common=(set(aggregate[place])&
                set(aggregate["Western_Sahara"])&
                set(aggregate[alternate]))
        jmain=np.mean([_jaccard(aggregate[place][plant],
                     aggregate["Western_Sahara"][plant]) for plant in common])
        jalt=np.mean([_jaccard(aggregate[place][plant],
                     aggregate[alternate][plant]) for plant in common])
        geo_sensitivity.append({
           "place":place,"distance_matched_other_island":alternate,
           "km_to_original_mainland":distance_mainland,
           "km_to_distance_matched_other_island":greatcircle(place,alternate),
           "n_same_plant_species_at_all_three_locations":len(common),
           "same_plant_mean_visitors_jaccard_with_mainland":float(jmain),
           "same_plant_mean_visitors_jaccard_with_other_island":float(jalt),
           "other_island_partner_overlap_exceeds_mainland":bool(jalt>jmain),
        })
    return {
        "status":"ORIGINAL_CANARY_SAME_PLANT_PARTNER_TURNOVER_ONLY",
        "source_reference":"Troejelsgaard_et_al_2015_EMLN60",
        "n_original_visitor_links":len(edges),
        "n_source_plant_species":int(nodes.type.eq("plant").sum()),
        "n_source_visitor_taxa":int(nodes.type.eq("pollinator").sum()),
        "n_geographical_islands":int(len(by_island)),
        "n_island_study_locations":len(x),
        "n_independent_mainland_study_locations":1,
        "n_total_replicate_source_sites":len(layers),
        "n_island_locations_where_strict_site_matched_within_exceeds_cross":int(
            x.strict_within_minus_mainland_partner_overlap.gt(0).sum()),
        "n_geographical_islands_where_strict_within_exceeds_cross":int(
            by_island.strict_within_minus_mainland_partner_overlap.gt(0).sum()),
        "mean_location_equal_within_island_site_partner_jaccard":float(
            x.mean_within_island_site_partner_jaccard.mean()),
        "mean_location_equal_between_island_mainland_single_site_partner_jaccard":float(
            x.mean_between_single_sites_partner_jaccard.mean()),
        "strict_match_location_equal_within_jaccard":float(
            x.strictly_matched_same_plant_within_site_jaccard.mean()),
        "strict_match_location_equal_cross_jaccard":float(
            x.strictly_matched_same_plant_island_mainland_jaccard.mean()),
        "strict_match_geographical_island_equal_within_jaccard":float(
            by_island.strictly_matched_same_plant_within_site_jaccard.mean()),
        "strict_match_geographical_island_equal_cross_jaccard":float(
            by_island.strictly_matched_same_plant_island_mainland_jaccard.mean()),
        "n_strict_plant_location_combinations":int(
            x.n_strict_same_plants_present_all_four_sites.sum()),
        "n_distance_matched_three_source_locations_with_same_plant_species":int(
            sum(v["n_same_plant_species_at_all_three_locations"]>=2
                for v in geo_sensitivity)),
        "n_distance_matched_island_location_triplets_where_mainland_visitor_overlap_lower":int(
            sum(v["other_island_partner_overlap_exceeds_mainland"]
                for v in geo_sensitivity)),
        "distance_matched_island_island_vs_mainland_source_sensitivity":
            geo_sensitivity,
        "locations":[{
            k:(int(row[k]) if k.startswith("n_") else
               float(row[k]) if k not in ("place","geographical_island") else row[k])
            for k in row.index}
            for _,row in x.iterrows()],
        "source_sampling_one_site_pair_each_not_independent_taxon_histories":True,
        "single_shared_mainland_and_taxonomic_regional_pool_dependence":True,
        "historical_pollinator_decline_adaptation_or_effectiveness_measured":False,
        "raw_original_plant_visitor_and_site_ids_exported":False,
    }


def mauritius_native_partner_rewiring(
    nodes:pd.DataFrame, edges:pd.DataFrame, layers:pd.DataFrame,
)->dict:
    """ONE managed/restored vs ONE invaded/unmanaged site over 12 periods."""
    if len(layers)!=EXPECTED_MAURITIUS_LAYER_COUNT or len(edges)!=1953:
        raise ValueError("Original Mauritius source panel drifted")
    if not {"site","period","perturbed"}.issubset(layers):
        raise ValueError("Original temporal/perturbation descriptors missing")
    if not {"origin","order"}.issubset(nodes):
        raise ValueError("Original pollinator order and plant native status absent")
    l=layers.copy()
    l["period"]=pd.to_numeric(l.period,errors="raise").astype(int)
    if (set(l.period)!=set(range(1,13)) or
        set(l.site)!={"CMA","Control"} or
        len(l.groupby(["site","period"]))!=24 or
        l.duplicated(["site","period"]).any() or
        not l.loc[l.site.eq("CMA"),"perturbed"].eq("FALSE").all() or
        not l.loc[l.site.eq("Control"),"perturbed"].eq("TRUE").all()):
        raise ValueError("Original one-managed one-control 12-period design drift")
    orig=nodes.set_index("node_name").origin
    order=nodes.set_index("node_name").order
    e=edges.merge(l[["name","site","period"]],left_on="layer_from",
        right_on="name",validate="many_to_one")
    e["plant_origin"]=e.plant.map(orig)
    e["visitor_order"]=e.visitor.map(order)
    native=e.loc[e.plant_origin.isin(("Native","Endemic"))].copy()
    if native.empty or native.visitor_order.isna().any():
        raise ValueError("Original native flower visitor source unsupported")
    rows=[]
    for period in range(1,13):
        per=native.loc[native.period.eq(period)]
        a=per.loc[per.site.eq("CMA")]
        b=per.loc[per.site.eq("Control")]
        def metrics(g):
            return dict(
              n_native_visited_plants=int(g.plant.nunique()),
              n_native_visitor_taxa=int(g.visitor.nunique()),
              n_native_visitor_links=int(len(g)),
              mean_native_visitors_per_visited_plant=float(
                  g.groupby("plant").visitor.nunique().mean()),
              frac_native_plants_with_2plus_visitor_orders=float(
                  g.groupby("plant").visitor_order.nunique().ge(2).mean()))
        m1=metrics(a);m2=metrics(b)
        same=set(a.plant)&set(b.plant)
        va=a.groupby("plant").visitor.agg(set).to_dict()
        vb=b.groupby("plant").visitor.agg(set).to_dict()
        oa=a.groupby("plant").visitor_order.agg(set).to_dict()
        ob=b.groupby("plant").visitor_order.agg(set).to_dict()
        common_rows=[dict(
           deg_dif=len(va[p])-len(vb[p]),
           order_dif=len(oa[p])-len(ob[p]),
           documented_partner_species_jaccard=_jaccard(va[p],vb[p]),
           documented_partner_order_jaccard=_jaccard(oa[p],ob[p]))
           for p in same]
        z=pd.DataFrame(common_rows)
        if len(z)<3:
            raise ValueError("Insufficient original same-native-plant comparison")
        rows.append(dict(period=period,
          n_same_native_plant_species=len(same),
          n_native_visited_plant_difference_CMA_minus_Control=(
            m1["n_native_visited_plants"]-m2["n_native_visited_plants"]),
          n_native_visitor_taxon_difference_CMA_minus_Control=(
            m1["n_native_visitor_taxa"]-m2["n_native_visitor_taxa"]),
          mean_native_visitor_per_plant_difference_CMA_minus_Control=(
            m1["mean_native_visitors_per_visited_plant"]-
            m2["mean_native_visitors_per_visited_plant"]),
          frac_native_plants_multiple_orders_difference_CMA_minus_Control=(
            m1["frac_native_plants_with_2plus_visitor_orders"]-
            m2["frac_native_plants_with_2plus_visitor_orders"]),
          same_plant_mean_visitor_degree_difference_CMA_minus_Control=float(z.deg_dif.mean()),
          same_plant_mean_order_degree_difference_CMA_minus_Control=float(z.order_dif.mean()),
          same_plant_median_partner_species_jaccard=float(
              z.documented_partner_species_jaccard.median()),
          same_plant_median_partner_order_jaccard=float(
              z.documented_partner_order_jaccard.median()),
        ))
    res=pd.DataFrame(rows)
    metrics=(
       "n_native_visited_plant_difference_CMA_minus_Control",
       "n_native_visitor_taxon_difference_CMA_minus_Control",
       "mean_native_visitor_per_plant_difference_CMA_minus_Control",
       "frac_native_plants_multiple_orders_difference_CMA_minus_Control",
       "same_plant_mean_visitor_degree_difference_CMA_minus_Control",
       "same_plant_mean_order_degree_difference_CMA_minus_Control",
       "same_plant_median_partner_species_jaccard",
       "same_plant_median_partner_order_jaccard"
    )
    return {
       "status":"ORIGINAL_MAURITIUS_SINGLE_PAIRED_SITE_TEMPORAL_VISITOR_TURNOVER",
       "source_reference":"Kaiser_Bunbury_et_al_2009_EMLN71",
       "n_original_visitor_links":int(len(edges)),
       "n_native_or_endemic_plant_species_in_original_nodes":int(
          nodes.loc[nodes.type.eq("plant"),"origin"].isin(("Native","Endemic")).sum()),
       "n_two_week_periods_with_two_survey_sites":12,
       "n_independent_managed_sites":1,
       "n_independent_unmanaged_invaded_sites":1,
       "n_periods_more_native_flower_species_in_managed_site":int(
           res.n_native_visited_plant_difference_CMA_minus_Control.gt(0).sum()),
       "n_periods_more_native_visitor_taxa_in_managed_site":int(
           res.n_native_visitor_taxon_difference_CMA_minus_Control.gt(0).sum()),
       "period_equal_mean_original_metrics":{
           v:float(res[v].mean()) for v in metrics},
       "periods":[{
           k:int(row[k]) if k in ("period","n_same_native_plant_species") else float(row[k])
           for k in row.index
       } for _,row in res.iterrows()],
       "repeated_periods_not_independent_experimental_replicates":True,
       "original_site_invasion_or_restoration_not_randomized":True,
       "partner_guild_record_not_pollen_delivery_effectiveness":True,
       "raw_original_plant_visitor_and_site_ids_exported":False,
    }
