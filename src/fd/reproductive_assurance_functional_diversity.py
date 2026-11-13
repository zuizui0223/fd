"""Independent reproductive-assurance × joint-flower-profile support/assembly audit.

Original species Direct high/medium only. No low-confidence genus imputation,
no inference that SC is autonomous selfing or that reproductive states confer
effective pollen service. Island and mainland source records remain private.
"""
from __future__ import annotations
from collections import Counter
import hashlib
import numpy as np
import pandas as pd

from .pollinator_functional_diversity import direct_triple

REPRO = {
 "self_compatibility":("self_incompatibility",{"SC"},{"SI"}),
 "autonomous_selfing":("autonomous_selfing_capacity",{"autonomous","delayed"},{"absent"}),
 "selfing_mating_system":("mating_system",{"predominantly_selfing","obligate_selfing"},{"predominantly_outcrossing"}),
}


def direct_binary_assurance(ledger:pd.DataFrame) -> pd.DataFrame:
    req={"accepted_species","trait_name","normalized_value",
         "evidence_scope","evidence_quality"}
    if not req.issubset(ledger):
        raise ValueError("Direct source evidence columns missing")
    x=ledger.loc[
        ledger.trait_name.isin([v[0] for v in REPRO.values()]) &
        ledger.evidence_scope.eq("species_direct") &
        ledger.evidence_quality.isin(("high","medium")) &
        ledger.normalized_value.notna() &
        ~ledger.normalized_value.astype(str).str.contains("|",regex=False),
        list(req)].copy()
    if x.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Conflicting Direct reproductive trait")
    scored={}
    for trait,(original,pos,neg) in REPRO.items():
        sub=x.loc[x.trait_name.eq(original),["accepted_species","normalized_value"]].copy()
        sub=sub.loc[sub.normalized_value.isin(pos|neg)]
        if sub.accepted_species.duplicated().any():
            raise ValueError("Duplicate Direct trait after recoding")
        scored[trait]=sub.assign(state=sub.normalized_value.isin(pos).astype(int)).set_index(
            "accepted_species").state
    return pd.DataFrame(scored).astype(float)


def _joint(arr):
    n=len(arr)
    counts=np.asarray(list(Counter(map(tuple,arr)).values()),dtype=float)
    return float(1-np.square(counts/n).sum())


def _seed(name,trait,seed):
    return int.from_bytes(hashlib.sha256(f"{seed}:{name}:{trait}".encode()).digest()[:8],"big")


def _sampling(binary:np.ndarray,morph:np.ndarray,*,n:int,repeats:int,rng):
    scores=[]
    for _ in range(repeats):
        ids=rng.choice(len(binary),size=n,replace=False)
        p=float(binary[ids].mean())
        scores.append((p,2*p*(1-p),_joint(morph[ids])))
    return np.mean(scores,axis=0)


def _block_ci(frame,key,draws,seed):
    b=frame.groupby("spatial_block")[key].mean().dropna().to_numpy(float)
    if len(b)<5:return (None,None,None)
    rng=np.random.default_rng(seed)
    boot=b[rng.integers(len(b),size=(draws,len(b)))].mean(axis=1)
    return tuple(map(float,(b.mean(),np.quantile(boot,.025),np.quantile(boot,.975))))


def audit_assurance_flora(
    original_ledger:pd.DataFrame, original_flora:pd.DataFrame,
    geo:pd.DataFrame, sources:pd.DataFrame, native_mainland:pd.DataFrame,
    guilds:pd.DataFrame, *,
    seed:int=20261113, sample_size:int=20, repeats:int=20,
    draws:int=1999,
)->dict:
    if sample_size!=20 or repeats<10 or draws<999:
        raise ValueError("Reproductive trait independent support frozen contract")
    morph=direct_triple(original_ledger)
    binary=direct_binary_assurance(original_ledger)
    required={
      "flora":(original_flora,{"island_id","accepted_species","origin_status"}),
      "geography":(geo,{"island_id","analysis_regime","spatial_block"}),
      "sources":(sources,{"island_id","source_mode","source_rank","entity_ID"}),
      "GIFT":(native_mainland,{"entity_ID","work_species"}),
      "guild":(guilds,{"island_id","channel_id","channel_state","source_state"}),
    }
    for k,(t,col) in required.items():
        if not col.issubset(t):raise ValueError("Original "+k+" contract missing")
    original=original_flora.loc[original_flora.accepted_species.isin(morph.index)].copy()
    if original.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Duplicate island plant")
    if not original.origin_status.isin(("native","endemic","introduced","unresolved")).all():
        raise ValueError("Unrecognized original origin status")
    counts=original.groupby("island_id").accepted_species.nunique()
    elig=set(counts[counts>=30].index)
    if len(elig)!=928:
        raise ValueError("Original 928 island Direct morphology universe changed")
    original=original.loc[original.island_id.isin(elig)]
    g=geo[["island_id","analysis_regime","spatial_block"]].drop_duplicates().set_index("island_id")
    if g.index.has_duplicates or not elig.issubset(g.index):
        raise ValueError("Original geography not available for all 928")
    donors=sources.loc[
       sources.island_id.isin(elig)&sources.source_mode.eq("geo_k5")&
       sources.source_rank.between(1,5),["island_id","source_rank","entity_ID"]]
    if donors.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicated fixed GIFT source")
    src=donors.sort_values(["island_id","source_rank"]).groupby("island_id").entity_ID.agg(tuple).to_dict()
    if any(len(src.get(i,()))!=5 or len(set(src[i]))!=5 for i in elig):
        raise ValueError("Missing five GIFT geography source entities")
    gift=(native_mainland.loc[native_mainland.work_species.isin(morph.index),
                              ["entity_ID","work_species"]].dropna().drop_duplicates()
          .groupby("entity_ID").work_species.agg(set).to_dict())
    ch=guilds[guilds.island_id.isin(elig)][["island_id","channel_id","channel_state","source_state"]].copy()
    if ch.duplicated(["island_id","channel_id"]).any():
        raise ValueError("Original channel IDs duplicated")
    if not (ch.loc[ch.channel_state.eq("disrupted"),"source_state"].eq("available")).all():
        raise ValueError("Qualified nondetection lacks source support")
    c=ch.pivot(index="island_id",columns="channel_id",values="channel_state")
    if not {"bombus","flower_visiting_birds"}.issubset(c.columns):
        raise ValueError("No original Bombus/bird opportunity")
    focal=set(c.index[c.bombus.eq("disrupted")&c.flower_visiting_birds.eq("retained")])
    if len(focal)!=21:
        raise ValueError("21 original focal candidate opportunity islands drifted")
    byisland={i:x for i,x in original.groupby("island_id")}
    rows=[]
    for island in sorted(elig):
        gi=byisland[island]
        src_names=set().union(*(gift.get(entity,set()) for entity in src[island]))
        src_names=src_names.intersection(morph.index)
        island_species=set(gi.accepted_species)
        for scope in ("all_recorded","source_native"):
            local=gi if scope=="all_recorded" else gi.loc[gi.origin_status.eq("native")]
            local_names=set(local.accepted_species)
            for trait in REPRO:
                values=binary[trait]
                known=set(values.dropna().index)
                isl=sorted(local_names & known)
                source=sorted(src_names & known)
                # Prevent different trait carrier count from changing effective
                # support under the metric itself. Never impute unknown states.
                if len(isl)<sample_size:
                    continue
                rand=np.random.default_rng(_seed(island,scope+trait,seed))
                local_morph=morph.loc[isl].to_numpy(str)
                local_bin=values.loc[isl].to_numpy(float)
                loc=_sampling(local_bin,local_morph,n=sample_size,repeats=repeats,rng=rand)
                src_val=None
                if len(source)>=sample_size:
                    rng=np.random.default_rng(_seed(island,"source"+trait,seed))
                    src_val=_sampling(values.loc[source].to_numpy(float),
                            morph.loc[source].to_numpy(str),
                            n=sample_size,repeats=repeats,rng=rng)
                d=dict(island_id=island,spatial_block=g.loc[island,"spatial_block"],
                    regime=g.loc[island,"analysis_regime"],source_signature=src[island][0],
                    trait=trait,scope=scope,n_island_known=len(isl),
                    n_source_known=len(source),candidate_bombus_bird=island in focal,
                    island_positive_share=float(loc[0]),island_binary_Q=float(loc[1]),
                    island_morph_joint_Q=float(loc[2]))
                for k,v in enumerate(("positive_share","binary_Q","morph_joint_Q")):
                    d["delta_"+v]=float(loc[k]-src_val[k]) if src_val is not None else np.nan
                rows.append(d)
    private=pd.DataFrame.from_records(rows)
    out=[]
    for (scope,trait,regime),q in private.groupby(["scope","trait","regime"]):
        if len(q)<10 or q.spatial_block.nunique()<5: continue
        pair=q.dropna(subset=["delta_positive_share"])
        rec=dict(scope=scope,trait=trait,regime=regime,n_islands=len(q),
            n_source_comparable=len(pair),
            n_spatial_blocks=q.spatial_block.nunique(),
            median_island_positive_share=float(q.island_positive_share.median()),
            median_island_reproductive_binary_Q=float(q.island_binary_Q.median()),
            median_island_morph_joint_Q_same_reproductive_carriers=float(q.island_morph_joint_Q.median()))
        if len(pair)>=10 and pair.spatial_block.nunique()>=5:
            for j,key in enumerate(("delta_positive_share","delta_binary_Q","delta_morph_joint_Q")):
                mean,lo,hi=_block_ci(pair,key,draws,seed+j)
                rec["block_mean_"+key]=mean
                rec["block_ci025_"+key]=lo
                rec["block_ci975_"+key]=hi
        out.append(rec)
    cand=[]
    for (scope,trait),q in private.loc[private.candidate_bombus_bird].groupby(["scope","trait"]):
        pair=q.dropna(subset=["delta_positive_share"])
        rec=dict(scope=scope,trait=trait,n_candidate_islands_with_20_known=len(q),
                 n_source_comparable=len(pair),n_spatial_blocks=q.spatial_block.nunique())
        if len(pair)>=10 and pair.spatial_block.nunique()>=5:
            for j,key in enumerate(("delta_positive_share","delta_morph_joint_Q")):
                rec["block_mean_"+key],rec["block_ci025_"+key],rec["block_ci975_"+key]=(
                    _block_ci(pair,key,draws,seed+100+j))
        cand.append(rec)
    return {
     "status":"DIRECT_REPRODUCTIVE_FUNCTIONAL_PROXIES_NOT_POLLINATION_SERVICE",
     "n_original_direct_triple_flora_islands":928,
     "n_original_focal_Bombus_nondetection_bird_detection_candidates":21,
     "n_species_with_direct_assurance_by_trait":{
        t:int(binary[t].notna().sum()) for t in REPRO},
     "reproductive_assurance_functional_assemblage_by_region":out,
     "candidate_potential_pollinator_shift_assurance_support":cand,
     "carrier_specific_support_min":sample_size,
     "rarefactions":repeats,
     "same_carrier_morphology_and_assurance_scoring":True,
     "unknown_reproductive_trait_is_not_state_zero":True,
     "source_geography_is_not_an_ancestor":True,
     "no_actual_pollinator_loss_or_shift_measured":True,
     "no_species_island_or_source_level_data_exported":True,
    }
