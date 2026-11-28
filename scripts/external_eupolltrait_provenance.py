#!/usr/bin/env python3
"""Contributor/DOI cluster dependence audit of study-disjoint bee ITD validation."""
import argparse, collections, json, pathlib, re, tempfile
import numpy as np
import pandas as pd
from external_eupolltrait_size_turnover import load_traits,report
from external_eupolltrait_size_qc_family import corrected_traits,network
from external_eupolltrait_disjoint_validation import scope
from external_eupolltrait_intake import prepare_network
from external_euppollnet_bee_loss import BASE,SOURCES,download
from external_emln_spatial_audit import doi_set
from external_euppollnet_spatial_scaling import assign_genus,fixed_weights

SEED=20261128
DRAWS=9999
ORCID=re.compile(r"(?<!\d)\d{4}-\d{4}-\d{4}-\d{3}[\dXx](?!\d)")

def valid_orcid(x):
    # ISO 7064 MOD 11-2 check digit.
    s=x.replace("-","").upper()
    t=0
    for c in s[:-1]:t=(t+int(c))*2
    v=(12-t%11)%11
    return s[-1]==("X" if v==10 else str(v))

def components(meta):
    # Public provenance identifiers used in-memory, never printed/exported.
    # Connect all publisher studies, including excluded scope studies, to
    # conservatively account for transitive DOI/ORCID links.
    ids=sorted(meta.Study_id.unique());parent={s:s for s in ids};tokens={}
    def root(s):
        while parent[s]!=s:
            parent[s]=parent[parent[s]];s=parent[s]
        return s
    def union(a,b):
        a,b=root(a),root(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    invalid=0;no_orcid=0;no_doi=0;none=0
    for s,g in meta.groupby("Study_id",sort=True):
        ds=doi_set(g.DOI)
        allids={x.upper() for v in g.Orcids.fillna("") for x in ORCID.findall(v)}
        os={x for x in allids if valid_orcid(x)}
        invalid+=len(allids-os)
        no_orcid+=not bool(os);no_doi+=not bool(ds);none+=not (os or ds)
        for tok in sorted({("doi",x) for x in ds}|{("orcid",x) for x in os}):
            if tok in tokens:union(s,tokens[tok])
            else:tokens[tok]=s
    mp={s:root(s) for s in ids}
    # No raw contributor IDs, names, or emails in result.
    sizes=collections.Counter(mp.values())
    return mp,{"all_studies":len(ids),"components":len(sizes),"component_sizes":sorted(sizes.values(),reverse=True),
              "studies_without_valid_orcid":int(no_orcid),"studies_without_doi":int(no_doi),
              "studies_without_either":int(none),"invalid_orcid_study_token_count":invalid}

def cluster_infer(rows,mp,equal_component=False):
    if rows.empty:return {"status":"no_rows","supported":False}
    r,audit=assign_genus(rows)
    if r.empty:return {"status":"no_exact_names","supported":False}
    q=r.groupby(["study","year","plant","genus"],as_index=False)[["size_obs","size_null","size_gain"]].mean()
    q["component"]=q.study.map(mp);assert q.component.notna().all()
    w=fixed_weights(q)
    if equal_component:
        nc=q[["study","component"]].drop_duplicates().groupby("component").size()
        w=w/q.component.map(nc).to_numpy()
    groups=sorted(q.component.unique());genera=sorted(q.genus.unique())
    si=q.component.map({k:i for i,k in enumerate(groups)}).to_numpy()
    gi=q.genus.map({k:i for i,k in enumerate(genera)}).to_numpy()
    a=np.zeros((len(groups),len(genera)));b=np.zeros_like(a)
    np.add.at(a,(si,gi),w*q.size_gain.to_numpy())
    np.add.at(b,(si,gi),w)
    obs=float(a.sum()/b.sum())
    ns,ng=a.shape
    size=q[["study","component"]].drop_duplicates().groupby("component").size()
    pairs=int(r[["study","year","site_a","site_b"]].drop_duplicates().shape[0])
    gate=ns>=10 and ng>=20 and q.plant.nunique()>=50 and pairs>=30
    out={"name_audit":audit,"comparisons":len(r),"site_pairs":pairs,"study_count":int(q.study.nunique()),
        "component_count":ns,"plant_species":int(q.plant.nunique()),"plant_genera":ng,
        "contributing_component_sizes":sorted(map(int,size),reverse=True),
        "largest_component_weight_fraction":float(b.sum(1).max()/b.sum()),
        "study_equal_with_component_bootstrap":not equal_component,
        "estimate":obs,"observed":float(np.average(q.size_obs,weights=w)),
        "null":float(np.average(q.size_null,weights=w)),"gate":bool(gate),"supported":False}
    if ns<2 or ng<2:return out
    rng=np.random.default_rng(SEED+503)
    S=rng.multinomial(ns,np.ones(ns)/ns,size=DRAWS)
    G=rng.multinomial(ng,np.ones(ng)/ng,size=DRAWS)
    def ratio(x,y):
        ok=y>1e-12
        return x[ok]/y[ok]
    v=ratio(np.einsum("ij,ij->i",S@a,G),np.einsum("ij,ij->i",S@b,G))
    c=ratio(S@a.sum(1),S@b.sum(1))
    g=ratio(G@a.sum(0),G@b.sum(0))
    for label,z in (("crossed_component_plant_genus",v),("component",c),("plant_genus",g)):
        out[label+"_valid_draws"]=len(z)
        out[label+"_ci"]=np.quantile(z,[.025,.975]).tolist() if len(z)>=.95*DRAWS else None
    lc=ratio(a.sum()-a.sum(1),b.sum()-b.sum(1))
    lg=ratio(a.sum()-a.sum(0),b.sum()-b.sum(0))
    out["leave_component_range"]=[float(lc.min()),float(lc.max())] if len(lc)==ns else None
    out["leave_plant_genus_range"]=[float(lg.min()),float(lg.max())] if len(lg)==ng else None
    ci=out["crossed_component_plant_genus_ci"]
    out["supported"]=bool(gate and ci and ci[0]>0 and len(lc)==ns and len(lg)==ng and lc.min()>0 and lg.min()>0)
    return out

def self_test():
    assert valid_orcid("0000-0002-1825-0097") and not valid_orcid("0000-0002-1825-0098")
    m=pd.DataFrame({"Study_id":["a","b","c","d"],"DOI":["10.1234/a","10.1234/a","10.1234/c","NA"],
        "Orcids":["","0000-0002-1825-0097","0000-0002-1825-0097",""]})
    mp,a=components(m)
    assert mp["a"]==mp["b"]==mp["c"] and mp["d"]!=mp["a"]
    assert a["component_sizes"]==[3,1]
    r=pd.DataFrame({"study":["a","b","c","d"],"year":[2020]*4,
        "plant":["Genus species","Genus other","Other species","Other other"],
        "site_a":[(0,0)]*4,"site_b":[(0,1)]*4,
        "size_obs":[1.,1.,1.,1.],"size_null":[2.,2.,2.,3.],"size_gain":[1.,1.,1.,2.]})
    x=cluster_infer(r,mp);y=cluster_infer(r,mp,True)
    assert abs(x["estimate"]-1.25)<1e-12 and abs(y["estimate"]-1.5)<1e-12
    assert not x["gate"]
    print("FD_EUPOLLTRAIT_PROVENANCE_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    tax,m,meta,pins=load_traits();tab,selection=corrected_traits(m,"bee")
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);npins={}
        for f,(sha,blob) in SOURCES.items():npins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
        raw=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False);assert len(raw)==623476
        md=pd.read_csv(w/"Metadata.csv",dtype=str)
        keys,sa=scope(md);mp,ma=components(md)
        discovery_components={mp[s] for s in sa["discovery_study_ids"]}
        v_ids=set(keys.Study_id);linked={s for s in v_ids if mp[s] in discovery_components}
        ma.update(discovery_components=len(discovery_components),
            validation_candidate_studies=len(v_ids),validation_candidate_components=len({mp[s] for s in v_ids}),
            validation_candidates_linked_to_discovery=len(linked))
        print("FD_EUPOLLTRAIT_PROVENANCE_SCOPE",json.dumps(ma,allow_nan=False),flush=True)
        d,na=prepare_network(raw,md);annual,_=network(d,tab,"bee",tax);ar=report(annual)
        assert len(annual)==1261 and abs(ar["study_equal"]["size_gain"]-.11064886953661242)<1e-10
        v,va=prepare_network(raw,md,eligible_study_year=keys)
        rows,diag=network(v,tab,"bee",tax);vr=report(rows)
        assert vr["comparisons"]==19618 and abs(vr["study_equal"]["size_gain"]-.10953097474723927)<1e-10
        primary=cluster_infer(rows,mp)
        assert abs(primary["estimate"]-.10953097474723927)<1e-10
        disjoint=rows.loc[~rows.study.isin(linked)]
        # Novel focal genus defined from all eligible exact discovery turnover
        # comparisons, not a taxonomic roster or fitted effect.
        disc,_=assign_genus(annual);val,_=assign_genus(rows)
        novel=val.loc[~val.genus.isin(set(disc.genus))]
        novel_disjoint=novel.loc[~novel.study.isin(linked)]
        out={"source_pins":pins,"network_source_pins":npins,"metadata_audit":ma,
          "validation_replay":{"comparisons":vr["comparisons"],"size_gain":vr["study_equal"]["size_gain"]},
          "primary":primary,
          "discovery_component_disjoint":cluster_infer(disjoint,mp),
          "descriptive_equal_component":cluster_infer(rows,mp,True),
          "descriptive_novel_focal_genus":cluster_infer(novel,mp),
          "descriptive_novel_genus_and_component_disjoint":cluster_infer(novel_disjoint,mp),
          "scope":"Post-discovery provenance audit. Metadata DOI/valid ORCID connected components across all 54 studies, not proof of field-project identity or independence. Study-equal point estimate retained, component x focal genus product bootstrap. Discovery-linked exclusion is a prespecified stricter sensitivity, novel genus refers only to original eligible comparisons. No personal identifiers exported, no new source or causal claim."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
        print("FD_EUPOLLTRAIT_PROVENANCE_RESULT",json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()
