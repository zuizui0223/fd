#!/usr/bin/env python3
"""Apply the frozen discovery month-matching audit to study-ID-disjoint validation."""
import argparse,json,pathlib,tempfile
import numpy as np
import pandas as pd
from external_eupolltrait_size_turnover import load_traits,report
from external_eupolltrait_size_qc_family import corrected_traits,network
from external_eupolltrait_disjoint_validation import scope
from external_eupolltrait_provenance import components
from external_eupolltrait_month_match import monthly_network,annual_common_support,panel,cell_means
from external_eupolltrait_intake import prepare_network
from external_euppollnet_bee_loss import BASE,SOURCES,download
from external_euppollnet_spatial_scaling import assign_genus,fixed_weights

def component_month(rows,mp):
    if rows.empty:return {"status":"no_rows","supported":False}
    r,audit=assign_genus(rows)
    if r.empty:return {"status":"no_exact_names","supported":False}
    q=cell_means(r,["size_gain"]);q["component"]=q.study.map(mp)
    assert q.component.notna().all()
    w=fixed_weights(q)
    comps=sorted(q.component.unique());gs=sorted(q.genus.unique())
    ci=q.component.map({x:i for i,x in enumerate(comps)}).to_numpy()
    gi=q.genus.map({x:i for i,x in enumerate(gs)}).to_numpy()
    A=np.zeros((len(comps),len(gs)));B=np.zeros_like(A)
    np.add.at(A,(ci,gi),w*q.size_gain.to_numpy());np.add.at(B,(ci,gi),w)
    est=float(A.sum()/B.sum())
    pairs=len(r[["study","year","site_a","site_b"]].drop_duplicates())
    ns,ng=A.shape;gate=ns>=10 and ng>=20 and q.plant.nunique()>=50 and pairs>=30
    out={"comparisons":len(r),"site_pairs":pairs,"plant_species":int(q.plant.nunique()),"plant_genera":ng,
         "studies":int(q.study.nunique()),"components":ns,"estimate":est,"support_gate":bool(gate),"supported":False}
    if ns<2 or ng<2:return out
    rng=np.random.default_rng(20261128+509);n=9999
    S=rng.multinomial(ns,np.ones(ns)/ns,size=n);G=rng.multinomial(ng,np.ones(ng)/ng,size=n)
    def ratio(a,b):
        keep=b>1e-12
        return a[keep]/b[keep]
    z=ratio(np.einsum("ij,ij->i",S@A,G),np.einsum("ij,ij->i",S@B,G))
    lc=ratio(A.sum()-A.sum(1),B.sum()-B.sum(1));lg=ratio(A.sum()-A.sum(0),B.sum()-B.sum(0))
    bounds=np.quantile(z,[.025,.975]).tolist() if len(z)>=.95*n else None
    out.update(crossed_component_plant_genus_ci=bounds,valid_draws=len(z),
      leave_component_range=[float(lc.min()),float(lc.max())] if len(lc)==ns else None,
      leave_plant_genus_range=[float(lg.min()),float(lg.max())] if len(lg)==ng else None,
      supported=bool(gate and bounds and bounds[0]>0 and len(lc)==ns and len(lg)==ng and lc.min()>0 and lg.min()>0))
    return out

def self_test():
    x=pd.DataFrame({"study":["a","a","b"],"year":[2020]*3,"plant":["Genus species"]*3,
        "month":[1,2,1],"size_gain":[1.,3.,4.],"site_a":[(0,0)]*3,"site_b":[(0,1)]*3})
    q=component_month(x,{"a":"a","b":"b"})
    assert abs(q["estimate"]-3.)<1e-12 and not q["support_gate"]
    print("FD_EUPOLLTRAIT_VALIDATION_MONTH_SELF_TEST_PASS",flush=True)

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
        disc={mp[s] for s in sa["discovery_study_ids"]}
        v,va=prepare_network(raw,md,eligible_study_year=keys)
        assert not set(v.Study_id)&set(sa["discovery_study_ids"])
        annual,ad=network(v,tab,"bee",tax);ar=report(annual)
        assert ar["comparisons"]==19618 and abs(ar["study_equal"]["size_gain"]-.10953097474723927)<1e-10
        rows,diag=monthly_network(v,tab,tax)
        joined=annual_common_support(rows,annual)
        common=int(joined.annual_size_gain.notna().sum()) if not joined.empty else 0
        excl=rows.loc[~rows.study.map(mp).isin(disc)] if not rows.empty else rows
        main=panel(rows)
        cluster=component_month(rows,mp)
        if not rows.empty:assert abs(cluster["estimate"]-main["study_equal"]["size_gain"])<1e-10
        out={"source_pins":pins,"network_source_pins":npins,"selection":selection,
            "scope_audit":sa,"network_name_audit":va,"metadata_components":ma,
            "annual_replay":{"comparisons":ar["comparisons"],"study_equal":ar["study_equal"]},
            "monthly_network":diag,"primary":main,"same_month_component_dependence":cluster,
            "annual_common_support_comparisons":common,"descriptive_annual_common_support":panel(joined,True),
            "descriptive_both_exclusive_ge2":panel(rows.loc[rows.min_exclusive.ge(2)] if not rows.empty else rows),
            "descriptive_discovery_component_disjoint":component_month(excl,mp),
            "scope":"Post-discovery validation-group month audit, same recorded-month protocol and gates as discovery. Study-ID disjointness not field-observation independence; component grouping is conservative metadata dependence audit. Recorded dates and within-month opportunity/effort are not independently certified. No causal or pollination effectiveness claim."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
        print("FD_EUPOLLTRAIT_VALIDATION_MONTH_RESULT",json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()
