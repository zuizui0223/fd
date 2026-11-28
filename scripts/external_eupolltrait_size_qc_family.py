#!/usr/bin/env python3
"""Source-defined bee F/worker ITD, hoverfly QC, and exact family-pair size accounting."""
import argparse, collections, itertools, json, math, pathlib, tempfile
import numpy as np
import pandas as pd
from external_eupolltrait_size_turnover import (load_traits,parse_range,trait_table,report,
    BEE_TRAIT,FLY_TRAIT,SEED,DRAWS,MIN_ACCEPT)
from external_eupolltrait_intake import prepare_network
from external_euppollnet_bee_loss import BASE,SOURCES,PLANT,download
from external_euppollnet_curveball import simulate
from external_euppollnet_spatial_scaling import infer,assign_genus,fixed_weights

FLY_MIN=0.5
FLY_MAX=50.0
MIN_STRATUM_DRAWS=10

def corrected_traits(m,group):
    if group=="bee":
        raw=m.loc[m.group.eq(group)&m.measurementType.eq(BEE_TRAIT)].copy()
        x=raw.loc[raw.sex.str.strip().str.casefold().eq("f")&raw.originalColumn.eq("ITD_F_mm")]
        count=collections.Counter();tab={}
        for name,q in x.groupby("name",sort=True):
            vals=set()
            for v in q.itertuples():
                if str(v.lifeStage).strip().casefold() not in {"adult","adults","imago","na","nan","none",""}:
                    count["stage_excluded"]+=1;continue
                if str(v.measurementUnit).strip().casefold() not in {"mm","na","nan",""}:
                    count["unit_excluded"]+=1;continue
                y=parse_range(v.measurementValue)
                if y is None or y[0]!=y[1]:
                    count["non_scalar_or_invalid"]+=1;continue
                vals.add(y[0])
            if len(vals)!=1:
                count["unresolved_or_conflicting_species"]+=1;continue
            tab[name]=math.log2(next(iter(vals)))
        old,_,old_audit=trait_table(m,group)
        assert set(old).issubset(tab) and all(old[k]==tab[k] for k in old)
        return tab,{"raw_rows":len(raw),"exact_f_rows":len(x),"excluded":dict(count),
                    "original_subset_species":len(old),"restored_source_f_species":len(set(tab)-set(old)),
                    "retained_species":len(tab),"range_mm":[2**min(tab.values()),2**max(tab.values())] if tab else None,
                    "definition":"Published female ITD, explicitly workers for social bees; exact F/ITD_F_mm. Species-level measurements, visit sex unknown."}
    raw=m.loc[m.group.eq(group)&m.measurementType.eq(FLY_TRAIT)].copy()
    old,bounds,audit=trait_table(m,group)
    tab={k:v for k,v in old.items() if FLY_MIN<=bounds[k][0]<=bounds[k][1]<=FLY_MAX}
    invalid=set(old)-set(tab)
    # Report aggregate raw value shapes, never guess what a corrupted number meant.
    yr=sum(1900<=bounds[k][0]<=2100 or 1900<=bounds[k][1]<=2100 for k in invalid)
    bad=raw.loc[raw.name.isin(invalid),"measurementValue"].value_counts()
    shapes=collections.Counter()
    for v in raw.measurementValue:
        a=parse_range(v)
        shapes["unparsed" if a is None else ("scalar" if a[0]==a[1] else "range")]+=1
    return tab,{"original_parsed_species":len(old),"quarantined_species":len(invalid),
                "year_like_quarantined_species":yr,"quarantined_value_counts":{str(k):int(v) for k,v in bad.items()},
                "raw_format_counts":dict(shapes),"guard_mm":[FLY_MIN,FLY_MAX],"retained_species":len(tab),
                "retained_range_mm":[min(bounds[k][0] for k in tab),max(bounds[k][1] for k in tab)] if tab else None,
                "definition":"Post-discovery QC guard; does not certify plausible values as valid. No attempt to restore suspected spreadsheet corruption."}

def stratified(a,b,T,F,K):
    u=a-b;v=b-a
    if not u or not v:return None
    i=np.array(sorted(u));j=np.array(sorted(v))
    fi=F[i,None];fj=F[j][None,:]
    code=(np.minimum(fi,fj)*K+np.maximum(fi,fj)).ravel()
    dist=np.abs(T[i,None]-T[j][None,:]).ravel()
    return np.bincount(code,minlength=K*K),np.bincount(code,weights=dist,minlength=K*K)

def decompose(c,d,C,D,hits,n,K):
    # Each draw has the same exclusive sizes, so pooled cross-pair weights
    # equal the mean weights of the retained null draws.
    N=int(c.sum())
    assert N>0 and int(C.sum())==N*n
    p=c/N;q=C/(N*n)
    observed=float(d.sum()/N);expected=float(D.sum()/(N*n))
    diag=np.arange(K)*(K+1)
    result={"size_obs":observed,"size_null":expected,"size_gain":expected-observed,
            "family_obs":float(p[diag].sum()),"family_null":float(q[diag].sum()),
            "family_gain":float(p[diag].sum()-q[diag].sum()),
            "observed_family_strata":int((c>0).sum()),
            "low_support_observed_family_mass":float(p[(c>0)&(hits<MIN_STRATUM_DRAWS)].sum())}
    if np.any((c>0)&(hits<MIN_STRATUM_DRAWS)):
        result["family_standardization_eligible"]=False
        return result
    mu=np.divide(D,C,out=np.zeros_like(D,dtype=float),where=C>0)
    standardized=float(np.dot(p,mu))
    comp=expected-standardized;within=standardized-observed
    assert abs((comp+within)-(expected-observed))<1e-10
    result.update(family_standardization_eligible=True,size_null_at_observed_family_mix=standardized,
                  family_composition_gain=comp,family_standardized_size_gain=within)
    return result

def conditional(a,b,sa,sb,T,F,K):
    original=stratified(a,b,T,F,K)
    if original is None:return None
    c,d=original;h=len(a&b);C=np.zeros(K*K,dtype=np.int64);D=np.zeros(K*K)
    hits=np.zeros(K*K,dtype=np.int64);halves=[0,0]
    for k,(x,y) in enumerate(zip(sa,sb)):
        assert len(x)==len(a) and len(y)==len(b)
        if len(x&y)!=h:continue
        cc,dd=stratified(x,y,T,F,K)
        assert cc.sum()==c.sum()
        C+=cc;D+=dd;hits+=cc>0;halves[int(k>=len(sa)//2)]+=1
    n=sum(halves)
    if n<MIN_ACCEPT or min(halves)<10:return {"accepted":n,"eligible":False}
    return dict(decompose(c,d,C,D,hits,n,K),accepted=n,eligible=True)

def network(d,tab,group,tax,pair_fn=conditional):
    fm=tax.set_index("name").family.to_dict()
    q=d.loc[d.group.eq(group)&d.name.isin(tab)&d.name.map(fm).eq(d.family)].copy()
    ids={k:i for i,k in enumerate(sorted(q.name.unique()))}
    families={k:i for i,k in enumerate(sorted({fm[k] for k in ids}))}
    T=np.array([tab[k] for k in ids]);F=np.array([families[fm[k]] for k in ids],dtype=int)
    rng=np.random.default_rng(SEED);rows=[];diag=collections.Counter();acc=[]
    for (study,year),g in q.groupby(["Study_id","year"],sort=True):
        sites={}
        for coord,s in g.groupby(["Latitude","Longitude"],sort=True):
            edges=s[[PLANT,"name"]].drop_duplicates()
            if edges[PLANT].nunique()<5 or edges.name.nunique()<5:continue
            p=edges.groupby(PLANT).name.agg(set).to_dict();names=sorted(p)
            z=[{ids[v] for v in p[k]} for k in names]
            sims,info=simulate(z,DRAWS,rng)
            sites[coord]=({k:i for i,k in enumerate(names)},z,list(zip(*sims)))
            diag["eligible_sites"]+=1;diag["frozen_sites"]+=int(info["sample_changes"]==0)
        for (ca,(ia,ra,sa)),(cb,(ib,rb,sb)) in itertools.combinations(sites.items(),2):
            for p in sorted(ia.keys()&ib.keys()):
                a,b=ra[ia[p]],rb[ib[p]]
                if not a-b or not b-a:continue
                diag["candidate_pairs"]+=1
                z=pair_fn(a,b,sa[ia[p]],sb[ib[p]],T,F,len(families))
                acc.append(z["accepted"])
                if not z["eligible"]:
                    diag["rejected_conditional_support"]+=1;continue
                rows.append(dict(z,study=study,year=int(year),plant=p,site_a=ca,site_b=cb,min_exclusive=min(len(a-b),len(b-a))))
    return pd.DataFrame(rows),{"network_visit_rows":len(q),"visitor_species":len(ids),
        "studies_with_visits":int(q.Study_id.nunique()),"family_count":len(families),**dict(diag),
        "accepted_draw_quantiles":np.quantile(acc,[0,.25,.5,.75,1]).tolist() if acc else []}

def decomposition_report(rows):
    if rows.empty:return {"status":"no_eligible_pairs","supported":False}
    support=rows.loc[rows.family_standardization_eligible].copy()
    out={"parent_pairs":len(rows),"family_common_support_pairs":len(support),
         "excluded_pairs":len(rows)-len(support),
         "mean_observed_family_mass_under_gate":float(rows.low_support_observed_family_mass.mean()),
         "supported":False}
    if support.empty:return out
    r,audit=assign_genus(support)
    if r.empty:return out|{"name_audit":audit}
    cols=["size_obs","size_null","size_gain","size_null_at_observed_family_mix","family_composition_gain","family_standardized_size_gain","family_gain"]
    q=r.groupby(["study","year","plant","genus"],as_index=False)[cols].mean()
    w=fixed_weights(q)
    a={k:float(np.average(q[k],weights=w)) for k in cols}
    assert abs(a["size_gain"]-a["family_composition_gain"]-a["family_standardized_size_gain"])<1e-10
    pairs=int(r[["study","year","site_a","site_b"]].drop_duplicates().shape[0])
    tests={k:infer(q,k) for k in ("size_gain","family_composition_gain","family_standardized_size_gain")}
    return out|{"name_audit":audit,"comparisons":len(r),"site_pairs":pairs,"study_equal":a,
                "crossed_intervals":tests,"supported":bool(tests["family_standardized_size_gain"]["supported"] and pairs>=30)}

def self_test():
    T=np.array([0.,1.,3.,4.]);F=np.array([0,0,1,1])
    c,d=stratified({0,1},{2,3},T,F,2)
    assert c.sum()==4 and d.sum()==12
    # Construct a fixed null whose family composition is partly same-family.
    cc,dd=stratified({0,2},{1,3},T,F,2)
    C=cc*20;D=dd*20;hits=(cc>0)*20
    z=decompose(c,d,C,D,hits,20,2)
    assert z["family_standardization_eligible"] and abs(z["size_gain"]-(z["family_composition_gain"]+z["family_standardized_size_gain"]))<1e-12
    assert z["family_obs"]==0 and z["family_null"]==.5
    z=decompose(c,d,C,D,hits*0,20,2)
    assert not z["family_standardization_eligible"]
    # Constant families imply zero composition component, total = within.
    F=np.zeros(4,dtype=int);c,d=stratified({0,1},{2,3},T,F,1);cc,dd=stratified({0,2},{1,3},T,F,1)
    z=decompose(c,d,cc*20,dd*20,(cc>0)*20,20,1)
    assert abs(z["family_composition_gain"])<1e-12
    print("FD_EUPOLLTRAIT_QC_FAMILY_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    tax,m,metadata,pins=load_traits()
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);npins={}
        for f,(sha,blob) in SOURCES.items():npins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
        raw=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False);assert len(raw)==623476
        d,na=prepare_network(raw,pd.read_csv(w/"Metadata.csv",dtype=str))
        results={}
        for g in ("bee","hoverfly"):
            tab,audit=corrected_traits(m,g)
            print("FD_EUPOLLTRAIT_QC_SELECTION",g,json.dumps(audit,allow_nan=False),flush=True)
            rows,diag=network(d,tab,g,tax)
            sensitivity=rows.loc[rows.min_exclusive.ge(2)] if not rows.empty else rows
            results[g]={"selection":audit,"network":diag,"primary":report(rows),"both_exclusive_ge2":report(sensitivity),
                        "exact_family_pair_accounting":decomposition_report(rows),
                        "family_accounting_both_exclusive_ge2":decomposition_report(sensitivity)}
            print("FD_EUPOLLTRAIT_QC_FAMILY_PANEL",g,json.dumps(results[g],allow_nan=False),flush=True)
        out={"source_pins":pins,"network_source_pins":npins,"network_name_audit":na,"panels":results,
             "scope":"Post-outcome source-definition correction and QC; no independent replication. Original bee gate unchanged; fly descriptive only. Standardization is exact arithmetic on common family-pair support, not causal mediation. No actual effectiveness or island claim."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
        print("FD_EUPOLLTRAIT_QC_FAMILY_COMPLETE",json.dumps({"bee_primary_supported":results["bee"]["primary"].get("supported",False),"bee_family_standardized_supported":results["bee"]["exact_family_pair_accounting"].get("supported",False),"hoverfly_descriptive_only":True}),flush=True)
if __name__=="__main__":main()
