#!/usr/bin/env python3
"""Recorded calendar-month matching of cross-site bee ITD turnover."""
import argparse, collections, itertools, json, pathlib, tempfile
import numpy as np
import pandas as pd
from external_eupolltrait_size_turnover import load_traits,report
from external_eupolltrait_size_qc_family import corrected_traits,network,conditional,decomposition_report
from external_eupolltrait_intake import prepare_network
from external_euppollnet_bee_loss import BASE,SOURCES,PLANT,download
from external_euppollnet_curveball import simulate
from external_euppollnet_spatial_scaling import assign_genus,infer,fixed_weights

SEED=20261128
DRAWS=399

def calendar_dates(q):
    # Require ISO calendar day (no year-only/month-only date or invented day).
    # The underlying publisher may have already imputed dates; cannot rule out.
    text=q.Date.astype(str)
    iso=text.str.fullmatch(r"\d{4}-\d{2}-\d{2}")
    dt=pd.to_datetime(text.where(iso),format="%Y-%m-%d",errors="coerce")
    good=dt.notna()&dt.dt.year.eq(q.year)
    out=q.loc[good].copy()
    out["month"]=dt.loc[good].dt.month.astype(int)
    out["record_day"]=dt.loc[good].dt.day.astype(int)
    audit={"input_rows":len(q),"calendar_day_rows":len(out),"dropped_rows":len(q)-len(out),
           "months":{str(k):int(v) for k,v in out.month.value_counts().sort_index().items()},
           "day_of_month_counts":{str(k):int(v) for k,v in out.record_day.value_counts().sort_index().items()},
           "year_month_study_cells":len(out[["Study_id","year","month"]].drop_duplicates()),
           "first_day_fraction":float(out.record_day.eq(1).mean()) if len(out) else None}
    return out,audit

def monthly_network(q,tab,tax):
    fm=tax.set_index("name").family.to_dict()
    q=q.loc[q.group.eq("bee")&q.name.isin(tab)&q.name.map(fm).eq(q.family)].copy()
    q,audit=calendar_dates(q)
    ids={k:i for i,k in enumerate(sorted(q.name.unique()))}
    fam={k:i for i,k in enumerate(sorted({fm[v] for v in ids}))}
    T=np.array([tab[k] for k in ids])
    F=np.array([fam[fm[k]] for k in ids],dtype=int)
    rng=np.random.default_rng(SEED+101);diag=collections.Counter();accepted=[];rows=[]
    for (study,year,month),g in q.groupby(["Study_id","year","month"],sort=True):
        sites={}
        for coord,s in g.groupby(["Latitude","Longitude"],sort=True):
            diag["candidate_site_months"]+=1
            # Require at least two recorded days to avoid single-date records;
            # not a guarantee of equal sampling intensity or actual dates.
            if s.record_day.nunique()<2:
                diag["one_record_day_site_months"]+=1;continue
            e=s[[PLANT,"name"]].drop_duplicates()
            if e[PLANT].nunique()<5 or e.name.nunique()<5:
                diag["small_site_months"]+=1;continue
            p=e.groupby(PLANT).name.agg(set).to_dict();names=sorted(p)
            rr=[{ids[v] for v in p[k]} for k in names]
            sims,info=simulate(rr,DRAWS,rng)
            sites[coord]=({k:i for i,k in enumerate(names)},rr,list(zip(*sims)))
            diag["eligible_site_months"]+=1
            diag["frozen_site_months"]+=int(info["sample_changes"]==0)
        for (ca,(ia,ra,sa)),(cb,(ib,rb,sb)) in itertools.combinations(sites.items(),2):
            for p in sorted(ia.keys()&ib.keys()):
                a,b=ra[ia[p]],rb[ib[p]]
                if not a-b or not b-a:continue
                diag["candidate_pairs"]+=1
                z=conditional(a,b,sa[ia[p]],sb[ib[p]],T,F,len(fam))
                accepted.append(z["accepted"])
                if not z["eligible"]:
                    diag["rejected_conditional_support"]+=1;continue
                rows.append(dict(z,study=study,year=int(year),month=int(month),plant=p,
                                 site_a=ca,site_b=cb,min_exclusive=min(len(a-b),len(b-a))))
    return pd.DataFrame(rows),{"calendar_audit":audit,"linked_visitor_species":len(ids),**dict(diag),
         "accepted_draw_quantiles":np.quantile(accepted,[0,.25,.5,.75,1]).tolist() if accepted else []}

def annual_common_support(monthly,annual):
    keys=["study","year","plant","site_a","site_b"]
    if monthly.empty:return monthly
    assert not annual.duplicated(keys).any()
    a=annual[keys+["size_gain"]].rename(columns={"size_gain":"annual_size_gain"})
    # Monthly comparisons can have different exclusive taxa, and some
    # monthly turnover vanishes when annual visits are pooled.
    return monthly.merge(a,on=keys,how="left",validate="many_to_one")

def cell_means(r,cols):
    # Equal spatial pairs within plant-month, then equal observed months
    # within plant-year. Annual effects on the matched keys use identical
    # weights. Avoid treating study-month as an independent study.
    keys=["study","year","plant","genus"]
    return r.groupby(keys+["month"],as_index=False)[cols].mean().groupby(keys,as_index=False)[cols].mean()

def panel(rows,paired=False):
    if rows.empty:return {"status":"no_eligible_pairs","supported":False}
    r,audit=assign_genus(rows)
    if r.empty:return {"status":"no_exact_focal_names","supported":False}
    cols=["size_obs","size_null","size_gain","family_gain"]
    if paired:
        r=r.loc[r.annual_size_gain.notna()].copy()
        if r.empty:return {"status":"no_annual_common_support","supported":False}
        r["monthly_minus_annual"]=r.size_gain-r.annual_size_gain
        cols+=["annual_size_gain","monthly_minus_annual"]
    q=cell_means(r,cols);w=fixed_weights(q)
    means={k:float(np.average(q[k],weights=w)) for k in cols}
    spatial_pairs=len(r[["study","year","site_a","site_b"]].drop_duplicates())
    site_month_pairs=len(r[["study","year","month","site_a","site_b"]].drop_duplicates())
    tests={k:infer(q,k) for k in ("size_gain","family_gain")+(("annual_size_gain","monthly_minus_annual") if paired else ())}
    return {"name_audit":audit,"comparisons":len(r),"site_pairs":spatial_pairs,
            "site_month_pairs":site_month_pairs,"focal_plant_month_cells":len(r[["study","year","month","plant"]].drop_duplicates()),
            "study_equal":means,"inference":tests,
            "supported":bool(spatial_pairs>=30 and tests["size_gain"]["supported"])}

def self_test():
    q=pd.DataFrame({"Date":["2020-02-29","2020-02-30","2020-02","2020-01-01","2021-01-01"],
        "year":[2020]*5,"Study_id":["s"]*5})
    d,a=calendar_dates(q);assert len(d)==2 and a["first_day_fraction"]==.5
    r=pd.DataFrame({"study":["s"]*3,"year":[2020]*3,"plant":["P a"]*3,"genus":["p"]*3,"month":[1,1,2],"x":[0.,2.,3.]})
    assert cell_means(r,["x"]).x.iloc[0]==2.
    print("FD_EUPOLLTRAIT_MONTH_MATCH_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    tax,m,meta,pins=load_traits();tab,selection=corrected_traits(m,"bee")
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);npins={}
        for f,(sha,blob) in SOURCES.items():npins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
        raw=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False);assert len(raw)==623476
        d,na=prepare_network(raw,pd.read_csv(w/"Metadata.csv",dtype=str))
        annual,ad=network(d,tab,"bee",tax)
        ar=report(annual)
        assert len(annual)==1261 and abs(ar["study_equal"]["size_gain"]-.11064886953661242)<1e-10
        rows,diag=monthly_network(d,tab,tax)
        joined=annual_common_support(rows,annual)
        matched=int(joined.annual_size_gain.notna().sum()) if not joined.empty else 0
        out={"source_pins":pins,"network_source_pins":npins,"selection":selection,"network_name_audit":na,
             "annual_replay":ar,"monthly_network":diag,"primary":panel(rows),
             "annual_common_support_comparisons":matched,"descriptive_annual_common_support":panel(joined,True),
             "descriptive_both_exclusive_ge2":panel(rows.loc[rows.min_exclusive.ge(2)] if not rows.empty else rows),
             "scope":"Post-discovery recorded-month matching, not independent replication. New monthly subweb and local fixed-degree/shared-count null. Exact annual replay and descriptive matched annual spatial keys. Recorded date precision and sampling coverage not certified; month matching does not eliminate within-month phenology or unequal effort. No actual pollination, causal mechanism or island claim."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
        print("FD_EUPOLLTRAIT_MONTH_MATCH_RESULT",json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()
