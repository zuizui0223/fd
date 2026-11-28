#!/usr/bin/env python3
"""Audit visitor-channel retention with fixed-margin local network nulls."""
import argparse, collections, json, pathlib, tempfile
import numpy as np
import pandas as pd
from external_euppollnet_bee_loss import BASE,SOURCES,PLANT,VISITOR,download
from external_euppollnet_channel_retention import prepare,collect,summarize

DRAWS=399
MIN_ACCEPT=30
SEED=20261128

def trade(rows,rng):
    i,j=rng.choice(len(rows),2,replace=False)
    a,b=rows[i],rows[j]
    shared=a&b
    x,y=a-shared,b-shared
    if not x or not y:return False
    pool=sorted(x|y)
    chosen=set(rng.choice(pool,size=len(x),replace=False).tolist())
    if chosen==x:return False
    rows[i]=shared|chosen
    rows[j]=shared|(set(pool)-chosen)
    return True

def margin(rows):
    return tuple(map(len,rows)),collections.Counter(v for r in rows for v in r)

def simulate(rows,draws,rng):
    rows=[set(r) for r in rows]
    before=margin(rows)
    n=len(rows)
    moved=sum(trade(rows,rng) for _ in range(10*n))
    out=[];moves=0
    for _ in range(draws):
        moves+=sum(trade(rows,rng) for _ in range(2*n))
        assert margin(rows)==before,"degree margin drift"
        out.append([frozenset(x) for x in rows])
    return out,{"burn_changes":moved,"sample_changes":moves,"plants":n}

def matches(a,b,labels):
    u=a-b;v=b-a
    if not u or not v:return None
    ca=collections.Counter(labels[x] for x in u)
    cb=collections.Counter(labels[x] for x in v)
    return sum(x*cb[k] for k,x in ca.items())/(len(u)*len(v))

def conditional_pair(a,b,draw_a,draw_b,labels,min_accept):
    shared=len(a&b)
    obs=matches(a,b,labels)
    if obs is None:return None
    vals=[];halves=[[],[]]
    for i,(x,y) in enumerate(zip(draw_a,draw_b)):
        assert len(x)==len(a) and len(y)==len(b)
        if len(x&y)!=shared:continue
        q=matches(x,y,labels)
        assert q is not None
        vals.append(q);halves[int(i>=len(draw_a)//2)].append(q)
    if len(vals)<min_accept or min(map(len,halves))<10:
        return {"accepted":len(vals),"eligible":False}
    mean=float(np.mean(vals))
    return {"accepted":len(vals),"eligible":True,"observed":obs,"null":mean,
            "residual":obs-mean,"half_gap":float(np.mean(halves[0])-np.mean(halves[1]))}

def build_sites(d):
    sites={};labels={}
    for (study,year),g in d.groupby(["Study_id","year"],sort=True):
        bad=set(g.groupby(VISITOR).family.nunique().loc[lambda x:x>1].index)
        g=g.loc[~g[VISITOR].isin(bad)]
        ids={s:i for i,s in enumerate(sorted(g[VISITOR].unique()))}
        labels[study,year]={ids[v]:c for v,c in zip(g[VISITOR],g.channel)}
        for coord,h in g.groupby(["Latitude","Longitude"],sort=True):
            edge=h[[PLANT,VISITOR]].drop_duplicates()
            if edge[PLANT].nunique()<5 or edge[VISITOR].nunique()<5:continue
            p=edge.groupby(PLANT)[VISITOR].agg(set).to_dict()
            names=sorted(p)
            sites[study,year,coord]={"index":{s:i for i,s in enumerate(names)},"rows":[{ids[v] for v in p[s]} for s in names]}
    return sites,labels

def estimate(records,d,draws=DRAWS,min_accept=MIN_ACCEPT):
    sites,labels=build_sites(d)
    needed=set()
    for r in records.itertuples(index=False):
        needed.update([(r.study,r.year,r.site_a),(r.study,r.year,r.site_b)])
    rng=np.random.default_rng(SEED);bank={};diagnostics=[]
    for key in sorted(needed):
        samples,info=simulate(sites[key]["rows"],draws,rng)
        bank[key]=list(zip(*samples))
        diagnostics.append(info)
    rows=[];accepted=[];rejected=0
    for r in records.itertuples(index=False):
        ka=(r.study,r.year,r.site_a);kb=(r.study,r.year,r.site_b)
        ia=sites[ka]["index"][r.plant];ib=sites[kb]["index"][r.plant]
        a=sites[ka]["rows"][ia];b=sites[kb]["rows"][ib]
        z=conditional_pair(a,b,bank[ka][ia],bank[kb][ib],labels[r.study,r.year],min_accept)
        assert z is not None
        accepted.append(z["accepted"])
        if not z["eligible"]:
            rejected+=1;continue
        assert abs(z["observed"]-r.channel_observed)<1e-12
        rows.append({"study":r.study,"year":r.year,"plant":r.plant,"site_a":r.site_a,"site_b":r.site_b,
                     "min_exclusive":r.min_exclusive,"original_residual":r.channel_residual,**z})
    return pd.DataFrame(rows),{"sites":len(needed),"frozen_sites":sum(x["sample_changes"]==0 for x in diagnostics),
          "median_site_sample_changes":float(np.median([x["sample_changes"] for x in diagnostics])) if diagnostics else None,
          "parent_pairs":len(records),"accepted_pairs":len(rows),"rejected_pairs":rejected,
          "accepted_draw_quantiles":np.quantile(accepted,[0,.25,.5,.75,1]).tolist() if accepted else []}

def summary(d):
    if d.empty:return {"support_gate":False,"primary_supported":False,"status":"no_eligible_pairs"}
    values=["observed","null","residual","original_residual","half_gap"]
    by=d.groupby(["study","year","plant"])[values].mean().groupby(level=["study","year"]).mean().groupby(level="study").mean()
    n=len(by)
    support={"comparisons":len(d),"species":int(d.plant.nunique()),"studies":n,
             "site_pairs":int(d[["study","year","site_a","site_b"]].drop_duplicates().shape[0])}
    gate=n>=10 and support["species"]>=50 and support["site_pairs"]>=30
    r={"support":support,"support_gate":bool(gate),"study_equal":{k:float(v) for k,v in by.mean().items()},
       "median_abs_pair_half_gap":float(d.half_gap.abs().median()),"primary_supported":False}
    if n>=2:
        rng=np.random.default_rng(SEED+1)
        v=by.residual.to_numpy();boot=v[rng.integers(n,size=(9999,n))].mean(1)
        loo=(v.sum()-v)/(n-1)
        r["study_bootstrap_ci"]=np.quantile(boot,[.025,.975]).tolist()
        r["leave_study_range"]=[float(loo.min()),float(loo.max())]
        r["primary_supported"]=bool(gate and r["study_bootstrap_ci"][0]>0 and loo.min()>0)
    return r

def self_test():
    rows=[{0,1,2},{1,2,3},{2,3,4},{3,4,5},{0,4,5}]
    sims,diag=simulate(rows,40,np.random.default_rng(123))
    assert all(margin(a)==margin(rows) for a in sims) and diag["sample_changes"]>0
    labels={0:"a",1:"a",2:"b",3:"a",4:"b"}
    a={0,2};b={1,2}
    da=[{0,2},{0,4},{2,3}]*20
    db=[{1,2},{1,3},{1,2}]*20
    q=conditional_pair(a,b,da,db,labels,30)
    assert q["eligible"] and q["accepted"]==40 and q["observed"]==1 and q["null"]==1
    assert not conditional_pair(a,b,da,db,labels,41)["eligible"]
    print("FD_EXTERNAL_CURVEBALL_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true")
    a=p.parse_args();self_test()
    if a.self_test:return
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);pins={}
        for f,(sha,blob) in SOURCES.items():pins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
        raw=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False)
        assert len(raw)==623476
        d=prepare(raw,pd.read_csv(w/"Metadata.csv",dtype=str));records,census=collect(d)
        orig=summarize(records)
        assert len(records)==3043 and abs(orig["study_equal"]["channel_residual"]-.039527962549189964)<1e-10
        rows,diag=estimate(records,d)
        result={"source_pins":pins,"parent_replay":orig["study_equal"]["channel_residual"],"diagnostics":diag,"primary":summary(rows),
                "sensitivity_min_two_exclusive":summary(rows.loc[rows.min_exclusive.ge(2)] if not rows.empty else rows),
                "scope":"Post-discovery conditional null audit, no new independent observations. Both local binary degree margins, species identities and observed cross-site shared count preserved; 399 curveball samples with fixed burn/thinning, minimum 30 accepted per pair and at least 10 each half. Chain mixing not certified. No island-specific effect or effectiveness."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False))
        print("FD_EXTERNAL_CURVEBALL_RESULT",json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":main()
