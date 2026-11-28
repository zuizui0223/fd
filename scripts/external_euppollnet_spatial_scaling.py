#!/usr/bin/env python3
"""Post-discovery spatial scaling of external visitor-channel retention."""
import argparse, collections, json, math, pathlib, tempfile
import numpy as np
import pandas as pd
from external_euppollnet_bee_loss import BASE,SOURCES,download
from external_euppollnet_channel_retention import prepare,collect,summarize
from external_euppollnet_curveball import estimate,summary
from external_floral_replication_intake import name

SEED=20261128
N_BOOT=9999
KEY=["study","year","plant"]

def km(a,b):
    la,lo,lb,lp=map(math.radians,(*a,*b))
    h=math.sin((lb-la)/2)**2+math.cos(la)*math.cos(lb)*math.sin((lp-lo)/2)**2
    return 6371.0088*2*math.asin(math.sqrt(min(1.,max(0.,h))))

def fixed_weights(q):
    # Equal species cells per year, equal years per study, equal studies.
    if q.empty:return np.array([])
    cells=q.groupby(["study","year"]).plant.transform("size").to_numpy()
    years=q[["study","year"]].drop_duplicates().groupby("study").size()
    return 1.0/(cells*q.study.map(years).to_numpy())

def support(q):
    return {"plant_study_year_cells":len(q),"species":int(q.plant.nunique()),
            "genera":int(q.genus.nunique()),"studies":int(q.study.nunique()),
            "study_years":int(q[["study","year"]].drop_duplicates().shape[0])}

def infer(q,value):
    if q.empty:return {"status":"no_eligible_cells","gate":False,"supported":False}
    w=fixed_weights(q);v=q[value].to_numpy(float)
    studies=sorted(q.study.unique());genera=sorted(q.genus.unique())
    si=q.study.map({k:i for i,k in enumerate(studies)}).to_numpy()
    gi=q.genus.map({k:i for i,k in enumerate(genera)}).to_numpy()
    ns,ng=len(studies),len(genera)
    a=np.zeros((ns,ng));b=np.zeros_like(a)
    np.add.at(a,(si,gi),w*v);np.add.at(b,(si,gi),w)
    assert np.allclose(b.sum(1),1)
    obs=float(a.sum()/b.sum())
    sup=support(q)
    gate=ns>=10 and ng>=20 and sup["species"]>=50
    out={"support":sup,"gate":bool(gate),"estimate":obs,"supported":False}
    if ns<2 or ng<2:return out
    rng=np.random.default_rng(SEED+10)
    S=rng.multinomial(ns,np.ones(ns)/ns,size=N_BOOT)
    G=rng.multinomial(ng,np.ones(ng)/ng,size=N_BOOT)
    def ratio(num,den):
        ok=den>1e-12
        return num[ok]/den[ok]
    sb=ratio(S@a.sum(1),S@b.sum(1))
    gb=ratio(G@a.sum(0),G@b.sum(0))
    # Pigeonhole product weighting, conditional fixed original cell weights.
    cb=ratio(np.einsum("ij,ij->i",S@a,G),np.einsum("ij,ij->i",S@b,G))
    for k,bo in [("study",sb),("genus",gb),("crossed_study_genus",cb)]:
        out[k+"_valid_draws"]=len(bo)
        out[k+"_bootstrap_ci"]=np.quantile(bo,[.025,.975]).tolist() if len(bo)>=.95*N_BOOT else None
    ls=ratio(a.sum()-a.sum(1),b.sum()-b.sum(1))
    lg=ratio(a.sum()-a.sum(0),b.sum()-b.sum(0))
    out["leave_study_range"]=[float(ls.min()),float(ls.max())] if len(ls)==ns else None
    out["leave_genus_range"]=[float(lg.min()),float(lg.max())] if len(lg)==ng else None
    ci=out["crossed_study_genus_bootstrap_ci"];ls=out["leave_study_range"];lg=out["leave_genus_range"]
    out["supported"]=bool(gate and ci and ls and lg and ci[0]>0 and ls[0]>0 and lg[0]>0)
    return out

def assign_genus(rows):
    r=rows.copy()
    r["exact_name"]=r.plant.map(name)
    valid=r.dropna(subset=["exact_name"])
    bad=set(valid.groupby("exact_name").plant.nunique().loc[lambda z:z>1].index)
    r=r.loc[r.exact_name.notna()&~r.exact_name.isin(bad)].copy()
    r["genus"]=r.exact_name.str.split().str[0]
    return r,{"input_pairs":len(rows),"exact_name_pairs":len(r),"colliding_exact_names":len(bad)}

def level_cells(r):
    return r.groupby(KEY+["genus"],as_index=False).residual.mean()

def spatial_cells(r):
    out=[];counts=collections.Counter()
    for key,g in r.groupby(KEY+["genus"],sort=True):
        counts["candidate_cells"]+=1
        g=g.loc[g.distance_km.gt(0)]
        if len(g)<3 or g.distance_km.nunique()<3:
            counts["fewer_than_three_distances"]+=1;continue
        x=np.log2(g.distance_km.to_numpy(float))
        if np.ptp(x)<1:
            counts["less_than_twofold_distance_span"]+=1;continue
        xc=x-x.mean();den=np.dot(xc,xc)
        if den<=0:continue
        slopes={k:float(np.dot(xc,g[k].to_numpy(float))/den) for k in ("residual","observed","null")}
        assert abs(slopes["residual"]-(slopes["observed"]-slopes["null"]))<1e-10
        out.append(dict(zip(KEY+["genus"],key))|{"slope":slopes["residual"],"observed_slope":slopes["observed"],"null_slope":slopes["null"],"pairs":len(g)})
    return pd.DataFrame(out),dict(counts)

def panel(rows):
    if rows.empty:return {"status":"no_rows"}
    r,audit=assign_genus(rows)
    if r.empty:return {"name_audit":audit,"status":"no_exact_names"}
    r["distance_km"]=[km(a,b) for a,b in zip(r.site_a,r.site_b)]
    lev=level_cells(r);s,census=spatial_cells(r)
    p={"name_audit":audit,"distance_km_quantiles":np.quantile(r.distance_km,[0,.1,.25,.5,.75,.9,1]).tolist(),
       "level":infer(lev,"residual"),"spatial_census":census,
       "distance_slope_per_doubling":infer(s,"slope")}
    if not s.empty:
        w=fixed_weights(s);p["distance_slope_components"]={k:float(np.average(s[k],weights=w)) for k in ("slope","observed_slope","null_slope")}
        p["spatial_pairs"]=int(s.pairs.sum())
    return p

def self_test():
    assert abs(km((0,0),(0,1))-111.195)<.01
    assert abs(km((0,179),(0,-179))-222.390)<.01
    q=pd.DataFrame({"study":["a"]*3+["b"]*3,"year":[2020]*6,"plant":["p"]*3+["q"]*3,
                    "genus":["g"]*3+["h"]*3,"distance_km":[1,2,4]*2,
                    "observed":[.2,.3,.4]*2,"null":[.2,.25,.3]*2,"residual":[0,.05,.1]*2})
    s,c=spatial_cells(q);assert np.allclose(s.slope,.05) and len(s)==2
    z=infer(s,"slope");assert abs(z["estimate"]-.05)<1e-12 and not z["gate"]
    assert np.allclose(fixed_weights(s),[1,1])
    json.dumps(z,allow_nan=False)
    print("FD_EXTERNAL_SPATIAL_SELF_TEST_PASS",flush=True)

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
        assert len(records)==3043 and abs(summarize(records)["study_equal"]["channel_residual"]-.039527962549189964)<1e-10
        rows,diag=estimate(records,d);parent=summary(rows)
        assert len(rows)==2537 and abs(parent["study_equal"]["residual"]-.04458227814019177)<1e-10
        primary=panel(rows);sensitivity=panel(rows.loc[rows.min_exclusive.ge(2)])
        assert abs(primary["level"]["estimate"]-.04458227814019177)<1e-10
        result={"source_pins":pins,"parent_replay":parent["study_equal"]["residual"],"primary":primary,
                "descriptive_both_exclusive_ge2":sensitivity,
                "scope":"Post-discovery, fixed-margin conditional null and original retained pairs replayed. Haversine separation; within plant-study-year slopes on log2 distance, >=3 distinct positive distances spanning >=2-fold. Equal plant cells within years, years within studies. Crossed study/genus product-weight bootstrap is an exploratory dependence audit with fixed original weights, not calibrated spatial inference. No island effect or causal/pollination claim."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False))
        print("FD_EXTERNAL_SPATIAL_RESULT",json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":main()
