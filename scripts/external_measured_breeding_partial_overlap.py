#!/usr/bin/env python3
"""Post-calibration partial-overlap audit and conditional two-group contrasts.

Original 4-cell insect/wind x self-code model is non-estimable. Do NOT
reinterpret a two-group association as its confirmation. Do not require a
second GIFT trait to be present when that trait is not in the estimand.
Same-row, same-metric measured SC and AF; exact normalized binomial only.
No own-project data and no island inference.
"""
import argparse,collections,json,pathlib
import numpy as np
import pandas as pd
from external_measured_breeding_calibration import (
    book,name,clean,traits,names_index,cohort,describe,fetch,refs,
    METRICS,SEED,DRAWS,BOOK_SHA)

DETAIL_PIN="29e088a18c738c761a0d41f67966c4b455bbe3e9d7bf792b7eecac440b34ee52"
ESTIMANDS={
    "self_code_all":{"axis":"self","lo":"absent","hi":"present"},
    "self_code_insect_only":{"axis":"self","lo":"absent","hi":"present","guild":"insect"},
    "guild_only":{"axis":"guild","lo":"wind","hi":"insect"}
}
STRATA=("raw","family","measured_source","gift_source","family_and_both_sources")
# Gates set after original four-cell failure, before any two-group fit.
MIN_CELL=20;MIN_GENUS_PER_CELL=10;MIN_TOTAL_GENERA=30
MIN_STRATA=5

def detailed():
    rows,p=fetch("traits",traitid="3.6.3",biasref=0,biasderiv=0,startat=0,limit=10000)
    assert p["sha256"]==DETAIL_PIN
    d={}
    for r in rows:
        a=r.get("agreement")
        if a not in (None,"") and float(a)<.66:continue
        v=str(r["trait_value"]).strip().casefold()
        if v in ("","none","na","nan","null"):continue
        i=str(r["work_ID"]);assert i not in d
        d[i]=(v,refs(r.get("references")))
    assert len(d)==3119
    return d,p

def all_measured(d,metric,index,T):
    sc,af=METRICS[metric]
    x=d.copy()
    x["name"]=x.Species.map(name)
    x["sc"]=pd.to_numeric(x[sc],errors="coerce")
    x["af"]=pd.to_numeric(x[af],errors="coerce")
    ok=x.sc.between(0,1)&x.af.between(0,1)&x.name.notna()
    out=[];loss=collections.Counter()
    for n,g in x[ok].groupby("name",sort=True):
        matches=index.get(n,set())
        if len(matches)!=1:loss["not_unique_work_id"]+=1;continue
        i=next(iter(matches))
        fam=sorted({clean(v).casefold() for v in g.Family if clean(v)})
        if len(fam)!=1:loss["missing_or_conflicting_family"]+=1;continue
        src=tuple(sorted({clean(v) for v in g.Data_source if clean(v)}))
        if not src:loss["missing_workbook_source"]+=1;continue
        r={"name":n,"work_id":i,"family":fam[0],"genus":n.split()[0],
           "sc":float(g.sc.mean()),"af":float(g.af.mean()),"paired_rows":len(g),
           "measured_source":json.dumps(src)}
        for key in T:
            v,p=T[key].get(i,(None,()))
            r[key]=v if v is not None and p else "unrecorded"
            r[key+"_source"]=json.dumps(p)
        out.append(r)
    a=pd.DataFrame(out)
    assert len({r["name"] for r in out})==len({r["work_id"] for r in out})==len(out)
    return a,{"paired_named_rows":int(ok.sum()),"paired_named_species":int(x.loc[ok,"name"].nunique()),
              "unique_measured_matched_species":len(a),"attrition":dict(loss)}

def support(d,c):
    out={"species":len(d),"genera":int(d.genus.nunique()),"families":int(d.family.nunique()),
         "measured_source_patterns":int(d.measured_source.nunique()),"cells":[]}
    for v in (c["lo"],c["hi"]):
        a=d[d[c["axis"]]==v]
        out["cells"].append({"state":v,"species":len(a),"genera":int(a.genus.nunique()),
              "families":int(a.family.nunique()),"SC_mean":float(a.sc.mean()) if len(a) else None,
              "AF_mean":float(a.af.mean()) if len(a) else None,
              "SC_range":[float(a.sc.min()),float(a.sc.max())] if len(a) else None,
              "AF_ge_half":int((a.af>=.5).sum()),"AF_positive":int((a.af>0).sum())})
    return out

def gate(s):
    return s["genera"]>=MIN_TOTAL_GENERA and all(c["species"]>=MIN_CELL and c["genera"]>=MIN_GENUS_PER_CELL for c in s["cells"])

def weighted_solve(X,y,g,w):
    mass=np.bincount(g,weights=w)
    assert (mass>0).all()
    z=np.column_stack([X,y])
    mu=np.column_stack([np.bincount(g,weights=w*z[:,i])/mass for i in range(z.shape[1])])
    q=(z-mu[g])*np.sqrt(w[:,None])
    b,_,rank,sv=np.linalg.lstsq(q[:,:X.shape[1]],q[:,X.shape[1]],rcond=None)
    if rank!=X.shape[1] or not len(sv) or sv[-1]<=1e-10:return None
    return b

def fit(d,c,kind):
    s0=support(d,c)
    if not gate(s0):return {"status":"insufficient_cell_support","support":s0}
    lo=max(z["SC_range"][0] for z in s0["cells"])
    hi=min(z["SC_range"][1] for z in s0["cells"])
    if hi-lo<.1:return {"status":"insufficient_common_SC_range","support":s0,"SC_range":[lo,hi]}
    a=d[d.sc.between(lo,hi)].copy()
    axis=c["axis"]
    a["stratum"]="all" if kind=="raw" else a.family if kind=="family" else a.measured_source if kind=="measured_source" else a[axis+"_source"] if kind=="gift_source" else a.family+"|"+a.measured_source+"|"+a[axis+"_source"]
    if kind!="raw":
        keep=[]
        for k,b in a.groupby("stratum",sort=True):
            n=b[axis].value_counts()
            if min(n.get(c["lo"],0),n.get(c["hi"],0))>=2:keep.append(k)
        a=a[a.stratum.isin(keep)].copy()
    s=support(a,c);s["strata"]=int(a.stratum.nunique())
    if not gate(s) or (kind!="raw" and s["strata"]<MIN_STRATA):
        return {"status":"insufficient_matched_support","support":s,"SC_range":[lo,hi]}
    # A common slope is primary for the conditional descriptive difference;
    # heterogeneous SC slopes are a separately reported sensitivity at the
    # pooled, fixed mean SC in the retained cohort. Linear model is not a
    # biological extrapolation to all self-compatible plants.
    center=float(a.sc.mean())
    z=a.sc.to_numpy(float)-center
    h=(a[axis]==c["hi"]).to_numpy(float)
    y=a.af.to_numpy(float)
    strata=pd.factorize(a.stratum,sort=True)[0]
    genus=pd.factorize(a.genus,sort=True)[0]
    X=np.column_stack([z,h])
    Xi=np.column_stack([z,h,z*h])
    def model(X):
        p=weighted_solve(X,y,strata,np.ones(len(a)))
        if p is None:return {"status":"rank_deficient"}
        rng=np.random.default_rng(SEED);draw=[]
        for j in range(DRAWS):
            w=np.maximum(rng.exponential(size=int(genus.max())+1),1e-12)[genus]
            q=weighted_solve(X,y,strata,w)
            if q is not None:draw.append(q)
        return {"status":"estimated","point":p.tolist(),"draws":len(draw),
              "genus_weight_interval":np.quantile(draw,[.025,.975],axis=0).T.tolist() if len(draw)>=.9*DRAWS else None}
    main=model(X);interaction=model(Xi)
    return {"status":"exploratory_estimable" if main["status"]=="estimated" else "rank_deficient",
            "support":s,"SC_range":[lo,hi],"SC_center":center,
            "common_slope":{"order":["SC","high_minus_low"],**main},
            "different_slopes":{"order":["SC","high_minus_low_at_center","SC_by_high"],**interaction},
            "warning":"Conditional descriptive association, positive shared genus weights not calibrated sampling CI. Source/family strata are not verified independent experiments."}

def census(a):
    out={}
    for k in ("guild","self","detailed"):
        out[k]=[{"state":str(s),"species":len(b),"genera":int(b.genus.nunique()),"families":int(b.family.nunique())}
            for s,b in a.groupby(k,sort=True)]
    out["guild_by_self"]=[{"guild":str(g),"self":str(s),"species":len(b)} for (g,s),b in a.groupby(["guild","self"],sort=True)]
    return out

def self_test():
    d=pd.DataFrame({"Species":["Aaa_aaa","Aaa_aaa","Bbb_bbb","Ccc_ccc"],"Family":["F"]*4,"Data_source":["1"]*4,
        METRICS["FS"][0]:[.5,None,.8,.6],METRICS["FS"][1]:[None,.7,.4,.2]})
    idx={"aaa aaa":{"1"},"bbb bbb":{"2"},"ccc ccc":{"3"}}
    T={"guild":{"1":("insect",("1",)),"2":("wind",("1",))},"self":{"1":("present",("2",)),"3":("absent",("2",))},"detailed":{}}
    a,c=all_measured(d,"FS",idx,T);assert len(a)==2
    assert list(a.guild)==["wind","unrecorded"] and list(a.self)==["unrecorded","absent"]
    z=np.arange(50,dtype=float)/50;h=np.arange(50)%2;X=np.column_stack([z,h])
    y=.1+.2*z+.3*h
    assert np.allclose(weighted_solve(X,y,np.zeros(50,int),np.ones(50)),[.2,.3])
    print("EXTERNAL_MEASURED_PARTIAL_OVERLAP_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    d,source=book()
    wanted={n for s in d.Species if (n:=name(s)) is not None}
    T,pins=traits();T["detailed"],pins["detailed"]=detailed()
    index,species_pins=names_index(wanted)
    assert sum(len(v)==1 for v in index.values())==1501
    panels={}
    for metric in METRICS:
        x,attrition=all_measured(d,metric,index,T)
        orig,old_counts=cohort(d,metric,index,T)
        assert old_counts["eligible_species"]==({"FS":66,"SFL":30}[metric])
        assert (orig.guild=="insect").all()
        q=x[(x.guild=="insect")&x.self.isin(["present","absent"])]
        assert set(q.name)==set(orig.name) and len(q)==len(orig)
        # Reproduce the earlier 66/30-species intersection before expanding.
        q=q.set_index("name").sort_index();oo=orig.set_index("name").sort_index()
        assert np.allclose(q[["sc","af"]],oo[["sc","af"]],atol=1e-12)
        result={"attrition":attrition,"available_trait_counts":census(x),"estimands":{}}
        for key,c in ESTIMANDS.items():
            a0=x[x[c["axis"]].isin([c["lo"],c["hi"]])].copy()
            if "guild" in c:a0=a0[a0.guild==c["guild"]]
            result["estimands"][key]={"definition":c,"all_support":support(a0,c),
                "fits":{k:fit(a0,c,k) for k in STRATA}}
        panels[metric]=result
        print("EXTERNAL_MEASURED_PARTIAL_PANEL",metric,json.dumps(result,ensure_ascii=False,allow_nan=False),flush=True)
    out={"source":source,"trait_pins":pins,"species_pins":species_pins,"panels":panels,
         "warning":"Post-outcome follow-up to non-estimable four-cell calibration, new two-group estimands. Detailed guild is only a coverage census, not a measured specialization index. No equivalence of unrecorded and absent. SC/AF same row, pollinator exclusion may include apomixis. No geographic or causal effect."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False))
    print("EXTERNAL_MEASURED_PARTIAL_COMPLETE",json.dumps({"metrics":len(panels),"estimands_per_metric":len(ESTIMANDS),"fits_per_estimand":len(STRATA),"source_sha256":BOOK_SHA}),flush=True)
if __name__=="__main__":main()
