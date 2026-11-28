#!/usr/bin/env python3
"""Local-pool-adjusted cross-site visitor-channel concordance, external inputs only."""
import argparse, collections, itertools, json, pathlib, tempfile
import numpy as np
import pandas as pd
from external_euppollnet_bee_loss import BASE,SOURCES,BEE,ORDERS,PLANT,VISITOR,REQ,download

SITE=["Latitude","Longitude"]
BEE_F={x.casefold() for x in BEE}
CHANNELS={"hymenoptera","diptera","lepidoptera","coleoptera"}

def channel(order,family):
    if order=="hymenoptera":
        return "bee" if family in BEE_F else "other_hymenoptera"
    return order if order in CHANNELS else "other_orders"

def prepare(d,meta,eligible_study_year=None):
    assert set(REQ).issubset(d)
    bad=d[REQ].fillna("").apply(lambda x:x.str.contains("\ufffd",regex=False)).any(axis=1)
    assert bad.mean()<=.01
    d=d.loc[~bad,REQ].copy()
    d["year"]=pd.to_numeric(d.Date.str.extract(r"^(\d{4})-",expand=False),errors="coerce")
    meta=meta.copy()
    meta["year"]=pd.to_numeric(meta.Year,errors="coerce")
    meta["broad"]=meta.Taxa_recorded.fillna("").str.lower().map(lambda v:all(k in v for k in ORDERS))
    ok=meta.groupby(["Study_id","year"]).broad.all().reset_index()
    eligible=ok.loc[ok.broad,["Study_id","year"]] if eligible_study_year is None else eligible_study_year[["Study_id","year"]].drop_duplicates()
    d=d.merge(eligible,on=["Study_id","year"],validate="many_to_one")
    for k in SITE:d[k]=pd.to_numeric(d[k],errors="coerce")
    for k in ["Plant_rank","Pollinator_rank","Pollinator_order","Pollinator_family"]:
        d[k]=d[k].str.strip().str.casefold()
    d=d.loc[d.Latitude.between(-90,90)&d.Longitude.between(-180,180)&
       d.Plant_rank.eq("species")&d.Pollinator_rank.eq("species")&
       d[PLANT].notna()&d[VISITOR].notna()&d.Pollinator_order.notna()&d.Pollinator_family.notna()&
       d.Pollinator_order.ne("")&d.Pollinator_family.ne("")&
       pd.to_numeric(d.Interaction,errors="coerce").gt(0)&d.Network_id.notna()].copy()
    d["channel"]=[channel(a,b) for a,b in zip(d.Pollinator_order,d.Pollinator_family)]
    d["family"]=[a+":"+b for a,b in zip(d.Pollinator_order,d.Pollinator_family)]
    return d

def concordance(a,b,labels):
    aa=collections.Counter(labels[x] for x in a)
    bb=collections.Counter(labels[x] for x in b)
    return sum(v*bb[k] for k,v in aa.items())/(len(a)*len(b))

def null(a,b,labels):
    """Exact cross-site weighted draws conditional on different species."""
    aa=collections.Counter();bb=collections.Counter()
    for x,v in a.items():aa[labels[x]]+=v
    for x,v in b.items():bb[labels[x]]+=v
    same=sum(v*b.get(x,0) for x,v in a.items())
    den=sum(a.values())*sum(b.values())-same
    if den<=0:return None
    num=sum(v*bb[k] for k,v in aa.items())-same
    assert 0<=num<=den
    return num/den

def background(degree,partners,weighted):
    # Remove focal plant's edges; no abundance or sampling effort inferred.
    return {s:(v-int(s in partners) if weighted else 1)
            for s,v in degree.items() if v>int(s in partners)}

def collect(d):
    counts=collections.Counter();records=[]
    for (study,year),g in d.groupby(["Study_id","year"],sort=True):
        # Exclude visitor identities with contradictory family/order within this study-year.
        bad=set(g.groupby(VISITOR).family.nunique().loc[lambda s:s>1].index)
        counts["conflicted_visitor_study_years"]+=len(bad)
        g=g.loc[~g[VISITOR].isin(bad)]
        labels={k:dict(zip(g[VISITOR],g[k])) for k in ["channel","family"]}
        sites={}
        for coord,h in g.groupby(SITE,sort=True):
            edges=h[[PLANT,VISITOR]].drop_duplicates()
            if edges[PLANT].nunique()<5 or edges[VISITOR].nunique()<5:continue
            sites[coord]=(edges.groupby(PLANT)[VISITOR].agg(set).to_dict(),
                          edges.groupby(VISITOR)[PLANT].nunique().to_dict())
        counts["eligible_sites"]+=len(sites)
        for (ia,(pa,da)),(ib,(pb,db)) in itertools.combinations(sites.items(),2):
            counts["eligible_site_pairs"]+=1
            for plant in sorted(pa.keys()&pb.keys()):
                counts["shared_plant_site_pairs"]+=1
                a,b=pa[plant],pb[plant]
                ea,eb=a-b,b-a
                if not ea or not eb:continue
                counts["two_sided_taxonomic_turnover"]+=1
                wa,wb=background(da,a,True),background(db,b,True)
                ua,ub=background(da,a,False),background(db,b,False)
                if min(len(wa),len(wb))<2:continue
                rec={"study":study,"year":int(year),"plant":plant,"site_a":ia,"site_b":ib,
                     "min_exclusive":min(len(ea),len(eb)),
                     "taxonomic_jaccard":len(a&b)/len(a|b)}
                for k,ll in labels.items():
                    o=concordance(ea,eb,ll)
                    nw,nu=null(wa,wb,ll),null(ua,ub,ll)
                    if nw is None or nu is None:break
                    rec.update({k+"_observed":o,k+"_degree_null":nw,k+"_uniform_null":nu,
                               k+"_residual":o-nw,k+"_uniform_residual":o-nu})
                else:records.append(rec)
    return pd.DataFrame(records),dict(counts)

def summarize(d):
    if d.empty:return {"status":"no_eligible_comparisons","primary_support_gate":False}
    values=["taxonomic_jaccard"]+[f"{k}_{v}" for k in ("channel","family") for v in ("observed","degree_null","uniform_null","residual","uniform_residual")]
    by=d.groupby(["study","year","plant"])[values].mean().groupby(level=["study","year"]).mean().groupby(level="study").mean()
    n=len(by)
    support={"comparisons":len(d),"species":int(d.plant.nunique()),"studies":n,
             "study_years":int(d[["study","year"]].drop_duplicates().shape[0]),
             "site_pairs":int(d[["study","year","site_a","site_b"]].drop_duplicates().shape[0])}
    gate=n>=10 and support["species"]>=50 and support["site_pairs"]>=30
    r={"support":support,"primary_support_gate":bool(gate),"study_equal":{k:float(v) for k,v in by.mean().items()}}
    if n>=2:
        v=by[["channel_residual","channel_uniform_residual","family_residual"]].to_numpy()
        rng=np.random.default_rng(20261128)
        boot=v[rng.integers(n,size=(9999,n))].mean(axis=1)
        loo=(v.sum(axis=0)-v)/(n-1)
        for i,k in enumerate(["channel_residual","channel_uniform_residual","family_residual"]):
            r[k+"_study_bootstrap_ci"]=np.quantile(boot[:,i],[.025,.975]).tolist()
            r[k+"_leave_study_range"]=[float(loo[:,i].min()),float(loo[:,i].max())]
    r["primary_supported"]=bool(gate and r.get("channel_residual_study_bootstrap_ci",[0])[0]>0
                               and r.get("channel_residual_leave_study_range",[0])[0]>0)
    return r

def self_test():
    labels={"a":"bee","b":"bee","c":"fly","d":"fly"}
    a={"a":2,"b":1,"c":3};b={"a":1,"d":4}
    pairs=[(x,y,v*w) for x,v in a.items() for y,w in b.items() if x!=y]
    exact=sum(w for x,y,w in pairs if labels[x]==labels[y])/sum(w for x,y,w in pairs)
    assert abs(null(a,b,labels)-exact)<1e-12
    assert null({"a":1},{"a":1},labels) is None
    assert concordance({"a"},{"b"},labels)==1
    assert background({"a":1,"b":2,"c":1},{"a","b"},True)=={"b":1,"c":1}
    rec=[]
    for si,coords in enumerate([(45,1),(46,1)]):
        for p in range(5):
            for v in (p,(p+1)%5):
                q=f"v{si}_{v}"
                rec.append({"Study_id":"s","year":2020,"Latitude":coords[0],"Longitude":coords[1],
                            PLANT:f"p{p}",VISITOR:q,"channel":"bee" if v<3 else "diptera",
                            "family":"hymenoptera:apidae" if v<3 else "diptera:syrphidae"})
    d,counts=collect(pd.DataFrame(rec))
    assert len(d)==5 and counts["two_sided_taxonomic_turnover"]==5
    r=summarize(d)
    assert r["support"]["species"]==5 and not r["primary_support_gate"]
    json.dumps(r,allow_nan=False)
    print("FD_EXTERNAL_REPLACEMENT_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true");p.add_argument("--out",type=pathlib.Path)
    a=p.parse_args();self_test()
    if a.self_test:return
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);pins={}
        for f,(sha,blob) in SOURCES.items():pins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
        d=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False)
        assert len(d)==623476
        meta=pd.read_csv(w/"Metadata.csv",dtype=str)
        filtered=prepare(d,meta);records,census=collect(filtered)
        output={"source_pins":pins,"source_rows":len(d),"qualified_rows":len(filtered),
                "qualified_studies":int(filtered.Study_id.nunique()),"census":census,
                "primary":summarize(records),
                "sensitivity_min_two_exclusive":summarize(records.loc[records.min_exclusive.ge(2)] if not records.empty else records),
                "interpretation":"Same species, different sites within study-year; conditional on observed two-sided taxonomic turnover. Coarse visitor channels are not measured functional effectiveness. No island-specific inference, historical change, or independent replication of other EuPPollNet results."}
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(output,indent=2,allow_nan=False))
        print("FD_EXTERNAL_REPLACEMENT_RESULT",json.dumps(output,allow_nan=False),flush=True)
if __name__=="__main__":main()
