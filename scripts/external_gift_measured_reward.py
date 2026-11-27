#!/usr/bin/env python3
"""Measured floral rewards in metadata-defined public native floras.
Use verified raw sugar and publisher log10 pollen; never predict reward from morphology.
Same 288-region discovery roster, exploratory two-endpoint family.
"""
import argparse,collections,hashlib,json,pathlib
import numpy as np
import pandas as pd
import external_gift_native_flora_roster as roster
from external_floral_reward_source import API,PINS,RECORD,verified_file,csv_table,resolve_names
from tropical_pollinator_source_preflight import get
from external_gift_length_holdout import test
PLAN_SHA="8945ce01bd1b5a9294fe18ee570ec6ed9d207e93363c70ec3732ffe013f3c935"
SPECS=(("loose","disjoint",5),("loose","inclusive",5),("loose","disjoint",10),("strict","disjoint",5))
PRIMARY="loose__disjoint__n5"
def table(d,col):
    assert col in d.columns,"expected publisher value column missing"
    c,keys,diagnostic=resolve_names(d)
    if keys is None:raise ValueError("unresolved publisher names")
    vals=pd.to_numeric(d[col],errors="coerce")
    by=collections.defaultdict(set)
    for k,v in zip(keys,vals):
        if k is not None and np.isfinite(v):by[k].add(float(v))
    # No averaging of conflicting observations or imputation.
    out={k:(0,next(iter(v))) for k,v in by.items() if len(v)==1}
    return out,{"rows":len(d),"name_column":c,"value_column":col,"name_audit":diagnostic,"finite_log_values":int(np.isfinite(vals).sum()),"distinct_names":len(by),"conflicts":sum(len(v)>1 for v in by.values()),"usable_species":len(out)}
def nectar_tables(d):
    col="Nectar.microg.per.day"
    assert {col,"Species","Genus","Source"}<=set(d.columns),"pinned nectar schema changed"
    c,keys,diag=resolve_names(d)
    assert keys is not None
    vals=pd.to_numeric(d[col],errors="coerce")
    assert np.isfinite(vals).all() and (vals>=0).all(),"invalid raw sugar production"
    by=collections.defaultdict(lambda:collections.defaultdict(list))
    for k,v,ref in zip(keys,vals,d["Source"]):
        if k is None:continue
        ref=" ".join(str(ref).strip().split()).casefold()
        if ref in {"","na","nan","none","null","unknown","not recorded"}:ref=None
        by[k][ref].append(float(v))
    def mean(xs):return float(np.mean(xs))
    def balanced(z):return mean([mean(x) for x in z.values()])
    def pooled(z):return mean([v for xs in z.values() for v in xs])
    def med(z):return float(np.median([np.median(xs) for xs in z.values()]))
    known={k:{ref:xs for ref,xs in z.items() if ref is not None} for k,z in by.items()}
    known={k:z for k,z in known.items() if z}
    signature={k:tuple(sorted(ref if ref is not None else "<unreported>" for ref in z)) for k,z in by.items()}
    si={z:i for i,z in enumerate(sorted(set(signature.values())))}
    primary={k:(0,float(np.log10(1+balanced(z)))) for k,z in by.items()}
    out={
        "nectar":primary,
        "nectar_source_matched":{k:(si[signature[k]],v[1]) for k,v in primary.items()},
        "nectar_pooled_mean":{k:(0,float(np.log10(1+pooled(z)))) for k,z in by.items()},
        "nectar_source_median":{k:(0,float(np.log10(1+med(z)))) for k,z in by.items()},
        "nectar_known_source":{k:(0,float(np.log10(1+balanced(z)))) for k,z in known.items()},
        "nectar_strict_positive":{k:(0,float(np.log10(balanced(z)))) for k,z in by.items() if min(v for xs in z.values() for v in xs)>0},
        "nectar_detected":{k:(0,float(any(v>0 for xs in z.values() for v in xs))) for k,z in by.items()},
    }
    source_counts=collections.Counter(ref if ref is not None else "<unreported>" for z in by.values() for ref in z)
    audit={"rows":len(d),"name_column":c,"name_audit":diag,"value_column":col,"unit":"micrograms sugar per flower per day","zero_rows":int(vals.eq(0).sum()),"species":len(by),"species_with_multiple_rows":sum(sum(map(len,z.values()))>1 for z in by.values()),"species_with_multiple_known_sources":sum(len(z)>1 for z in known.values()),"species_with_unreported_source":sum(None in z for z in by.values()),"species_all_zero":sum(max(v for xs in z.values() for v in xs)==0 for z in by.values()),"species_mixed_zero_positive":sum(min(v for xs in z.values() for v in xs)==0<max(v for xs in z.values() for v in xs) for z in by.values()),"source_signature_count":len(si),"source_species_counts":dict(source_counts),"variants":{k:len(t) for k,t in out.items()},"aggregation":"primary mean raw sugar within each reported source and then equal source mean (unknown source is one distinct stratum); log10(1+sugar / 1 microgram per flower per day). Repeated values retained, no assumed independent records. Pooled row mean matches paper's stated averaging rule; median and known-source sensitivities. Strict-positive requires all recorded values positive; detected means at least one positive measurement, not obligate reward status."}
    return out,audit

def sources():
    meta=json.loads(get(API,maxbytes=3_000_000));assert int(meta["id"])==RECORD
    f={x["key"]:x for x in meta["files"]}
    out={};audit={}
    for endpoint,file in (("nectar","Nectar_dataset.csv"),("pollen","Pollen_dataset_no_duplicates.csv")):
        try:
            raw,receipt=verified_file(f[file]);d=csv_table(raw,".csv")
            if endpoint=="nectar":
                tabs,z=nectar_tables(d);out.update(tabs)
            else:
                t,z=table(d,"log10Pollen");out[endpoint]=t
            audit[endpoint]={"file":file,"receipt":receipt,**z}
        except Exception as e:
            audit[endpoint]={"file":file,"status":"unavailable","error":str(e)[:400]}
    print("FD_GIFT_REWARD_SOURCE",json.dumps(audit,allow_nan=False),flush=True)
    return out,audit
def selection(plan,data,tab,native,mode,minimum):
    sc=roster.SCOPES[0 if native=="loose" else 1];rows=[]
    for z in plan:
        i=z["target"];m=z["donor"]
        if i not in data or m not in data:continue
        t=data[i][sc]&tab.keys();s=data[m][sc]&tab.keys();d=s-data[i][roster.SCOPES[0]]
        if min(len(t),len(d))<minimum:continue
        rows.append(z|{"t":t,"s":d if mode=="disjoint" else s})
    return rows
def self_test():
    d=pd.DataFrame({"Species":["Abelia_biflora","Abelia biflora","Abelia other","Abelia third","Abelia third"],"Genus":["Abelia"]*5,"log10Pollen":["1","1","-1","2","3"]})
    t,a=table(d,"log10Pollen");assert len(t)==2 and a["conflicts"]==1 and t["abelia other"]==(0,-1.)
    d2=pd.DataFrame({"Species":["Abelia_biflora"]*3+["Abelia_other"],"Genus":["Abelia"]*4,"Source":["A","A","B",""],"Nectar.microg.per.day":["0","2","9","0"]})
    z,q=nectar_tables(d2)
    assert np.isclose(z["nectar"]["abelia biflora"][1],np.log10(6))
    assert np.isclose(z["nectar_pooled_mean"]["abelia biflora"][1],np.log10(1+11/3))
    assert z["nectar"]["abelia other"][1]==0 and "abelia biflora" not in z["nectar_strict_positive"]
    assert q["species_mixed_zero_positive"]==1 and q["species_all_zero"]==1
    print("FD_GIFT_REWARD_SELF_TEST_PASS",flush=True)
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--nectar-column",required=True);p.add_argument("--self-test",action="store_true");a=p.parse_args();self_test()
    if a.self_test:return
    assert a.nectar_column=="auto","fixed verified raw nectar column, no alternate column guessing"
    tabs,audit=sources()
    assert {"nectar","pollen"}<=set(tabs),"both pinned measured reward endpoints required"
    roster.LIMIT=288
    plan,data,prov=roster.cohort()
    assert len(plan)==288 and prov["plan_sha256"]==PLAN_SHA,"geography plan drift"
    print("FD_GIFT_REWARD_ROSTER",json.dumps({k:v for k,v in prov.items() if k!="checklist_sources"},allow_nan=False),flush=True)
    panels={}
    for endpoint,tab in tabs.items():
        for native,mode,n in (SPECS if endpoint in ("nectar","pollen") else (SPECS[0],)):
            key=endpoint+"__"+"__".join([native,mode,"n"+str(n)])
            rows=selection(plan,data,tab,native,mode,n)
            # test() fixes genus x symmetry; all dummy symmetry states zero
            # so this is a whole-species genus-only reward permutation.
            panels[key]=test(rows,tab,False,999 if endpoint in ("nectar","pollen") and key.endswith(PRIMARY) else 199)
            print("FD_GIFT_REWARD_PROGRESS",key,json.dumps({"eligible_pairs":len(rows)},allow_nan=False),flush=True)
    # Retain both endpoints in multiplicity even if unavailable.
    tests=[]
    for endpoint in ("nectar","pollen"):
        z=panels.get(endpoint+"__"+PRIMARY,{}).get("group_weightings",{})
        B=z.get("block",{}).get("tropical_minus_north",{});D=z.get("donor",{}).get("tropical_minus_north",{})
        x=B.get("metrics",{}).get("total",{});y=D.get("metrics",{}).get("total",{})
        gate=bool(B.get("null_gate") and D.get("null_gate"))
        tests.append((x.get("coherent_p",1.) if gate else 1.,endpoint,x,y,gate))
    last=0;verdict={}
    for i,(pv,endpoint,x,y,gate) in enumerate(sorted(tests,key=lambda q:q[0])):
        last=max(last,min(1.,pv*(2-i)))
        def signed(e):
            c=e.get("residual_ci");l=e.get("residual_leave_range")
            if not c or not l:return 0
            return 1 if c[0]>0 and l[0]>0 else -1 if c[1]<0 and l[1]<0 else 0
        verdict[endpoint]={"gate":gate,"holm_p":last,"supported":bool(gate and last<.05 and y.get("coherent_p",1)<.05 and signed(x)!=0 and signed(x)==signed(y))}
    if {"nectar","pollen"}<=set(tabs):
        # Compare endpoints on exactly identical species and geographic carriers.
        # Difference of transformed endpoints is not a physical reward ratio.
        names=tabs["nectar"].keys()&tabs["pollen"].keys()
        matched={endpoint:{s:tabs[endpoint][s] for s in names} for endpoint in ("nectar","pollen")}
        matched["difference"]={s:(0,matched["nectar"][s][1]-matched["pollen"][s][1]) for s in names}
        rows=selection(plan,data,matched["nectar"],"loose","disjoint",5)
        for endpoint,tab in matched.items():
            panels["common_measured__"+endpoint+"__"+PRIMARY]=test(rows,tab,False,199)
        # Pollen primary is the previously observed endpoint; check whether
        # fixed-plan retrieval reproduces it, not a new independent replicate.
        old=panels["pollen__"+PRIMARY].get("group_weightings",{})
        observed=old.get("block",{}).get("tropical_minus_north",{}).get("metrics",{}).get("total",{}).get("observed")
        print("FD_GIFT_REWARD_POLLEN_REPRODUCTION",json.dumps({"previous":-0.0662076285437096,"current":observed,"same_point":bool(observed is not None and np.isclose(observed,-0.0662076285437096,atol=1e-12))}),flush=True)
    for k,z in panels.items():print("FD_GIFT_REWARD_PANEL",k,json.dumps(z,allow_nan=False),flush=True)
    print("FD_GIFT_REWARD_COMPLETE",json.dumps({"panels":len(panels),"verdict":verdict},allow_nan=False),flush=True)
    out={"rosters":prov,"sources":audit,"panels":panels,"verdict":verdict,"warning":"Post-discovery public native flora reward comparison. Nectar endpoint fixed raw micrograms sugar per flower per day aggregated within and equally across reported source strata, log10(1+x) to retain explicit zeros. Pollen publisher log10 grains per flower. Not abundance-weighted productivity, landscape nectar supply, or realized pollination. Missing-question flags allowed primary, strict sensitivity. 288 metadata-selected existing discovery regions; not independent holdout. Two primary regional genus-null tests Holm, both group weights and fixed support gate (>=50 pairs, >=30 blocks, >=20 donors, >=30 movable species and >=10 movable genera). Unavailable endpoint assigned p=1. Same-carrier measured resource contrast descriptive only; difference of transformed endpoints is not a physical sugar/pollen ratio. Source-profile constrained null diagnostic and other nectar variants descriptive; do not rescue failed primary. No imputation or raw taxon export."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
if __name__=="__main__":main()
