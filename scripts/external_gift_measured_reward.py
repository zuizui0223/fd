#!/usr/bin/env python3
"""Measured floral rewards in metadata-defined public native floras.
Use publisher log-scale values only; never predict reward from morphology.
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
def sources(nectar_col):
    meta=json.loads(get(API,maxbytes=3_000_000));assert int(meta["id"])==RECORD
    f={x["key"]:x for x in meta["files"]}
    out={};audit={}
    for endpoint,file,col in (("nectar","Nectar_dataset.csv",nectar_col),("pollen","Pollen_dataset_no_duplicates.csv","log10Pollen")):
        try:
            raw,receipt=verified_file(f[file]);d=csv_table(raw,".csv")
            if col not in d.columns:raise ValueError("unverified log10 reward column; actual columns: "+str(list(d.columns)))
            t,z=table(d,col);out[endpoint]=t
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
    print("FD_GIFT_REWARD_SELF_TEST_PASS",flush=True)
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--nectar-column",required=True);p.add_argument("--self-test",action="store_true");a=p.parse_args();self_test()
    if a.self_test:return
    assert a.nectar_column in {"log10Nectar","log10Sugar","log10NectarSugar"},"requires verified explicit log10 nectar column"
    tabs,audit=sources(a.nectar_column)
    assert tabs,"no measured reward table retrieved"
    roster.LIMIT=288
    plan,data,prov=roster.cohort()
    assert len(plan)==288 and prov["plan_sha256"]==PLAN_SHA,"geography plan drift"
    print("FD_GIFT_REWARD_ROSTER",json.dumps({k:v for k,v in prov.items() if k!="checklist_sources"},allow_nan=False),flush=True)
    panels={}
    for endpoint,tab in tabs.items():
        for native,mode,n in SPECS:
            key=endpoint+"__"+"__".join([native,mode,"n"+str(n)])
            rows=selection(plan,data,tab,native,mode,n)
            # test() fixes genus x symmetry; all dummy symmetry states zero
            # so this is a whole-species genus-only reward permutation.
            panels[key]=test(rows,tab,False,999 if key.endswith(PRIMARY) else 199)
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
    if set(tabs)=={"nectar","pollen"}:
        # Same species for both traits, contrast of their log changes.
        # This is a dimensionless ratio-of-ratios, NOT a sugar/pollen mass ratio.
        names=tabs["nectar"].keys()&tabs["pollen"].keys()
        common={s:(0,tabs["nectar"][s][1]-tabs["pollen"][s][1]) for s in names}
        rows=selection(plan,data,common,"loose","disjoint",5)
        panels["common_measured_reward_log_contrast__"+PRIMARY]=test(rows,common,False,199)
    for k,z in panels.items():print("FD_GIFT_REWARD_PANEL",k,json.dumps(z,allow_nan=False),flush=True)
    print("FD_GIFT_REWARD_COMPLETE",json.dumps({"panels":len(panels),"verdict":verdict},allow_nan=False),flush=True)
    out={"rosters":prov,"sources":audit,"panels":panels,"verdict":verdict,"warning":"Post-discovery public native flora reward comparison, measured log10 values per-flower reward from publisher, not abundance-weighted productivity, landscape nectar supply, or realized pollination. Missing-question flags allowed primary, strict sensitivity. 288 metadata-selected existing discovery regions; not independent holdout. Two primary regional genus-null tests Holm, both group weights and fixed support gate (>=50 pairs, >=30 blocks, >=20 donors, >=30 movable species and >=10 movable genera). Unavailable endpoint assigned p=1. Same-carrier measured resource contrast descriptive only. No imputation or raw taxon export."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
if __name__=="__main__":main()
