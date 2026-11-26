#!/usr/bin/env python3
"""Metadata-defined official native roster overlap with public Wang/NHM traits."""
import argparse,collections,hashlib,io,json,pathlib,subprocess,tempfile
import numpy as np
import pandas as pd
from tropical_pollinator_source_preflight import get,blobsha
from external_wang_nhm_trait_bridge import ROOT,FILES,NHM_ROOT,NHM,strict_nhm,wang_tables
from external_gift_native_flora_roster import cohort,SCOPES,self_test as roster_test

def source_tables():
    b=get(NHM_ROOT+"combined.csv",maxbytes=25_000_000)
    assert blobsha(b)==NHM["combined.csv"]
    nhm,na=strict_nhm(b)
    pins={"nhm_combined":{"git_blob":blobsha(b),"sha256":hashlib.sha256(b).hexdigest()}}
    with tempfile.TemporaryDirectory() as td:
        p=pathlib.Path(td)
        for f,expected in FILES.items():
            raw=get(ROOT+f+"?download=1",maxbytes=15_000_000);assert hashlib.md5(raw).hexdigest()==expected
            (p/f).write_bytes(raw)
            pins[f]={"md5":expected,"sha256":hashlib.sha256(raw).hexdigest()}
        expr='a<-commandArgs(TRUE);e<-new.env(parent=emptyenv());nm<-load(a[1],envir=e);stopifnot(length(nm)==1,is.data.frame(e[[nm[1]]]));write.csv(e[[nm[1]]],a[2],row.names=FALSE,na="")'
        subprocess.run(["Rscript","--vanilla","-e",expr,str(p/next(iter(FILES))),str(p/"wang.csv")],check=True)
        wang,wa=wang_tables(pd.read_csv(p/"wang.csv",dtype=str,keep_default_na=False))
    fc=set(nhm["flower colour"]);cc=set(nhm["corolla colour"]);fs=set(nhm["flower shape"]);tube=set(nhm["corolla_tube_range_cm"])
    base={"flower_color":fc,"flower_form":fs,"corolla_color":cc,"tube_range":tube,"flower_color_form":fc&fs,"corolla_color_tube":cc&tube}
    for scope,sym in wang.items():
        sy=set(sym)
        for k,v in {"symmetry":sy,"flower_color_symmetry":fc&sy,"flower_form_symmetry":fs&sy,"tube_symmetry":tube&sy,"flower_color_form_symmetry":fc&fs&sy,"corolla_color_tube_symmetry":cc&tube&sy}.items():
            base[scope+"__"+k]=v
    return base,{"source_pins":pins,"nhm":na,"wang":wa,"global_carrier_counts":{k:len(v) for k,v in base.items()}}
def census(plan,data,traits):
    panels={}
    for scope in SCOPES:
        for mode in ["inclusive","native_name_disjoint"]:
            rows=[]
            for p in plan:
                i=p["target"];m=p["donor"]
                if i not in data or m not in data:continue
                t=data[i][scope];s=data[m][scope]
                if mode=="native_name_disjoint":s=s-data[i][SCOPES[0]]
                rows.append((p,t,s))
            key=scope+"__"+mode
            panels[key]={}
            for reg in ("north","tropical","south"):
                rr=[(p,t,s) for p,t,s in rows if p["region"]==reg]
                info={"downloaded_pairs":len(rr),"distinct_donors":len({p["donor"] for p,t,s in rr}),"blocks":len({p["block"] for p,t,s in rr}),"endpoints":{}}
                for label,v in traits.items():
                    cases=[(p,t&v,s&v) for p,t,s in rr]
                    k={}
                    for n in (5,10,30):
                        z=[(p,t,s) for p,t,s in cases if min(len(t),len(s))>=n]
                        k[str(n)]={"eligible_pairs":len(z),"blocks":len({p["block"] for p,t,s in z}),"donors":len({p["donor"] for p,t,s in z}),"distinct_target_carriers":len(set().union(*(t for p,t,s in z))),"distinct_donor_carriers":len(set().union(*(s for p,t,s in z)))}
                    info["endpoints"][label]={"support":k,"target_carrier_median":float(np.median([len(t) for p,t,s in cases])) if cases else None,"donor_carrier_median":float(np.median([len(s) for p,t,s in cases])) if cases else None}
                panels[key][reg]=info
    return panels
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args()
    roster_test()
    plan,data,rosters=cohort()
    assert data,"no checklists retrieved"
    traits,source=source_tables()
    panels=census(plan,data,traits)
    out={"rosters":rosters,"traits":source,"panels":panels,"warning":"Publisher-only coverage pilot; GIFT entities may be archipelagos, nearest mainland centroid not confirmed source, no isolation measurement. Strict native flag cohort and missing-question-flag sensitivity. Name disjoint removes all eligible literal native target binomials in loose cohort, not all possible island species. Exact binomial matching conservative, no synonym resolution. Wang is publisher-coded, unannotated not direct evidence. Source text dependence unknown. No biological effect or raw plant rows exported."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("FD_EXTERNAL_NATIVE_TAXONOMY",json.dumps({k:v for k,v in rosters.items() if k!="checklist_sources"},allow_nan=False),flush=True)
    print("FD_EXTERNAL_NATIVE_GLOBAL_CARRIERS",json.dumps(source["global_carrier_counts"],allow_nan=False),flush=True)
    for k,v in panels.items():print("FD_EXTERNAL_NATIVE_COVERAGE_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
    print("FD_EXTERNAL_NATIVE_COVERAGE_COMPLETE",json.dumps({"panels":len(panels),"planned_pairs":len(plan),"successful_lists":rosters["successful_lists"]}),flush=True)
if __name__=="__main__":main()
