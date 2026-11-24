#!/usr/bin/env python3
"""Exploratory independent bird matrix HWI versus realised within-archipelago occupancy.
Not colonisation probabilities or a cross-guild isolation test."""
import argparse,csv,hashlib,io,json,math,pathlib,urllib.parse
from collections import defaultdict
import numpy as np
from scipy.stats import spearmanr
from tropical_pollinator_source_preflight import BASE,get,blobsha
from tropical_pollinator_source_intake import MATRICES,matrix,get_diet,norm_name,write_csv
URL="https://raw.githubusercontent.com/stephchia/bird-morphology-tda/757d16598fa00ce0e5696a22a54807ad2d299e01/data/raw/AVONET_birdtree.csv"
BLOB="ee00cd347179844041f490e6eb3e3f902cd0ed94"
def key(x):return "".join(c for c in x.casefold() if c.isalnum())
def load_hwi():
    raw=get(URL,maxbytes=5_000_000)
    assert blobsha(raw)==BLOB
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    cols={key(c):c for c in reader.fieldnames}
    for c in ("species3","handwingindex","family3"):
        if c not in cols:raise ValueError("AVONET schema missing "+c+": "+str(reader.fieldnames))
    hwi={}
    for r in reader:
        n=norm_name(r[cols["species3"]])
        if not n:continue
        try:h=float(r[cols["handwingindex"]])
        except (TypeError,ValueError):continue
        if not math.isfinite(h) or h<0 or h>100:continue
        fam=r[cols["family3"]]
        if n in hwi and hwi[n]!=(h,fam):raise ValueError("conflicting AVONET species")
        hwi[n]=(h,fam)
    assert len(hwi)>=9000
    return hwi,hashlib.sha256(raw).hexdigest()

def correlation(rows):
    if len(rows)<5 or len({r["hwi"] for r in rows})<2 or len({r["occupancy"] for r in rows})<2:return None
    return float(spearmanr([r["hwi"] for r in rows],[r["occupancy"] for r in rows]).statistic)
def slope(rows):
    # Equal archipelago x family weight; within-cell demeaned HWI and occupancy.
    cells=defaultdict(list)
    for r in rows:cells[(r["group"],r["family"])].append(r)
    x=[];y=[]
    for v in cells.values():
        if len(v)<3:continue
        h=np.array([r["hwi"] for r in v]);o=np.array([r["occupancy"] for r in v])
        if np.ptp(h)==0 or np.ptp(o)==0:continue
        x.extend((h-h.mean())/math.sqrt(len(v)));y.extend((o-o.mean())/math.sqrt(len(v)))
    if not x:return None
    return float(np.dot(x,y)/np.dot(x,x))*10
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    hwi,sha=load_hwi();diet,dmeta=get_diet();rows=[];coverage={}
    for group,(path,expected) in MATRICES.items():
        raw=get(BASE+urllib.parse.quote(path));assert blobsha(raw)==expected
        islands,species,area,footer=matrix(raw)
        total=matched=0
        for s,v in species.items():
            if not sum(v):continue
            total+=1
            n=norm_name(s)
            if n not in hwi or n not in diet:continue
            matched+=1
            h,fam=hwi[n]
            rows.append({"group":group,"species":n,"family":fam,"hwi":h,"diet_nect":diet[n],"occupied_islands":sum(v),"island_count":len(islands),"occupancy":sum(v)/len(islands)})
        coverage[group]={"source_species":total,"joint_hwi_diet_species":matched,"islands":len(islands)}
    write_csv(a.out/"bird_hwi_occupancy.csv",rows)
    trop=[r for r in rows if r["group"]!="Ryukyu sensitivity"]
    out={}
    for label,cond in [("all_birds",lambda x:True),("nectar_gt0",lambda x:x["diet_nect"]>0),("nectar_ge20",lambda x:x["diet_nect"]>=20),("nectar_ge50",lambda x:x["diet_nect"]>=50),("nectar_zero",lambda x:x["diet_nect"]==0)]:
        sel=[r for r in trop if cond(r)]
        groups={g:{"n":sum(r["group"]==g for r in sel),"spearman":correlation([r for r in sel if r["group"]==g])} for g in sorted(coverage) if g!="Ryukyu sensitivity"}
        supported=[x["spearman"] for x in groups.values() if x["spearman"] is not None]
        drop={g:slope([r for r in sel if r["group"]!=g]) for g in groups}
        out[label]={"archipelago_species_rows":len(sel),"distinct_species":len({r["species"] for r in sel}),"supported_groups":len(supported),"positive_groups":sum(x>0 for x in supported),"median_group_rho":float(np.median(supported)) if supported else None,"family_group_conditioned_occupancy_change_per_10_hwi":slope(sel),"leave_one_group_out":drop,"groups":groups}
    result={"avonet_git_blob":BLOB,"avonet_sha256":sha,"avonet_names":len(hwi),"diet_sha256":dmeta["sha256"],"coverage":coverage,"results":out,"warning":"Realised occupancy among species already recorded in archipelago, not actual dispersal, colonisation, isolation or full guild substitution. No inferential CI; grouped slopes descriptive. HWI only a morphological proxy; seabirds included in all-bird sensitivity."}
    (a.out/"hwi_summary.json").write_text(json.dumps(result,indent=2))
    print("HWI_OCCUPANCY",json.dumps(result),flush=True)
if __name__=="__main__":main()
