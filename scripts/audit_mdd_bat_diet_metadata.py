#!/usr/bin/env python3
"""Metadata-only MDD mapped bat species × original EltonTraits mammal diet audit.
No range maps, island incidence, or ecological regression. Anonymous counts only.
"""
import argparse,collections,hashlib,io,json,pathlib,re,urllib.request,time
import pandas as pd
import numpy as np
MDD="https://zenodo.org/records/6644198/files/mdd_spList_wFamilieswOrders_mapped_6362species.csv?download=1"
MDD_MD5="4a800c367f7f8d767779d2f68df40f93"
DIET="https://figshare.com/ndownloader/files/5631084"
def fetch(url):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"fd research provenance audit"}),timeout=60) as r:
                b=r.read(5_000_001)
            if len(b)>5_000_000:raise ValueError("metadata too large")
            return b
        except Exception:
            if attempt==2:raise
            time.sleep(2**attempt)
def norm(v):
    if not isinstance(v,str):return None
    s=v.strip().replace("_"," ").split()
    if len(s)!=2 or any(not re.fullmatch("[A-Za-z-]+",x) for x in s):return None
    return " ".join(s).casefold()
def col(d,choices):
    keys={re.sub("[^a-z0-9]","",x.lower()):x for x in d.columns}
    for x in choices:
        if x in keys:return keys[x]
    raise ValueError("required source schema field missing: "+",".join(choices))
def table(raw,sep):
    for enc in ("utf-8-sig","latin-1"):
        try:return pd.read_csv(io.BytesIO(raw),sep=sep,encoding=enc,dtype=str,keep_default_na=False)
        except UnicodeDecodeError:pass
    raise ValueError("unknown encoding")
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args();a.out.parent.mkdir(parents=True,exist_ok=True)
    out={"status":"metadata_only","no_spatial_or_ecological_model":True}
    try:
        b=fetch(MDD)
        if hashlib.md5(b).hexdigest()!=MDD_MD5:raise ValueError("mapped species file hash mismatch")
        d=table(b,",");out["mapped_source"]={"rows":len(d),"sha256":hashlib.sha256(b).hexdigest(),"md5":MDD_MD5,"columns":list(d.columns)}
        assert len(d)==6362
        order=col(d,["order","ordername","mddorder"]);family=col(d,["family","familyname","mddfamily"])
        try:
            species=col(d,["scientificname","scientific","sciname","binomial","msw3name","mddname","speciesname"])
            names=d[species].map(norm)
            if names.notna().sum()<len(d)*.95:raise ValueError("not binomial")
        except ValueError:
            gen=col(d,["genus","genusname","mddgenus"]);sp=col(d,["species","specificepithet","epithet","mddspecies"])
            names=(d[gen]+" "+d[sp]).map(norm)
        bats=d[d[order].str.lower().eq("chiroptera")].copy()
        bats["name"]=names.loc[bats.index]
        if bats.name.isna().any() or bats.name.duplicated().any() or not len(bats):raise ValueError("invalid mapped bat names")
        out["mapped_bat_taxa"]=len(bats)
    except Exception as e:
        out["mapped_source_status"]="unavailable_or_schema_failed";out["mapped_error_type"]=type(e).__name__;out["mapped_error"]=str(e)[:180]
        bats=None
    try:
        b=fetch(DIET);d=table(b,"\t")
        n=col(d,["scientific"]);nect=col(d,["dietnect"]);fruit=col(d,["dietfruit"]);fam=col(d,["mswfamilylatin","family"])
        d["name"]=d[n].map(norm);d["nectar"]=pd.to_numeric(d[nect],errors="coerce");d["fruit"]=pd.to_numeric(d[fruit],errors="coerce")
        assert len(d)==5400
        valid=d.name.notna()&d.nectar.between(0,100)&d.fruit.between(0,100)
        z=d[valid].copy()
        conflicted=set(z.loc[z.name.duplicated(False),"name"])
        z=z[~z.name.isin(conflicted)].set_index("name")
        out["elton_mammals"]={"rows":len(d),"sha256":hashlib.sha256(b).hexdigest(),"valid_unique_diet_names":len(z),"duplicate_names_excluded":len(conflicted),"diet_provenance_fields":[x for x in d.columns if re.search("source|certainty|speclevel",x,re.I) and re.search("diet",x,re.I)]}
        if bats is not None:
            q=bats.join(z[["nectar","fruit",fam]],on="name",validate="one_to_one")
            eligible=q.nectar.notna()
            out["exact_name_match"]={"mapped_bats":len(q),"matched":int(eligible.sum()),"unmatched":int((~eligible).sum()),"unmatched_share":float((~eligible).mean()),"diet_classes_matched":{"nectar_ge50":int(q.nectar.ge(50).sum()),"nectar_partial":int((q.nectar.gt(0)&q.nectar.lt(50)).sum()),"nectar_zero":int(q.nectar.eq(0).sum()),"fruit_positive":int(q.fruit.gt(0).sum())},"within_5pct_global_list":bool((~eligible).mean()<=.05)}
            out["status"]="taxonomy_eligibility_only"
        else:out["status"]="partial_diet_metadata_only"
    except Exception as e:
        out["diet_source_status"]="unavailable_or_schema_failed";out["diet_error_type"]=type(e).__name__;out["diet_error"]=str(e)[:180]
        out["status"]="partial_or_unavailable_metadata"
    a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print("BAT_METADATA_AUDIT",json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()
