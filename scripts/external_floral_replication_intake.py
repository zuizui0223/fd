#!/usr/bin/env python3
"""Public publisher-only floral source readiness and exact-name overlap audit.

NHM's public GitHub release is distinct from its unavailable portal API.
AusTraits 5.0.0 is a version pin, not a latest-version claim.
No project biological input, trait imputation, or geographic effect estimate.
"""
import argparse, hashlib, io, json, pathlib, re, zipfile
from collections import Counter
import pandas as pd
from tropical_pollinator_source_preflight import get,blobsha

NHM_HEAD="4843c22bb2fc91018d9828dcddebf067c83b1097"
NHM_ROOT="https://raw.githubusercontent.com/NaturalHistoryMuseum/plant-traits/"+NHM_HEAD+"/output/raw/ocr-traits/"
NHM={"combined.csv":"e1cfd2d5c518c2061ec91d1a9e8b64fc7ca69389",
     "efloras.flora_of_china.csv":"ec87a99eee2cc1bf9a8868e1d5425654e7974897",
     "efloras.flora_of_north_america.csv":"3710b08d7b041205a6cb5e7b811c08de317e006d",
     "ecoflora.csv":"dbdc6b4469583743b8b35e26912d21d05b550f32"}
AUS_URL="https://zenodo.org/records/10156222/files/austraits-5.0.0.zip?download=1"
AUS_MD5="a4301cd000153533c0391abd9f60b889"
AUS_TRAITS=("flower_colour","perianth_colour","flower_perianth_symmetry","flower_perianth_fusion","flower_length","flower_diameter")
NHM_TRAITS=("flower colour","corolla colour","flower symmetry","flower shape","flower architecture","petal fusion","corolla tube min. length [cm]","corolla tube max. length [cm]")
def sha(b):return hashlib.sha256(b).hexdigest()
def norm(s):return re.sub(r"[^a-z0-9]","",str(s).casefold())
def valid(s):
    return isinstance(s,str) and s.strip().casefold() not in {"","na","nan","none","null","unknown","not recorded"}
def name(s):
    # Exact two-token scientific names only; do not truncate authors,
    # hybrids, infra-specific names, or unidentified taxa into species.
    if not valid(s):return None
    s=" ".join(s.strip().split())
    return s.casefold() if re.fullmatch(r"[A-Z][a-z-]+ [a-z][a-z-]+",s) else None
def column(cols,aliases):
    match=[c for c in cols if norm(c) in {norm(a) for a in aliases}]
    return match[0] if len(match)==1 else None
def consensus(records):
    z={}
    for i,v in records:
        if i is not None and valid(v):z.setdefault(i,set()).add(v.strip().casefold())
    return {i:next(iter(v)) for i,v in z.items() if len(v)==1},sum(len(v)>1 for v in z.values())
def count_states(t):
    n=Counter(t.values())
    return {"distinct_states":len(n),"top_values":dict(n.most_common(15))}
def overlap(t):
    ks=list(t)
    return {a+"__"+b:len(t[a].keys()&t[b].keys()) for i,a in enumerate(ks) for b in ks[i+1:]}
def read_nhm(raw):
    d=pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False)
    cols=list(d.columns)
    idcol=column(cols,["wfo_id","wfoID","wfo identifier","taxonID","taxon_id","wfo"])
    namecol=column(cols,["species","species_name","scientific_name","scientificName","accepted_name","acceptedScientificName","taxon_name","taxon","name","canonical_name","binomial","wfo_name"])
    keys=d[namecol].map(name) if namecol else pd.Series([None]*len(d),index=d.index)
    result={"rows":len(d),"columns":cols,"id_column":idcol,"name_column":namecol,"exact_binomial_rows":int(keys.notna().sum()),"traits":{}}
    tables={}
    for tr in NHM_TRAITS:
        c=column(cols,[tr])
        if c is None:result["traits"][tr]={"status":"column_not_resolved"};continue
        vals=d[c]
        m=vals.map(valid)
        t,conf=consensus(zip(keys,vals));tables[tr]=t
        result["traits"][tr]={"status":"retrieved","nonmissing_rows":int(m.sum()),"distinct_exact_name_consensus":len(t),"conflicting_names":conf,**count_states(t)}
    result["exact_name_joint_coverage"]=overlap(tables)
    return result,tables

def read_aus(raw):
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        files=z.namelist()
        # The versioned CSV archive includes a traits table. No guessed columns.
        choices=[x for x in files if pathlib.PurePosixPath(x).name.casefold()=="traits.csv"]
        if len(choices)!=1:return {"status":"traits_member_not_resolved","members":files[:80]},{}
        d=pd.read_csv(io.BytesIO(z.read(choices[0])),dtype=str,keep_default_na=False,low_memory=False)
    cols=list(d.columns)
    tc=column(cols,["trait_name"]);vc=column(cols,["value"]);nc=column(cols,["taxon_name"]);dc=column(cols,["dataset_id","dataset"])
    info={"status":"retrieved","member":choices[0],"rows":len(d),"columns":cols,"trait_column":tc,"value_column":vc,"name_column":nc,"dataset_column":dc,"traits":{}}
    if not tc or not vc or not nc:return info,{}
    tab={}
    for tr in AUS_TRAITS:
        x=d[d[tc].eq(tr)]
        keys=x[nc].map(name)
        t,conf=consensus(zip(keys,x[vc]));tab[tr]=t
        info["traits"][tr]={"raw_rows":len(x),"exact_binomial_rows":int(keys.notna().sum()),"exact_name_consensus":len(t),"conflicting_names":conf,"source_dataset_counts":dict(Counter(x[dc]).most_common(15)) if dc else {},**count_states(t)}
    info["exact_name_joint_coverage"]=overlap(tab)
    return info,tab

def cross(a,b):
    out={}
    for x,y in [("flower colour","flower_colour"),("corolla colour","perianth_colour"),("flower symmetry","flower_perianth_symmetry"),("petal fusion","flower_perianth_fusion")]:
        aa=a.get(x,{});bb=b.get(y,{});k=aa.keys()&bb.keys()
        # Exact-string identity is a schema diagnostic, NOT harmonized agreement.
        out[x+"__"+y]={"shared_exact_names":len(k),"identical_raw_casefolded_strings":sum(aa[i]==bb[i] for i in k),"candidate_mapping":dict(Counter(aa[i]+" | "+bb[i] for i in k).most_common(20))}
    return out

def self_test():
    assert name("Plant species")=="plant species" and name("Plant species L.") is None
    assert consensus([("a","red"),("a","Red"),("b","white"),("b","red")])==({"a":"red"},1)
    raw=b"species,flower colour,flower shape\nPlant species,red,tubular\nPlant other,white,open\n"
    a,t=read_nhm(raw);assert a["exact_name_joint_coverage"]["flower colour__flower shape"]==2
    print("FD_EXTERNAL_TRAIT_SOURCE_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true");p.add_argument("--out",type=pathlib.Path);a=p.parse_args()
    if a.self_test:self_test();return
    out={"nhm":{},"austraits":{},"cross_source":{}}
    combined={}
    for f,expected in NHM.items():
        try:
            b=get(NHM_ROOT+f,maxbytes=35_000_000);assert blobsha(b)==expected
            info,t=read_nhm(b)
            out["nhm"][f]={"status":"retrieved","url":NHM_ROOT+f,"git_blob":expected,"sha256":sha(b),"bytes":len(b),**info}
            if f=="combined.csv":combined=t
        except Exception as e:out["nhm"][f]={"status":"unavailable_or_invalid","error":type(e).__name__+": "+str(e)[:200]}
        print("FD_NHM_SOURCE_AUDIT",f,json.dumps(out["nhm"][f],allow_nan=False),flush=True)
    try:
        b=get(AUS_URL,maxbytes=45_000_000);assert hashlib.md5(b).hexdigest()==AUS_MD5
        info,t=read_aus(b)
        out["austraits"]={"url":AUS_URL,"version":"5.0.0","doi":"10.5281/zenodo.10156222","sha256":sha(b),"md5":AUS_MD5,**info}
        out["cross_source"]=cross(combined,t)
    except Exception as e:out["austraits"]={"status":"unavailable_or_invalid","error":type(e).__name__+": "+str(e)[:200]}
    print("FD_AUSTRAITS_SOURCE_AUDIT",json.dumps(out["austraits"],allow_nan=False),flush=True)
    print("FD_EXTERNAL_NAME_OVERLAP",json.dumps(out["cross_source"],allow_nan=False),flush=True)
    out["warning"]="Public source bytes and schema coverage only. Exact names not taxonomic harmonization. Raw label agreement not validated measurement agreement. NHM NLP may share underlying flora texts with AusTraits or other sources. Perianth fusion not equivalent to flower form, and tube length not automatically depth or accessibility. No ecological result."
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    assert any(x.get("status")=="retrieved" for x in out["nhm"].values()) or out["austraits"].get("status")=="retrieved","all source acquisition failed"
if __name__=="__main__":main()
