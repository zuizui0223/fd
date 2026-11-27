#!/usr/bin/env python3
"""Public measured reward source intake for FD, not a pollinator-identity inference.
Keasar & Bodner (2025), doi:10.1002/ajb2.70125; fixed Zenodo record 17157828.
Only data files, never execute the publisher's code. Export anonymous aggregate audit.
"""
import argparse,collections,hashlib,io,json,pathlib,re,urllib.parse,zipfile
import numpy as np
import pandas as pd
from tropical_pollinator_source_preflight import get
from external_floral_replication_intake import name,valid
from external_gift_floral_replication_coverage import source_tables
RECORD=17157828
API=f"https://zenodo.org/api/records/{RECORD}"
CAP=40_000_000
ALLOWED={".csv",".tsv",".xlsx",".xls"}
def digest(raw):return hashlib.sha256(raw).hexdigest()
def csv_table(raw,ext):
    if ext==".tsv":return pd.read_csv(io.BytesIO(raw),sep="\t",dtype=str,keep_default_na=False)
    return pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False,sep=None,engine="python")
def tables(raw,filename):
    ext=pathlib.PurePosixPath(filename).suffix.lower()
    if ext in {".csv",".tsv"}:
        yield filename,csv_table(raw,ext)
    elif ext in {".xlsx",".xls"}:
        f=pd.ExcelFile(io.BytesIO(raw));assert len(f.sheet_names)<=40
        for sheet in f.sheet_names:
            yield filename+"::"+sheet,pd.read_excel(f,sheet_name=sheet,dtype=str,keep_default_na=False)
    elif ext==".zip":
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            assert sum(x.file_size for x in z.infolist())<=CAP
            for f in z.infolist():
                if pathlib.PurePosixPath(f.filename).suffix.lower() in ALLOWED and not f.is_dir():
                    yield from tables(z.read(f),filename+"::"+f.filename)
def norm(s):return re.sub(r"[^a-z0-9]","",str(s).casefold())
def resolve_names(d):
    candidates=[]
    for c in d.columns:
        if norm(c) in {"species","speciesname","scientificname","taxon","taxonname","plant","plantname","plantspecies","binomial","acceptedname","speciesid","sp"}:
            k=d[c].map(name)
            candidates.append((str(c),k))
    if not candidates:return None,None,{}
    counts={c:int(k.notna().sum()) for c,k in candidates}
    good=[(c,k) for c,k in candidates if counts[c]>0]
    if len(good)!=1:
        if good and all(good[0][1].equals(k) for _,k in good[1:]):return good[0][0],good[0][1],counts
        return None,None,counts
    return good[0][0],good[0][1],counts
def summary(d,carriers):
    c,k,counts=resolve_names(d)
    z={"rows":len(d),"columns":[str(x) for x in d.columns],"name_column":c,"candidate_name_column_counts":counts,"numeric_fields":{}}
    if k is not None:
        ids=set(k.dropna());z.update(exact_binomial_rows=int(k.notna().sum()),exact_distinct_names=len(ids),repeated_name_rows=int(k.notna().sum())-len(ids),carrier_overlap={a:len(ids&v) for a,v in carriers.items()})
    for c in d.columns:
        if re.search(r"nectar|sugar|pollen|tube|length|symm|source|reference|unit",str(c),re.I):
            x=d[c];y=pd.to_numeric(x,errors="coerce")
            good=y[np.isfinite(y)]
            z["numeric_fields"][str(c)]={"nonempty":int(x.map(valid).sum()),"numeric":len(good),"positive":int(good.gt(0).sum()),"zero":int(good.eq(0).sum()),"negative":int(good.lt(0).sum()),"quantiles":good.quantile([0,.25,.5,.75,1]).tolist() if len(good) else None,"distinct_nonempty":int(x[x.map(valid)].nunique())}
    return z
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    d=pd.DataFrame({"Species":["Abelia biflora","Abelia biflora","Abelia sp."],"Sugar mg per flower per day":["1","2","0"]})
    v=summary(d,{"tube":{"abelia biflora"}});assert v["exact_distinct_names"]==1 and v["repeated_name_rows"]==1 and v["carrier_overlap"]["tube"]==1
    print("FD_REWARD_SOURCE_SELF_TEST_PASS",flush=True)
    if a.self_test:return
    # Retrieve the fixed record, verify the published file checksums before reading.
    raw=get(API,maxbytes=3_000_000);meta=json.loads(raw)
    assert int(meta["id"])==RECORD
    title=str(meta["metadata"]["title"]);assert "nectar" in title.casefold() and "pollen" in title.casefold()
    desc=meta["metadata"]
    files=meta.get("files",[]);assert isinstance(files,list) and files
    out={"source":{"record":RECORD,"doi":meta.get("doi"),"title":title,"metadata_sha256":digest(raw),"license":desc.get("license"),"access_right":desc.get("access_right"),"files":[{k:f.get(k) for k in ("key","size","checksum")} for f in files]},"files":{},"tables":{}}
    print("FD_REWARD_SOURCE_MANIFEST",json.dumps(out["source"],ensure_ascii=False),flush=True)
    _,pinned,nhm,wang=source_tables(include_values=True)
    sy=wang["unannotated_without_four_families"];tube=nhm["corolla_tube_range_cm"]
    common=sy.keys()&tube.keys()
    carriers={"NHM_tube":set(tube),"Wang_conservative_symmetry":set(sy),"common_tube_symmetry":common,"radial_tube":{x for x in common if sy[x]=="radial"},"bilateral_tube":{x for x in common if sy[x]=="bilateral"}}
    for f in files:
        key=f["key"];ext=pathlib.PurePosixPath(key).suffix.lower()
        if ext not in ALLOWED|{".zip"}:continue
        if int(f["size"])>CAP:continue
        url=f"https://zenodo.org/records/{RECORD}/files/"+urllib.parse.quote(key,safe="")+"?download=1"
        try:
            raw=get(url,maxbytes=CAP);assert len(raw)==int(f["size"])
            alg,h=f["checksum"].split(":",1);assert alg in {"md5","sha256"} and hashlib.new(alg,raw).hexdigest()==h
            out["files"][key]={"status":"retrieved","size":len(raw),"sha256":digest(raw),"published_checksum":f["checksum"]}
            for table,d in tables(raw,key):
                out["tables"][table]=summary(d,carriers)
                print("FD_REWARD_TABLE",json.dumps({"table":table,"audit":out["tables"][table]},ensure_ascii=False,allow_nan=False),flush=True)
        except Exception as e:
            out["files"][key]={"status":"failed","error_type":type(e).__name__,"detail":str(e)[:300]}
    out["trait_bridge"]={"source_pins":pinned["source_pins"],"common_carriers":{k:len(v) for k,v in carriers.items()}}
    out["warning"]="Publisher-only readiness; exact binomial names, no taxonomic synonyms or unknown units guessed; tables may contain repeated source measurements, not independent species. No raw names exported. Actual measured nectar and pollen must be harmonized with source units and references before community tests. Reward quantity is not accessibility, pollinator identity or realized pollination."
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False))
    print("FD_REWARD_SOURCE_COMPLETE",json.dumps({"retrieved_files":sum(v["status"]=="retrieved" for v in out["files"].values()),"tables":len(out["tables"])}),flush=True)
    assert out["tables"],"no usable published data tables"
if __name__=="__main__":main()
