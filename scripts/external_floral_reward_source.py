#!/usr/bin/env python3
"""Public measured reward source intake for FD, not a pollinator-identity inference.
Keasar & Bodner (2025), doi:10.1002/ajb2.70125; fixed Zenodo record 17157828.
Only data files, never execute the publisher's code. Export anonymous aggregate audit.
"""
import argparse,collections,hashlib,io,json,pathlib,re,time,urllib.parse,urllib.error,zipfile
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
    # Files are checksum-verified before parsing. Legacy 0xA0 bytes are
    # compatible with Windows-1252; do not silently replace invalid bytes.
    try:
        raw.decode("utf-8-sig",errors="strict");enc="utf-8-sig"
    except UnicodeDecodeError:
        raw.decode("cp1252",errors="strict");enc="cp1252"
    d=pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False,sep=",",encoding=enc)
    d.attrs["decode_encoding"]=enc
    return d
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
PINS={
    "Flower_trait_lookup.csv":(384273,"a1ed24e2509be4293b90bf95c5d6c786"),
    "Nectar_dataset.csv":(88095,"eeb83133bc0a8205bb59798f1eb3add7"),
    "Pollen_dataset_no_duplicates.csv":(87156,"e65c34116c424d1aee677e96686bba1a"),
}
def normalize_reward_name(raw,genus=None):
    # Publisher names can be underscored or use different capitalization.
    # Only two complete binomial tokens, not authors, hybrids or infraspecific names.
    if not valid(raw):return None
    raw=" ".join(raw.strip().replace("_"," ").split())
    m=re.fullmatch(r"([A-Za-z][A-Za-z-]+) ([A-Za-z][A-Za-z-]+)",raw)
    if not m:return None
    g,e=m.groups()
    if e.casefold() in {"sp","spp","cf","aff","unknown","unidentified"}:return None
    if genus is not None and valid(genus) and g.casefold()!=genus.strip().casefold():return None
    return name(g.capitalize()+" "+e.lower())
def resolve_names(d):
    candidates=[]
    genus_cols=[c for c in d.columns if norm(c)=="genus"]
    genus=d[genus_cols[0]] if len(genus_cols)==1 else [None]*len(d)
    for c in d.columns:
        if norm(c) in {"species","speciesname","scientificname","taxon","taxonname","plant","plantname","plantspecies","binomial","acceptedname","speciesid","sp"}:
            k=pd.Series([normalize_reward_name(x,g) for x,g in zip(d[c],genus)],index=d.index,dtype=object)
            candidates.append((str(c),k))
    if not candidates:return None,None,{}
    counts={c:{"normalized_binomials":int(k.notna().sum()),"literal_binomials":int(d[c].map(name).notna().sum()),"underscore_rows":int(d[c].str.contains("_",regex=False).sum())} for c,k in candidates}
    good=[(c,k) for c,k in candidates if counts[c]["normalized_binomials"]>0]
    if len(good)!=1:
        if good and all(good[0][1].equals(k) for _,k in good[1:]):return good[0][0],good[0][1],counts
        return None,None,counts
    return good[0][0],good[0][1],counts
def verified_file(f):
    key=f["key"];size,expected=PINS[key]
    assert int(f["size"])==size and f["checksum"]=="md5:"+expected,"fixed manifest drift"
    url=f"https://zenodo.org/records/{RECORD}/files/"+urllib.parse.quote(key,safe="")+"?download=1"
    errors=[]
    for attempt in range(3):
        try:
            raw=get(url,maxbytes=CAP)
            if len(raw)!=size:raise ValueError(f"length mismatch {len(raw)} != {size}")
            if hashlib.md5(raw).hexdigest()!=expected:raise ValueError("published MD5 mismatch")
            return raw,{"status":"retrieved","size":len(raw),"sha256":digest(raw),"published_checksum":f["checksum"],"prior_attempt_errors":errors}
        except Exception as e:
            errors.append({"type":type(e).__name__,"detail":str(e)[:300]})
            if isinstance(e,urllib.error.HTTPError) and e.code in {401,403,404}:break
            if attempt<2:time.sleep(2*(attempt+1))
    raise RuntimeError(json.dumps(errors))
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
    assert normalize_reward_name("Abelia_biflora","Abelia")=="abelia biflora"
    assert normalize_reward_name("abelia BIFLORA","Abelia")=="abelia biflora"
    assert normalize_reward_name("Abelia_biflora","Other") is None
    assert normalize_reward_name("Abelia biflora L.") is None
    d2=csv_table(b"Species,Genus,log10Nectar\nAbelia\xa0biflora,Abelia,0.5\n",".csv")
    assert d2.attrs["decode_encoding"]=="cp1252" and resolve_names(d2)[1].iloc[0]=="abelia biflora"
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
    assert {f["key"] for f in files}>=set(PINS)
    frames={}
    for f in files:
        key=f["key"]
        if key not in PINS:continue
        try:
            raw,receipt=verified_file(f)
            out["files"][key]=receipt
            # Separate download verification from parsing; report the actual failure.
            try:
                frames[key]=csv_table(raw,".csv")
                out["files"][key]["parse_status"]="parsed"
                out["files"][key]["decode_encoding"]=frames[key].attrs.get("decode_encoding")
            except Exception as e:
                out["files"][key]["parse_status"]="failed"
                out["files"][key]["parse_error"]=str(e)[:300]
        except Exception as e:
            out["files"][key]={"status":"failed","error_type":type(e).__name__,"detail":str(e)[:400]}
        print("FD_REWARD_FILE",json.dumps({"file":key,"receipt":out["files"][key]},ensure_ascii=False),flush=True)
    _,pinned,nhm,wang=source_tables(include_values=True)
    sy=wang["unannotated_without_four_families"];tube=nhm["corolla_tube_range_cm"]
    common=sy.keys()&tube.keys()
    carriers={"NHM_tube":set(tube),"Wang_conservative_symmetry":set(sy),"common_tube_symmetry":common,"radial_tube":{x for x in common if sy[x]=="radial"},"bilateral_tube":{x for x in common if sy[x]=="bilateral"}}
    name_sets={}
    for key,d in frames.items():
        out["tables"][key]=summary(d,carriers)
        c,k,_=resolve_names(d)
        if k is not None:name_sets[key]=set(k.dropna())
        print("FD_REWARD_TABLE",json.dumps({"table":key,"audit":out["tables"][key]},ensure_ascii=False,allow_nan=False),flush=True)
    out["exact_name_pair_overlap"]={a+"__"+b:len(name_sets[a]&name_sets[b]) for i,a in enumerate(name_sets) for b in list(name_sets)[i+1:]}
    print("FD_REWARD_NAME_OVERLAP",json.dumps(out["exact_name_pair_overlap"]),flush=True)
    # Audit the authors\u0027 code as text only for the reward variables and units.
    codefile=next((f for f in files if f["key"]=="Flower_shape_and_nectar_pollen_rewards_V5.R"),None)
    if codefile is not None:
        try:
            assert codefile["checksum"]=="md5:73405584049dbe3ae3413b0adc81bb8e"
            b=get(f"https://zenodo.org/records/{RECORD}/files/"+urllib.parse.quote(codefile["key"],safe="")+"?download=1",maxbytes=100_000)
            assert hashlib.md5(b).hexdigest()=="73405584049dbe3ae3413b0adc81bb8e" and len(b)==13437
            # Extract only lines defining log transformation or explicitly naming units.
            lines=[x.strip() for x in b.decode("cp1252").splitlines() if re.search(r"log10|microgram|milligram|mg/|ug/|sugar.*day",x,re.I)]
            out["publisher_code_units_audit"]={"sha256":digest(b),"lines":lines[:30]}
            print("FD_REWARD_PUBLISHER_UNITS",json.dumps(out["publisher_code_units_audit"]),flush=True)
        except Exception as e:out["publisher_code_units_audit"]={"status":"unavailable","detail":str(e)[:200]}
    out["trait_bridge"]={"source_pins":pinned["source_pins"],"common_carriers":{k:len(v) for k,v in carriers.items()}}
    out["warning"]="Publisher-only readiness; exact binomial names, no taxonomic synonyms or unknown units guessed; tables may contain repeated source measurements, not independent species. No raw names exported. Actual measured nectar and pollen must be harmonized with source units and references before community tests. Reward quantity is not accessibility, pollinator identity or realized pollination."
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False))
    print("FD_REWARD_SOURCE_COMPLETE",json.dumps({"retrieved_files":sum(v["status"]=="retrieved" for v in out["files"].values()),"tables":len(out["tables"])}),flush=True)
    assert set(frames)==set(PINS),"not all fixed publisher data tables parsed"
    assert all(x.get("exact_distinct_names",0)>0 for x in out["tables"].values()),"unresolved binomial names"
if __name__=="__main__":main()
