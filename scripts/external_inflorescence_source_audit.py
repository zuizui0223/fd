#!/usr/bin/env python3
"""Read-only audit of the published Wang et al. 2025 Dryad RData.

Official public file, no project database, no source code execution from
archive. Publication totals are targets for reconciliation, not results.
Stop on denied download; no alternate access or credentials.
"""
import argparse,hashlib,json,pathlib,re
from collections import Counter
import pandas as pd
import pyreadr
from tropical_pollinator_source_preflight import get

DOI="10.5061/dryad.73n5tb35w"
META="https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.73n5tb35w"
URL="https://datadryad.org/downloads/file_stream/4049686"
EXPECTED_FILE="Combine_data_uploaded_Final_20250411.Rdata"

def norm(x):return re.sub(r"[^a-z0-9]","",str(x).lower())
def clean(v):return v.astype("string").str.strip().replace({"":pd.NA,"NA":pd.NA,"NaN":pd.NA})
def census(v):
    a=clean(v);c=a.value_counts(dropna=True)
    return {"nonmissing":int(a.notna().sum()),"distinct":int(len(c)),
            "top_values":{str(k):int(x) for k,x in c.head(25).items()}}
def audit_table(df):
    columns={norm(c):c for c in df.columns}
    out={"rows":len(df),"columns":list(map(str,df.columns))}
    wanted=["order","family","genus","genusspecies","inflorescencetype","note","splistinphyloanalysis","polmodbinary","datasourceofpolmod"]
    out["fields"]={n:census(df[columns[n]]) for n in wanted if n in columns and n not in ("genusspecies","genus","family","order")}
    out["taxonomic_counts"]={n:{"distinct":int(clean(df[columns[n]]).nunique()),"nonmissing":int(clean(df[columns[n]]).notna().sum())} for n in wanted[:4] if n in columns}
    if not all(k in columns for k in ("genusspecies","inflorescencetype","polmodbinary")):
        out["status"]="schema_not_confirmed";return out
    s=clean(df[columns["genusspecies"]])
    a=clean(df[columns["inflorescencetype"]]);b=clean(df[columns["polmodbinary"]])
    sub=pd.DataFrame({"species":s,"inflorescence":a,"pollination":b})
    sub=sub[sub.species.notna()]
    by=sub.groupby("species",sort=False)
    out["duplicate_species_rows"]=int(len(sub)-sub.species.nunique())
    conflicts=(by[["inflorescence","pollination"]].nunique(dropna=True)>1).any(axis=1)
    out["conflicting_species"]=int(conflicts.sum())
    good=sub[~sub.species.isin(conflicts[conflicts].index)]
    # Take each nonmissing state independently, only for nonconflicting
    # species; no genus imputation or interpretation of numeric codes.
    q=good.groupby("species",sort=False)[["inflorescence","pollination"]].first()
    p=q.dropna()
    out["unique_species_with_inflorescence"]=int(q.inflorescence.notna().sum())
    out["unique_species_with_pollination"]=int(q.pollination.notna().sum())
    out["unique_species_with_both"]=int(len(p))
    out["joint_raw_state_counts"]=[{"inflorescence":str(i),"pollination":str(j),"species":int(n)} for (i,j),n in p.value_counts().items()]
    if "datasourceofpolmod" in columns:
        v=clean(df[columns["datasourceofpolmod"]]);out["pollination_source_nonmissing_among_coded_rows"]=int((v.notna()&b.notna()).sum())
    out["status"]="raw_schema_and_overlap_only"
    return out

def self_test():
    x=pd.DataFrame({"Genus_species":["A a","A a","B b","C c","C c"],"InflorescenceType":["a","a","b","a","b"],"PolModbinary":["0","0","1","1","1"]})
    a=audit_table(x)
    assert a["unique_species_with_both"]==2 and a["conflicting_species"]==1
    print("EXTERNAL_INFLORESCENCE_AUDIT_TEST_PASS")

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    a.out.parent.mkdir(parents=True,exist_ok=True)
    out={"doi":DOI,"publisher_file":EXPECTED_FILE,"publisher_file_url":URL,"reported_species":194476,"reported_pollination_mode_species":179263,
         "warning":"Source schema audit only. Wind/animal is NOT detailed pollinator guild. Inflorescence determinacy is NOT flower tube depth or size. Publication totals are not independently verified until file read. Raw trait/source codes not interpreted before dictionary validation. No ecological contrast or own-project data."}
    try:
        raw=get(META,maxbytes=2_000_000);m=json.loads(raw)
        out["metadata"]={"sha256":hashlib.sha256(raw).hexdigest(),"title":m.get("title"),"version":m.get("versionNumber"),"license":m.get("license")}
        assert "inflorescence" in str(m.get("title","")).lower()
        b=get(URL,maxbytes=10_000_000)
        out["download_bytes"]=len(b);out["file_sha256"]=hashlib.sha256(b).hexdigest()
        file=a.out.parent/"source.RData";file.write_bytes(b)
        try:
            tables=pyreadr.read_r(str(file))
            out["objects"]={str(k):audit_table(v) for k,v in tables.items() if isinstance(v,pd.DataFrame)}
            out["status"]="retrieved" if out["objects"] else "no_supported_data_frame"
        finally:file.unlink(missing_ok=True)
    except Exception as e:
        out["status"]="unavailable_or_invalid";out["error_type"]=type(e).__name__;out["error"]=str(e)[:400]
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False))
    print("EXTERNAL_INFLORESCENCE_AUDIT",json.dumps(out,ensure_ascii=False,allow_nan=False),flush=True)
if __name__=="__main__":main()
