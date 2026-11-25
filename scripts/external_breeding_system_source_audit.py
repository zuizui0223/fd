#!/usr/bin/env python3
"""Public publisher supplementary-file intake for measured breeding systems.
No assembled project data. Do not execute supplied scripts or redistribute
source files. No interpretation of unknown labels before dictionary review.
"""
import argparse,collections,hashlib,io,json,pathlib,re,zipfile
import pandas as pd
from tropical_pollinator_source_preflight import get
SOURCES={
 "Razanajatovo2016":{"doi":"10.1038/ncomms13313","pmc":"PMC5095580","reported_taxa":1752,
 "scope":"Quantitative self-compatibility and autofertility indices from literature experiments; separate fruit and seed metrics."},
 "Ferrer2024":{"doi":"10.1093/aob/mcae056","pmc":"PMC11805948","reported_taxa":5685,
 "scope":"Species self-incompatibility compilation; partial/conflicting states and genus-inferred reaction mechanism not direct autonomous selfing."}
}
API="https://www.ebi.ac.uk/europepmc/webservices/rest/"

def sha(b):return hashlib.sha256(b).hexdigest()
def scalar(v):
    if pd.isna(v):return ""
    return str(v).strip()
def norm(v):return re.sub(r"[^a-z0-9]","",str(v).lower())
def names(v):
    # Exact literal binomial only, no fuzzy or synonym resolution.
    v=scalar(v)
    return v if re.fullmatch(r"[A-Z][a-zA-Z-]+ [a-z][a-z-]+",v) else None

def audit_sheet(df):
    columns=[str(c) for c in df.columns]
    report={"rows":len(df),"columns":columns}
    candidates={}
    for c in df.columns:
        s=df[c].map(scalar)
        k=norm(c)
        if re.search("species|taxon|binomial|scientificname|acceptedname",k):
            v={n for x in s if (n:=names(x)) is not None}
            candidates[str(c)]=len(v)
        if re.search("self|autofert|mating|compatib|isi|index|^af|^sc|reference|citation|^ref",k):
            vc=s[s!=""].value_counts()
            d={"nonmissing":int((s!="").sum()),"distinct":len(vc)}
            if len(vc)<=20:d["raw_state_counts"]={str(k):int(n) for k,n in vc.items()}
            report.setdefault("reproduction_or_reference_fields",{})[str(c)]=d
    report["literal_binomial_column_counts"]=candidates
    return report

def audit_workbook(raw):
    out={}
    xl=pd.ExcelFile(io.BytesIO(raw),engine="openpyxl")
    for sh in xl.sheet_names:
        df=xl.parse(sh)
        out[str(sh)]=audit_sheet(df)
    return out

def acquire(spec):
    url=API+spec["pmc"]+"/supplementaryFiles"
    out={**spec,"url":url}
    raw=get(url,maxbytes=40_000_000)
    out["archive_sha256"]=sha(raw);out["archive_bytes"]=len(raw)
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        infos=z.infolist()
        if sum(x.file_size for x in infos)>100_000_000:raise ValueError("expanded archive exceeds cap")
        out["archive_members"]=[{"name":x.filename,"bytes":x.file_size} for x in infos]
        files={}
        for x in infos:
            if not x.filename.lower().endswith(".xlsx"):continue
            if x.file_size>20_000_000:raise ValueError("workbook exceeds cap")
            b=z.read(x)
            files[x.filename]={"sha256":sha(b),"sheets":audit_workbook(b)}
        out["workbooks"]=files
    out["status"]="retrieved" if files else "no_xlsx_in_archive"
    return out

def self_test():
    d=pd.DataFrame({"Species":["Aaa aaa","Aaa aaa","Bbb bbb","Unknown"],"AF":[0,.5,1,None],"Self_compatibility":["SC","SC","SI",None]})
    a=audit_sheet(d)
    assert a["literal_binomial_column_counts"]["Species"]==2
    assert a["reproduction_or_reference_fields"]["AF"]["nonmissing"]==3
    print("EXTERNAL_BREEDING_SYSTEM_INTAKE_TEST_PASS")

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    out={"sources":{},"warning":"Published open supplementary materials via official Europe PMC archive service. Source-level schema audit only, not harmonized trait values, independent studies or biological effects. Reported taxa are publication targets. Do not map compatibility to autofertility, realized mating rate, or GIFT self-fertilization without reviewing measured definitions and primary references."}
    for name,s in SOURCES.items():
        try:v=acquire(s)
        except Exception as e:v={**s,"status":"unavailable_or_invalid","error_type":type(e).__name__,"detail":str(e)[:250]}
        out["sources"][name]=v
        print("EXTERNAL_BREEDING_SYSTEM_SOURCE",name,json.dumps(v,ensure_ascii=False,allow_nan=False),flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False))
if __name__=="__main__":main()
