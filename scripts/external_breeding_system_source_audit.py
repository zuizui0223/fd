#!/usr/bin/env python3
"""Public publisher supplementary-file intake for measured breeding systems.
No assembled project data. Do not execute supplied scripts or redistribute
source files. No interpretation of unknown labels before dictionary review.
"""
import argparse,collections,hashlib,io,json,pathlib,re,zipfile,urllib.parse,xml.etree.ElementTree as ET
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

CLOUD="https://pmc-oa-opendata.s3.amazonaws.com/"
S3NS={"s":"http://s3.amazonaws.com/doc/2006-03-01/"}

def publisher_cloud_versions(pmc):
    url=CLOUD+"?"+urllib.parse.urlencode({"list-type":"2","prefix":pmc+".","delimiter":"/","max-keys":"100"})
    b=get(url,maxbytes=500000)
    root=ET.fromstring(b)
    assert root.findtext("s:IsTruncated",namespaces=S3NS)=="false"
    v=[e.text for e in root.findall("s:CommonPrefixes/s:Prefix",S3NS)]
    assert all(re.fullmatch(re.escape(pmc)+r"\.\d+/",x or "") for x in v)
    return sorted(v),{"url":url,"sha256":sha(b),"versions":sorted(v)}

def cloud_media_url(u,prefix):
    p=urllib.parse.urlsplit(u)
    # Only objects explicitly exposed by the official article-version
    # metadata. No host substitutions except documented S3 -> HTTPS.
    if p.scheme=="s3" and p.netloc=="pmc-oa-opendata":path=p.path
    elif p.scheme=="https" and p.netloc=="pmc-oa-opendata.s3.amazonaws.com":path=p.path
    else:raise ValueError("unreviewed cloud media host")
    key=urllib.parse.unquote(path).lstrip("/")
    assert key.startswith(prefix) and ".." not in key.split("/")
    checksum=urllib.parse.parse_qs(p.query).get("md5",[])
    assert len(checksum)==1 and re.fullmatch(r"[a-fA-F0-9]{32}",checksum[0]),"official object MD5 missing"
    return CLOUD+urllib.parse.quote(key,safe="/"),key,checksum[0].lower()

def cloud_acquire(spec):
    out={**spec,"service":"NIH NLM NCBI PMC Article Datasets on AWS","service_documentation":"https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/","accessed_utc_date":"2026-11-25"}
    versions,listpin=publisher_cloud_versions(spec["pmc"]);out["version_listing"]=listpin
    records=[];files={}
    for prefix in versions:
        metaurl=CLOUD+prefix+prefix.rstrip("/")+".json"
        raw=get(metaurl,maxbytes=2000000);m=json.loads(raw)
        assert str(m.get("doi","")).lower()==spec["doi"].lower(),"article DOI mismatch"
        assert str(m.get("pmcid","")).upper()==spec["pmc"],"PMCID mismatch"
        license_code=str(m.get("license_code","")).lower()
        rec={"prefix":prefix,"metadata_url":metaurl,"metadata_sha256":sha(raw),
             "title":m.get("title"),"doi":m.get("doi"),"license_code":m.get("license_code"),
             "is_pmc_openaccess":m.get("is_pmc_openaccess"),"is_manuscript":m.get("is_manuscript"),
             "is_retracted":m.get("is_retracted"),"media_count":len(m.get("media_urls") or [])}
        records.append(rec)
        if str(m.get("is_pmc_openaccess")).lower() not in ("yes","true","1"):continue
        if str(m.get("is_retracted")).lower() in ("yes","true","1"):continue
        if not license_code.startswith(("cc","cc0")):continue
        for u in m.get("media_urls") or []:
            if not isinstance(u,str):raise ValueError("unexpected media URL schema")
            keypath=urllib.parse.urlsplit(u).path.lower()
            if not keypath.endswith((".xlsx",".xls",".csv",".tsv")):continue
            url,key,md5=cloud_media_url(u,prefix)
            b=get(url,maxbytes=20000000)
            assert hashlib.md5(b).hexdigest()==md5,"official media MD5 mismatch"
            ext=pathlib.Path(key).suffix.lower()
            if ext==".xlsx":
                sheets=audit_workbook(b)
            elif ext==".xls":
                xl=pd.ExcelFile(io.BytesIO(b),engine="xlrd")
                sheets={str(sh):audit_sheet(xl.parse(sh)) for sh in xl.sheet_names}
            else:
                sheets={"table":audit_sheet(pd.read_csv(io.BytesIO(b),sep="\t" if ext==".tsv" else ","))}
            files[key]={"url":url,"source_md5":md5,"sha256":sha(b),"bytes":len(b),"sheets":sheets,"license_code":license_code}
    out["article_versions"]=records
    out["workbooks"]=files
    out["status"]="retrieved" if files else "no_licensed_tabular_supplement_available"
    return out

def self_test():
    d=pd.DataFrame({"Species":["Aaa aaa","Aaa aaa","Bbb bbb","Unknown"],"AF":[0,.5,1,None],"Self_compatibility":["SC","SC","SI",None]})
    a=audit_sheet(d)
    assert a["literal_binomial_column_counts"]["Species"]==2
    assert a["reproduction_or_reference_fields"]["AF"]["nonmissing"]==3
    u,k,h=cloud_media_url("s3://pmc-oa-opendata/PMC1.1/a.xlsx?md5="+"a"*32,"PMC1.1/")
    assert k=="PMC1.1/a.xlsx" and h=="a"*32 and u==CLOUD+k
    try:cloud_media_url("https://example.com/PMC1.1/a.xlsx?md5="+"a"*32,"PMC1.1/")
    except ValueError:pass
    else:raise AssertionError("unreviewed host accepted")
    print("EXTERNAL_BREEDING_SYSTEM_INTAKE_TEST_PASS")

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");p.add_argument("--service",choices=("europe-pmc","official-pmc-cloud"),default="official-pmc-cloud");a=p.parse_args()
    if a.self_test:self_test();return
    out={"sources":{},"warning":"Published open supplementary materials via official PMC Cloud Service or Europe PMC archive service. Source-level schema audit only, not harmonized trait values, independent studies or biological effects. Reported taxa are publication targets. Do not map compatibility to autofertility, realized mating rate, or GIFT self-fertilization without reviewing measured definitions and primary references."}
    for name,s in SOURCES.items():
        try:v=cloud_acquire(s) if a.service=="official-pmc-cloud" else acquire(s)
        except Exception as e:v={**s,"status":"unavailable_or_invalid","error_type":type(e).__name__,"detail":str(e)[:250]}
        out["sources"][name]=v
        print("EXTERNAL_BREEDING_SYSTEM_SOURCE",name,json.dumps(v,ensure_ascii=False,allow_nan=False),flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False))
if __name__=="__main__":main()
