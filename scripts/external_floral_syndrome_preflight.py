#!/usr/bin/env python3
"""Independent publisher source intake only; no fitted biological contrasts.

Official GIFT 3.1 metadata and independently distributed NHM ADEPT
resource. Source independence from our assembled database does not imply
independent original botanical descriptions. Never use project flora/trait
rows or their derived aggregates. No authentication bypass or inferred
checklist completeness.
"""
import argparse,csv,hashlib,io,json,pathlib,re,urllib.parse
from collections import Counter
from tropical_pollinator_source_preflight import get

GIFT="https://gift.uni-goettingen.de/api/extended/index3.1.php"
NHM="https://data.nhm.ac.uk/api/3/action/resource_show?id=2552d157-9519-476c-995d-906e2433cc19"
DRYAD="https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.hqbzkh1nh"

def sha(b):return hashlib.sha256(b).hexdigest()
def key(x):return re.sub(r"[^a-z0-9]","",x.casefold())
def tab(raw):
    x=json.loads(raw)
    if isinstance(x,dict):
        if "data" in x:x=x["data"]
        elif "result" in x:x=x["result"]
    if not isinstance(x,list) or not all(isinstance(v,dict) for v in x):raise ValueError("unexpected table structure")
    return x

def describe(rows):
    cols=sorted({c for r in rows for c in r})
    return {"rows":len(rows),"columns":cols,"nonempty":{c:sum(r.get(c) not in (None,"") for r in rows) for c in cols}}

def gift():
    out={}
    for q in ("versions","regions","lists","traits_meta","reference_traits"):
        try:
            url=GIFT+"?"+urllib.parse.urlencode({"query":q})
            raw=get(url,maxbytes=80_000_000)
            rows=tab(raw);d=describe(rows)
            if q=="traits_meta":
                # Metadata only: published trait definitions and IDs.
                d["floral_trait_metadata"]=[r for r in rows if re.search(r"flower|floral|corolla|pollinat|selfing|self.compat", " ".join(map(str,r.values())),re.I)]
            if q in ("regions","lists"):
                # Report coded fields, not inferred native completeness.
                d["categorical_field_counts"]={}
                for c in d["columns"]:
                    if re.search(r"island|mainland|complete|native|endemic|floristic|taxon|list.type|entity.type",c,re.I):
                        counts=Counter(str(r.get(c)) for r in rows)
                        if len(counts)<=30:d["categorical_field_counts"][c]=dict(counts)
            out[q]={"status":"retrieved","url":url,"sha256":sha(raw),**d}
        except Exception as e:
            out[q]={"status":"unavailable_or_invalid","error_type":type(e).__name__,"detail":str(e)[:180]}
    return out

def nhm():
    raw=get(NHM,maxbytes=1_000_000);a=json.loads(raw)
    assert a["success"] is True
    r=a["result"];assert r["id"]=="2552d157-9519-476c-995d-906e2433cc19"
    assert "combined.csv" in str(r.get("name","")).casefold() or "combined.csv" in str(r.get("url","")).casefold()
    url=r["url"]
    parsed=urllib.parse.urlparse(url)
    if parsed.scheme!="https" or parsed.hostname!="data.nhm.ac.uk":raise ValueError("unreviewed source host")
    b=get(url,maxbytes=40_000_000)
    reader=csv.DictReader(io.StringIO(b.decode("utf-8-sig")))
    cols=reader.fieldnames
    assert cols and len(cols)==len(set(cols))
    relevant=[c for c in cols if re.search(r"flower|floral|corolla|petal|pollinat|taxon|species|name|wfo|genus|family",c,re.I)]
    # Explicit coverage census only, no mapping generic corolla/petal
    # dimensions to whole-flower size without the source dictionary.
    counts=Counter();overlap=Counter();n=0
    color=[c for c in relevant if re.search(r"flower|corolla|petal",c,re.I) and re.search(r"colou?r",c,re.I)]
    size=[c for c in relevant if re.search(r"flower|corolla|petal",c,re.I) and re.search(r"size|length|width|diameter|depth",c,re.I)]
    def valid(v):return v is not None and v.strip().casefold() not in ("","na","nan","null","none")
    for row in reader:
        if None in row:raise ValueError("malformed record")
        n+=1
        for c in relevant:counts[c]+=int(valid(row[c]))
        for c in color:
            if valid(row[c]):
                for s in size:
                    if valid(row[s]):overlap[(c,s)]+=1
    return {"status":"retrieved","dataset_doi":"10.5519/p3dm31kc","resource_id":r["id"],"resource_url":url,"metadata_sha256":sha(raw),"csv_sha256":sha(b),"bytes":len(b),"rows":n,"columns":cols,"floral_and_taxon_nonmissing":dict(counts),
            "color_size_column_record_overlap":[{"color_column":c,"size_column":s,"records":v} for (c,s),v in sorted(overlap.items())],
            "warning":"Record overlap is not verified same-species, same-structure measurements or comparable whole-flower size. Extraction and source-description overlap must be audited before any inference."}

def antilles():
    raw=get(DRYAD,maxbytes=2_000_000);d=json.loads(raw)
    return {"status":"retrieved","doi":"10.5061/dryad.hqbzkh1nh","metadata_sha256":sha(raw),
        "title":d.get("title"),"version":d.get("versionNumber"),"links":d.get("_links"),"license":d.get("license"),
        "warning":"Documented observed-versus-inferred pollination labels in publisher README; later measured morphology validation must use confirmed observations only. No raw images or inferred-label scores downloaded."}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args()
    r={"official_gift_3_1":gift()}
    for k,fn in (("nhm_adept",nhm),("antillean_gesneriaceae",antilles)):
        try:r[k]=fn()
        except Exception as e:r[k]={"status":"unavailable_or_invalid","error_type":type(e).__name__,"detail":str(e)[:180]}
    r["scope"]="Source-only readiness audit. No statistical estimates or new biological result."
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False,allow_nan=False))
    print("EXTERNAL_FLORAL_PREFLIGHT",json.dumps(r,ensure_ascii=False,allow_nan=False),flush=True)
if __name__=="__main__":main()
