#!/usr/bin/env python3
"""Tropical Asian butterfly original occurrence intake: inventory and checksum/schema audit only."""
import argparse,csv,hashlib,json,pathlib,urllib.request,io
from tropical_pollinator_source_preflight import get

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=pathlib.Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    meta=json.loads(get("https://api.figshare.com/v2/articles/25037645"))
    fs=meta.get("files",[])
    if not fs:
        raise RuntimeError("missing Figshare files")
    inv=[{k:f.get(k) for k in ("id","name","size","computed_md5","supplied_md5","downloadUrl")} for f in fs]
    (a.out/"figshare_inventory.json").write_text(json.dumps(inv,indent=2))
    print("BUTTERFLY_INVENTORY",json.dumps(inv)[:20000],flush=True)
    candidates=[f for f in fs if f["name"].lower().endswith(".csv") and "occurrence" in f["name"].lower() and "metadata" not in f["name"].lower()]
    if len(candidates)!=1:
        raise RuntimeError("unambiguous original occurrence CSV not found")
    f=candidates[0]
    if not 0<f["size"]<=250_000_000:
        raise RuntimeError("source occurrence file outside download cap")
    raw=get(f["downloadUrl"],maxbytes=250_000_000)
    assert len(raw)==f["size"],"size mismatch"
    md5=hashlib.md5(raw).hexdigest()
    expected=f.get("computed_md5") or f.get("supplied_md5")
    assert expected and md5==expected,"Figshare MD5 mismatch"
    for enc in ("utf-8-sig","latin-1"):
        try:
            txt=raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    reader=csv.reader(io.StringIO(txt))
    header=next(reader)
    if len(header)<4:
        raise RuntimeError("unexpected CSV columns")
    rows=0
    malformed=0
    for row in reader:
        if not row or not any(x.strip() for x in row):
            continue
        rows+=1
        malformed+=len(row)!=len(header)
    result={"file":f["name"],"bytes":len(raw),"md5":md5,"sha256":hashlib.sha256(raw).hexdigest(),
            "encoding":enc,"columns":header,"data_rows":rows,"malformed_rows":malformed,
            "published_rows":730190,"source":"Yau et al. 2025 Figshare 25037645",
            "warning":"source observation rows, not verified island presences or absences; no ecological analysis"}
    (a.out/"occurrence_schema.json").write_text(json.dumps(result,indent=2))
    print("BUTTERFLY_OCCURRENCE_SCHEMA",json.dumps(result),flush=True)
    if malformed or rows<500_000:
        raise RuntimeError("source occurrence census requires review")
if __name__=="__main__":
    main()
