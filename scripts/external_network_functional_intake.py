#!/usr/bin/env python3
"""Audit independent published interaction sources before any ecological fit.

No own plant data, no source code execution, no guessed island labels.
The fixed publisher records are distinct candidates, not interchangeable tests.
"""
import argparse
import csv
import hashlib
import io
import json
import pathlib
import subprocess
import tempfile
import urllib.parse
import zipfile

import openpyxl

from tropical_pollinator_source_preflight import get

RECORDS = {
    "vollstadt_island_mainland": (15276347, "10.5281/zenodo.15276347"),
    "dore_global_networks": (4300427, "10.5281/zenodo.4300427"),
}
DATA_SUFFIX = {".csv", ".tsv", ".txt", ".xlsx", ".xls", ".rdata", ".rda", ".rds"}
MAX_FILE = 80_000_000
MAX_MEMBER = 50_000_000
MAX_TOTAL = 150_000_000

def sha(b):
    return hashlib.sha256(b).hexdigest()

def metadata_url(id):
    return "https://zenodo.org/api/records/" + str(id)

def public_url(id, key):
    # Documented publisher file route; never try private endpoints or tokens.
    return "https://zenodo.org/records/" + str(id) + "/files/" + urllib.parse.quote(key, safe="") + "?download=1"

def decode(raw):
    for enc in ("utf-8-sig", "cp1252"):
        try:
            return raw.decode(enc), enc
        except UnicodeDecodeError:
            continue
    raise ValueError("unsupported text encoding")

def table_info(path, raw):
    suffix = pathlib.PurePosixPath(path).suffix.lower()
    if suffix in {".csv", ".tsv", ".txt"}:
        s, enc = decode(raw)
        if pathlib.PurePosixPath(path).name.lower().startswith(("readme", "license")):
            return {"status":"documentation", "encoding":enc, "lines":len(s.splitlines())}
        sample = s[:8192]
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            return {"status":"unrecognized_delimiter","encoding":enc,"lines":len(s.splitlines())}
        reader = csv.reader(io.StringIO(s), dialect)
        header = next(reader, [])
        n = 0
        widths = {}
        for row in reader:
            if not row:
                continue
            n += 1
            widths[str(len(row))] = widths.get(str(len(row)), 0) + 1
        return {"status":"schema","encoding":enc,"columns":header[:250],"column_count":len(header),"rows":n,"row_widths":widths,"delimiter":dialect.delimiter}
    if suffix == ".xlsx":
        wb = openpyxl.load_workbook(io.BytesIO(raw), read_only=True, data_only=True)
        out = {}
        for ws in wb.worksheets:
            first = next(ws.iter_rows(min_row=1,max_row=1,values_only=True),())
            out[ws.title] = {"rows_including_header":ws.max_row,"columns":ws.max_column,"header":[str(x) if x is not None else "" for x in first[:250]]}
        wb.close()
        return {"status":"schema","sheets":out}
    if suffix in {".rds",".rdata",".rda"}:
        # Load data into a new environment; never source publisher R scripts.
        # Only structural metadata is emitted; nested lists receive a bounded census.
        expr = r'''
a <- commandArgs(TRUE)
inspect <- function(x, name, depth=0L) {
  if (is.data.frame(x) || is.matrix(x)) {
    cat("TABLE\t",name,"\t",nrow(x),"\t",ncol(x),"\t",paste(head(colnames(x),250),collapse="|"),"\n",sep="")
  } else if (is.list(x) && depth < 2L) {
    cat("LIST\t",name,"\t",length(x),"\t",paste(head(names(x),100),collapse="|"),"\n",sep="")
    nm <- names(x)
    for (i in seq_len(min(length(x),8L))) inspect(x[[i]],paste0(name,"/",if(is.null(nm)) i else nm[i]),depth+1L)
  } else cat("OBJECT\t",name,"\t",paste(class(x),collapse="|"),"\t",length(x),"\n",sep="")
}
if (a[2] == ".rds") inspect(readRDS(a[1]),"root") else {
  e <- new.env(parent=emptyenv())
  n <- load(a[1],envir=e)
  for (k in n) inspect(e[[k]],k)
}
'''
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d)/("data"+suffix)
            p.write_bytes(raw)
            r = subprocess.run(["Rscript","--vanilla","-e",expr,str(p),suffix],capture_output=True,text=True,timeout=90)
        if r.returncode:
            return {"status":"r_schema_failed","stderr":r.stderr[-1200:]}
        return {"status":"schema","objects":r.stdout.splitlines()[:250]}
    return {"status":"unsupported_format"}

def inspect_bytes(name, raw):
    if not name.lower().endswith(".zip"):
        return {name:table_info(name,raw)}
    out = {}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        infos=z.infolist()
        if sum(x.file_size for x in infos)>MAX_TOTAL:
            return {"archive":{"status":"uncompressed_cap","members":len(infos),"total_bytes":sum(x.file_size for x in infos)}}
        for item in infos:
            if item.is_dir():
                continue
            p=pathlib.PurePosixPath(item.filename)
            if p.is_absolute() or ".." in p.parts:
                raise ValueError("unsafe archive member name")
            # Do not extract arbitrary paths or run source scripts.
            if p.suffix.lower() not in DATA_SUFFIX or item.file_size>MAX_MEMBER:
                out[item.filename]={"status":"not_parsed","bytes":item.file_size}
                continue
            try:
                out[item.filename]=table_info(item.filename,z.read(item))
            except Exception as exc:
                out[item.filename]={"status":"schema_failed","error":str(exc)[:500]}
    return out

def verify_bytes(raw, info):
    alg, expected=info["checksum"].lower().split(":",1)
    if alg not in ("md5","sha256","sha1"):
        raise ValueError("unrecognized publisher checksum")
    if len(raw)!=int(info["size"]) or hashlib.new(alg,raw).hexdigest()!=expected:
        raise ValueError("publisher size/checksum mismatch")
    return {"bytes":len(raw),"publisher_checksum":info["checksum"],"sha256":sha(raw)}

def audit(id,doi):
    md=get(metadata_url(id),maxbytes=5_000_000)
    meta=json.loads(md)
    if int(meta["id"])!=id or meta.get("doi",meta.get("metadata",{}).get("doi"))!=doi:
        raise ValueError("publisher record identity mismatch")
    res={"doi":doi,"metadata_sha256":sha(md),"title":meta["metadata"]["title"],"license":meta["metadata"].get("license"),"files":[]}
    total=0
    for f in meta.get("files",[]):
        key=f["key"]
        row={"key":key,"declared_size":f.get("size"),"publisher_checksum":f.get("checksum")}
        res["files"].append(row)
        suffix=pathlib.PurePosixPath(key).suffix.lower()
        if suffix not in DATA_SUFFIX|{".zip"}:
            row["status"]="not_data_format";continue
        size=int(f["size"])
        if size>MAX_FILE or total+size>MAX_TOTAL:
            row["status"]="size_cap";continue
        try:
            raw=get(public_url(id,key),maxbytes=MAX_FILE)
            row.update(verify_bytes(raw,f))
            total+=size
            row["schema"]=inspect_bytes(key,raw)
            row["status"]="verified_bytes"
        except Exception as exc:
            row["status"]="unavailable_or_invalid"
            row["error"]=str(exc)[:500]
        print("FD_EXTERNAL_NETWORK_FILE",json.dumps(row,default=str),flush=True)
    res["verified_file_count"]=sum(x["status"]=="verified_bytes" for x in res["files"])
    return res

def self_test():
    b=b"plant,bird,visits\nA a,B b,2\nC c,B b,1\n"
    assert table_info("x.csv",b)["rows"]==2
    assert verify_bytes(b,{"checksum":"md5:"+hashlib.md5(b).hexdigest(),"size":len(b)})["bytes"]==len(b)
    try:
        verify_bytes(b,{"checksum":"md5:"+"0"*32,"size":len(b)})
    except ValueError:
        pass
    else:
        raise AssertionError("invalid digest accepted")
    w=io.BytesIO()
    with zipfile.ZipFile(w,"w") as z:
        z.writestr("network.csv",b)
        z.writestr("analysis.R","stop('must never execute')")
    out=inspect_bytes("source.zip",w.getvalue())
    assert out["network.csv"]["rows"]==2 and out["analysis.R"]["status"]=="not_parsed"
    print("FD_NETWORK_INTAKE_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=pathlib.Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    self_test()
    if a.self_test:
        return
    if a.out is None:
        p.error("--out is required")
    results={}
    for key,(id,doi) in RECORDS.items():
        try:
            results[key]=audit(id,doi)
        except Exception as exc:
            results[key]={"doi":doi,"status":"source_unavailable","error":str(exc)[:500]}
        print("FD_EXTERNAL_NETWORK_SOURCE",key,json.dumps(results[key],default=str),flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps({"sources":results,"ecological_effect_fitted":False,"warning":"Source identity and schema only. No own-data inputs, no inferred island status, no pollinator effectiveness claims."},indent=2,default=str))
    if not any(x.get("verified_file_count",0)>0 for x in results.values()):
        raise SystemExit("No publisher data file verified; no biological result.")
if __name__=="__main__":
    main()
