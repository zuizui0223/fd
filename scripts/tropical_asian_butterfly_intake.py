#!/usr/bin/env python3
"""Audit pinned revised tropical Asian butterfly occurrence archive (v1.1, 2026)."""
import argparse,csv,hashlib,json,pathlib,urllib.request,io,tempfile,math,collections,re
from tropical_pollinator_source_preflight import get

FILE_ID=66997676
MD5="b99c4f94e3781da5cae25035703c8b45"
BYTES=347922325

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=pathlib.Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    meta=json.loads(get("https://api.figshare.com/v2/articles/25037645"))
    fs=meta.get("files",[])
    f=next(f for f in fs if f["id"]==FILE_ID)
    assert f["size"]==BYTES and (f.get("computed_md5") or f.get("supplied_md5"))==MD5
    print("BUTTERFLY_PIN",json.dumps(f),flush=True)
    md5=hashlib.md5();sha=hashlib.sha256();n=0
    url=f.get("downloadUrl") or f.get("download_url") or ("https://figshare.com/ndownloader/files/"+str(FILE_ID))
    with tempfile.TemporaryDirectory() as td:
        path=pathlib.Path(td)/"occurrences.csv"
        req=urllib.request.Request(url,headers={"User-Agent":"fd-research/1.0"})
        with urllib.request.urlopen(req,timeout=120) as r,path.open("wb") as w:
            while True:
                chunk=r.read(1024*1024)
                if not chunk: break
                n+=len(chunk)
                if n>BYTES: raise RuntimeError("file too large")
                md5.update(chunk);sha.update(chunk);w.write(chunk)
        assert n==BYTES and md5.hexdigest()==MD5
        # Inspect metadata and schema without fabricating an occurrence-to-island join.
        with path.open("r",encoding="utf-8-sig",errors="replace",newline="") as stream:
            reader=csv.reader(stream);header=next(reader)
            keys={re.sub(r"[^a-z0-9]","",x.casefold()):i for i,x in enumerate(header)}
            def find(options): return next((keys[x] for x in options if x in keys),None)
            si=find(["scientificname","species","acceptedscientificname","speciesname","speciesbinomial"])
            la=find(["decimallatitude","latitude","lat"])
            lo=find(["decimallongitude","longitude","lon","long"])
            src=find(["source","datasource","recordsource","sourcetype","dataorigin"])
            rows=bad=coord=0;sp=set();sources=collections.Counter()
            for row in reader:
                if not row or not any(x.strip() for x in row): continue
                rows+=1
                if len(row)!=len(header): bad+=1;continue
                if si is not None and row[si].strip(): sp.add(row[si].strip())
                if src is not None: sources[row[src].strip()]+=1
                if la is not None and lo is not None:
                    try:
                        x=float(row[la]);y=float(row[lo])
                        if math.isfinite(x) and math.isfinite(y) and -90<=x<=90 and -180<=y<=180:coord+=1
                    except ValueError: pass
        result={"file":f["name"],"file_id":FILE_ID,"bytes":n,"md5":md5.hexdigest(),"sha256":sha.hexdigest(),
                "columns":header,"data_rows":rows,"malformed_rows":bad,"species_column":header[si] if si is not None else None,
                "unique_species_strings":len(sp) if si is not None else None,
                "coordinate_columns":[header[la],header[lo]] if la is not None and lo is not None else None,
                "valid_coordinate_rows":coord if la is not None and lo is not None else None,
                "source_column":header[src] if src is not None else None,"source_counts_top30":sources.most_common(30),
                "warning":"revised source data, not 2025 published row count; no inferred island presence or absence"}
        (a.out/"occurrence_schema.json").write_text(json.dumps(result,indent=2,ensure_ascii=False))
        print("BUTTERFLY_OCCURRENCE_SCHEMA",json.dumps(result,ensure_ascii=False),flush=True)
        if bad or rows<500000: raise RuntimeError("source occurrence census requires review")
if __name__=="__main__":
    main()
