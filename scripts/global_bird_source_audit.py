#!/usr/bin/env python3
"""Audit Baldwin/Myers global bird island source; no unverified use of eBird-derived data."""
import argparse,csv,hashlib,io,json,pathlib,re,zipfile
from tropical_pollinator_source_preflight import get

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    meta=json.loads(get("https://zenodo.org/api/records/12741516"))
    f=next(x for x in meta["files"] if x["key"]=="zenodo.zip")
    assert f["checksum"]=="md5:ce918e45ac29d7381bbc7f53b7f03923"
    raw=get(f["links"].get("content") or f["links"]["self"],maxbytes=180_000_000)
    assert len(raw)==f["size"] and hashlib.md5(raw).hexdigest()==f["checksum"].split(":")[1]
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        files=[{"path":m.filename,"bytes":m.file_size} for m in z.infolist() if not m.is_dir()]
        evidence=[];schemas=[]
        for m in z.infolist():
            if m.is_dir() or m.filename.startswith("__MACOSX/"):continue
            ext=pathlib.PurePosixPath(m.filename).suffix.lower()
            if ext in (".r",".rmd",".txt",".md") and m.file_size<5_000_000:
                s=z.read(m).decode("utf-8",errors="replace")
                hits=[{"line":i+1,"text":line[:400]} for i,line in enumerate(s.splitlines()) if re.search(r"ebird|birdlife|range.map|range map|avonet|island.*species|species.*island|read[._]csv|readRDS|load\(",line,re.I)]
                if hits:evidence.append({"path":m.filename,"matches":hits[:100]})
            if ext in (".csv",".tsv") and m.file_size<200_000_000:
                with z.open(m) as file:
                    stream=io.TextIOWrapper(file,encoding="utf-8-sig",errors="replace",newline="")
                    reader=csv.reader(stream,delimiter="\t" if ext==".tsv" else ",")
                    header=next(reader,[]);n=sum(1 for row in reader if row)
                schemas.append({"path":m.filename,"rows":n,"column_count":len(header),"columns":header[:80]})
        result={"zip_sha256":hashlib.sha256(raw).hexdigest(),"zip_bytes":len(raw),"inventory":files,"source_provenance_code":evidence,"table_schema":schemas,"warning":"Inventory only. Audit origin and native/resident filters before use; no eBird source assumed or used."}
        (a.out/"global_bird_archive_audit.json").write_text(json.dumps(result,indent=2))
        print("GLOBAL_BIRD_INVENTORY",json.dumps(files)[:24000],flush=True)
        print("GLOBAL_BIRD_SCHEMAS",json.dumps(schemas)[:24000],flush=True)
        print("GLOBAL_BIRD_PROVENANCE",json.dumps(evidence)[:30000],flush=True)
if __name__=="__main__":main()
