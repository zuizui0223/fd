#!/usr/bin/env python3
"""Source intake: complete candidate bird checklist x EltonTraits diet; bee archive schemas.
No geographic crosswalk, no pollinator effect. Missing taxonomy stays unknown.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import pathlib
import re
import unicodedata
import urllib.parse
import zipfile
from tropical_pollinator_source_preflight import BASE, PILOT, blobsha, get

EXTRA = {
    "Cape Verde": ("Baiser et al (2017) Cape Verde Birds_current_noAliens.csv", "098cd8b634df81e7a902f173e74e45fde319dbc1"),
    "Galapagos": ("Baiser et al (2017) Galapagos Birds_current.csv", "a3e6829301dda75d12637bec6ab0ed72dac68e09"),
    "Riau Lingga": ("Sin et al (2022) Riau Lingaa.csv", "b66b80974258c0c41273f504790e590b4508a12c"),
    "West Sumatra": ("Sin et al (2022) West Sumatra.csv", "86f6f1c5aa8e83cabb39c8ed379688ab914498eb"),
    "Ryukyu sensitivity": ("Kubota et al. Birds Ryukus_noAliens.csv", "e1083501ad6196d43c141a2d466ceff7fbeaa14d"),
}
MATRICES = {g: (x[0], x[4]) for g,x in PILOT.items()} | EXTRA
DIET_URLS = [
    "https://figshare.com/ndownloader/files/5631081",
    "https://figshare.com/ndownloader/files/5631084",
    "https://www.esapubs.org/archive/ecol/E095/178/BirdFuncDat.txt",
]
def norm_name(s):
    return " ".join(s.replace("_"," ").split()).casefold()

def namekey(s):
    s = unicodedata.normalize("NFKD", str(s)).casefold()
    return "".join(x for x in s if x.isalnum() and not unicodedata.combining(x))

def read_diet(raw):
    for encoding in ("utf-8-sig", "latin-1"):
        try:
            text=raw.decode(encoding)
        except UnicodeDecodeError:
            continue
        for delim in ("\t", ",", ";"):
            reader=csv.DictReader(io.StringIO(text), delimiter=delim)
            if not reader.fieldnames:
                continue
            cols={namekey(x):x for x in reader.fieldnames}
            if "scientific" not in cols or "dietnect" not in cols:
                continue
            data={}
            conflicting=set()
            for row in reader:
                name=norm_name(row.get(cols["scientific"]) or "")
                if not name:
                    continue
                try:
                    v=float(row.get(cols["dietnect"]) or "nan")
                except ValueError:
                    continue
                if not math.isfinite(v) or not 0<=v<=100:
                    continue
                if name in data and data[name]!=v:
                    conflicting.add(name)
                data[name]=v
            for name in conflicting:
                data.pop(name,None)
            if len(data)>=9000:
                return data, {"species_with_diet":len(data),"conflicts":len(conflicting),"encoding":encoding,"delimiter":repr(delim)}
    raise ValueError("not validated EltonTraits bird Scientific/Diet.Nect table")

def get_diet():
    errors=[]
    for url in DIET_URLS:
        try:
            raw=get(url)
            data, meta=read_diet(raw)
            return data, meta | {"url":url,"sha256":hashlib.sha256(raw).hexdigest()}
        except Exception as e:
            errors.append(str(e)[:180])
    raise RuntimeError("Cannot retrieve source diet: "+str(errors))

def matrix(raw):
    rows=list(csv.reader(io.StringIO(raw.decode("utf-8-sig"))))
    assert rows[0][0].strip().lower()=="species"
    islands=rows[0][1:]
    assert len(set(islands))==len(islands)
    species={}
    area=footer=None
    for row in rows[1:]:
        if not row or not row[0].strip():
            continue
        label=row[0].strip()
        assert len(row)==len(islands)+1
        if label=="Area (ha)":
            area=[float(x) for x in row[1:]]
        elif label=="sp.r":
            footer=[int(x) for x in row[1:]]
        else:
            assert label not in species and set(row[1:])<={"0","1"}
            species[label]=[int(x) for x in row[1:]]
    assert footer==[sum(x[j] for x in species.values()) for j in range(len(islands))]
    assert area is not None and all(x>0 for x in area)
    return islands,species,area,footer

def bird_intake(out):
    diet, meta=get_diet()
    rows=[]
    summaries={}
    unmatched=[]
    for group,(path,expected) in MATRICES.items():
        raw=get(BASE+urllib.parse.quote(path))
        if blobsha(raw)!=expected:
            raise ValueError("source blob mismatch "+group)
        islands,species,area,footer=matrix(raw)
        matches={s:diet[norm_name(s)] for s in species if norm_name(s) in diet}
        for s in species.keys()-matches.keys():
            unmatched.append({"archipelago":group,"species":s,"island_occurrences":sum(species[s])})
        for j,island in enumerate(islands):
            observed=[s for s,v in species.items() if v[j]]
            known=[matches[s] for s in observed if s in matches]
            unknown=len(observed)-len(known)
            rec={"archipelago":group,"source_island_name":island,"area_ha":area[j],
                 "all_birds":len(observed),"diet_known":len(known),"diet_unknown":unknown,
                 "diet_coverage":len(known)/len(observed) if observed else None,
                 "nectar_gt0_lb":sum(x>0 for x in known),
                 "nectar_ge20_lb":sum(x>=20 for x in known),
                 "nectar_ge50_lb":sum(x>=50 for x in known),
                 "nectar_ge50_ub":sum(x>=50 for x in known)+unknown,
                 "nectar_ge50_confirmed_zero":unknown==0 and not any(x>=50 for x in known)}
            rows.append(rec)
        sub=[x for x in rows if x["archipelago"]==group]
        summaries[group]={"islands":len(sub),"species":len(species),"matched_species":len(matches),
            "bird_occurrences":sum(x["all_birds"] for x in sub),"matched_occurrences":sum(x["diet_known"] for x in sub),
            "positive_ge50":sum(x["nectar_ge50_lb"]>0 for x in sub),
            "positive_ge20":sum(x["nectar_ge20_lb"]>0 for x in sub),
            "positive_gt0":sum(x["nectar_gt0_lb"]>0 for x in sub),
            "zero_ge50_complete":sum(x["nectar_ge50_confirmed_zero"] for x in sub),
            "unknown_ge50_presence":sum(x["nectar_ge50_lb"]==0 and x["diet_unknown"]>0 for x in sub),
            "sha256":hashlib.sha256(raw).hexdigest()}
    write_csv(out/"bird_island_diet_bounds.csv",rows)
    write_csv(out/"unmatched_bird_names.csv",unmatched)
    summary={"diet":meta,"groups":summaries,"island_columns":len(rows),
        "primary_candidate_island_columns":sum(x["archipelago"]!="Ryukyu sensitivity" for x in rows),
        "warning":"diet-based opportunity, not confirmed pollination; unresolved taxonomy bounded; island columns not necessarily unique geographic entities"}
    (out/"bird_diet_summary.json").write_text(json.dumps(summary,indent=2))
    print("BIRD_DIET_SUMMARY",json.dumps(summary),flush=True)

def write_csv(path,rows):
    if rows:
        with path.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

def table_schema(name,raw):
    result={"file":name,"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()}
    if name.lower().endswith(".xlsx"):
        from openpyxl import load_workbook
        book=load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
        result["sheets"]=[{"sheet":s.title,"rows":s.max_row,"columns":s.max_column,
            "first_rows":[[str(v)[:100] if v is not None else "" for v in row] for row in s.iter_rows(min_row=1,max_row=min(3,s.max_row),values_only=True)]} for s in book.worksheets]
    else:
        for enc in ("utf-8-sig","latin-1"):
            try:
                txt=raw.decode(enc)
                break
            except UnicodeDecodeError:
                continue
        try:
            delim=csv.Sniffer().sniff(txt[:20000],delimiters=",\t;").delimiter
        except csv.Error:
            delim="\t" if name.endswith(".tsv") else ","
        r=csv.reader(io.StringIO(txt),delimiter=delim)
        preview=[]
        count=0
        for row in r:
            count+=1
            if len(preview)<3:
                preview.append([v[:100] for v in row])
        result.update({"rows":count,"delimiter":repr(delim),"first_rows":preview})
    return result

def bee_intake(out):
    meta=json.loads(get("https://zenodo.org/api/records/23222636"))
    files=meta.get("files",[])
    if not files:
        raise ValueError("bee record has no files")
    inventory=[{k:x.get(k) for k in ("key","size","checksum","links")} for x in files]
    (out/"bee_file_inventory.json").write_text(json.dumps(inventory,indent=2))
    print("BEE_INVENTORY",json.dumps([{k:x.get(k) for k in ("key","size","checksum")} for x in files]),flush=True)
    schemas=[]
    cap=100_000_000
    for f in files:
        name=f["key"]
        size=f.get("size") or 0
        if not name.lower().endswith((".csv",".tsv",".xlsx",".zip",".txt")) or not 0<size<=cap:
            continue
        try:
            url=f.get("links",{}).get("content") or f.get("links",{}).get("self")
            if not url:
                raise ValueError("missing download link")
            raw=get(url,maxbytes=cap)
            if len(raw)!=size:
                raise ValueError("size mismatch")
            check=f.get("checksum","")
            if ":" in check:
                alg,expected=check.split(":",1)
                if hashlib.new(alg,raw).hexdigest()!=expected:
                    raise ValueError("archive checksum mismatch")
            if name.lower().endswith(".zip"):
                with zipfile.ZipFile(io.BytesIO(raw)) as z:
                    schemas.append({"file":name,"members":[{"name":x.filename,"bytes":x.file_size} for x in z.infolist()][:150]})
                    for m in z.infolist():
                        if m.file_size<=20_000_000 and m.filename.lower().endswith((".csv",".tsv",".xlsx")):
                            schemas.append(table_schema(name+"::"+m.filename,z.read(m)))
            elif name.lower().endswith((".csv",".tsv",".xlsx")):
                schemas.append(table_schema(name,raw))
            else:
                schemas.append({"file":name,"preview":raw[:1500].decode("utf-8",errors="replace")})
        except Exception as e:
            schemas.append({"file":name,"error":str(e)[:400]})
    (out/"bee_table_schemas.json").write_text(json.dumps(schemas,indent=2))
    print("BEE_SCHEMAS",json.dumps(schemas)[:24000],flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=pathlib.Path)
    p.add_argument("--stage",required=True,choices=["bird","bee","test"])
    a=p.parse_args()
    a.out.mkdir(exist_ok=True,parents=True)
    if a.stage=="bird":
        bird_intake(a.out)
    elif a.stage=="bee":
        bee_intake(a.out)
    else:
        assert norm_name(" A_b ")=="a b"
        assert matrix(b"species,X,Y\nA_b,1,0\nC_d,1,1\nArea (ha),5,10\nsp.r,2,1\n")[3]==[2,1]
        assert namekey("Diet.Nect")=="dietnect"

if __name__=="__main__":
    main()
