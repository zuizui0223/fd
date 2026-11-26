#!/usr/bin/env python3
"""Pinned Wang symmetry source + NHM independent publisher trait feasibility audit."""
import argparse,hashlib,io,json,pathlib,re,subprocess,tempfile
from collections import Counter
import pandas as pd
from tropical_pollinator_source_preflight import get,blobsha
from external_floral_replication_intake import NHM_ROOT,NHM,name,valid,column,consensus

ROOT="https://zenodo.org/records/6547106/files/"
FILES={"AngSplist_Symmetry20230625_Updated.rdata":"548c12e053d348696f4481d7d2c7104f","README.txt":"e3cea11344201ef4d385be0c6ecb8b75"}
COLORS={"red","pink","blue","yellow","orange","green","white","violet","purple","black","grey","gray","brown"}
FORMS={"campanulate","tubular","urceolate","ligulate","rotate","calyciform","infundibuliform"}
def color(x):
    if not valid(x):return None
    x=x.strip().casefold()
    return {"gray":"grey","purple":"violet"}.get(x,x) if x in COLORS else None
def length(x):
    if not valid(x):return None
    m=re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)\s*cm",x.strip().casefold())
    if not m:return None
    v=float(m[1]);return v if 0<v<=100 else None
def strict_nhm(raw):
    d=pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False);keys=d["taxon"].map(name)
    tab={};audit={}
    for c,fn in [("flower colour",color),("corolla colour",color),("flower shape",lambda x:x.strip().casefold() if x.strip().casefold() in FORMS else None)]:
        v=d[c].map(fn)
        tab[c],conf=consensus((k,x) for k,x in zip(keys,v) if x is not None)
        audit[c]={"accepted_species":len(tab[c]),"raw_nonmissing":int(d[c].map(valid).sum()),"conflict":conf,"states":dict(Counter(tab[c].values()))}
    a=d["corolla tube min. length [cm]"].map(length);b=d["corolla tube max. length [cm]"].map(length)
    keep=keys.notna()&a.notna()&b.notna()&a.le(b)
    length_table={}
    for i,x,y in zip(keys[keep],a[keep],b[keep]):
        if i in length_table and length_table[i]!=(x,y):raise ValueError("nonunique NHM length species")
        length_table[i]=(float(x),float(y))
    tab["corolla_tube_range_cm"]=length_table
    audit["corolla_tube_range_cm"]={"accepted_species":len(length_table),"raw_min":int(d["corolla tube min. length [cm]"].map(valid).sum()),"raw_max":int(d["corolla tube max. length [cm]"].map(valid).sum()),"parsed_both":int((a.notna()&b.notna()).sum()),"reversed_bounds":int((a.notna()&b.notna()&a.gt(b)).sum()),"max_length_quantiles":[float(z) for z in b[keep].quantile([0,.25,.5,.75,1])]}
    return tab,audit
def parse_symmetry(x):
    if not valid(x):return None
    z=x.strip().casefold()
    if z in {"actinomorphic","actinomorphy","actinomorph","radial","polysymmetric"}:return "radial"
    if z in {"zygomorphic","zygomorphy","zygomorph","bilateral","monosymmetric"}:return "bilateral"
    # Never interpret numeric codes or abbreviations without the actual publisher dictionary.
    return None
def audit_table(d,nhm):
    cols=list(d.columns)
    c=column(cols,["species","species_name","scientific_name","scientificName","accepted_name","taxon","taxon_name","name","binomial","sp","species_full","speciesname","species_name_accepted"])
    s=column(cols,["floral_symmetry","flower_symmetry","symmetry","floral symmetry","symmetry_type","floral_symmetry_type"])
    if c is None or s is None:
        return {"rows":len(d),"columns":cols,"status":"schema_unresolved","candidate_low_cardinality":{k:dict(Counter(d[k]).most_common(20)) for k in cols if d[k].nunique()<=20}}
    k=d[c].str.replace("_"," ",regex=False).map(name)
    raw=d[s]
    v=raw.map(parse_symmetry)
    t,conf=consensus((x,y) for x,y in zip(k,v) if y is not None)
    out={"status":"audited","rows":len(d),"columns":cols,"name_column":c,"symmetry_column":s,"raw_state_counts":dict(Counter(raw).most_common(30)),"exact_name_rows":int(k.notna().sum()),"accepted_word_label_species":len(t),"conflicting_species":conf,"source_fields":{},"joint_coverage":{}}
    for col in cols:
        if re.search(r"source|refer|origin|assign|infer|level|imput|method|evidence|data.?type|record",col,re.I):
            out["source_fields"][col]={"distinct_values":int(d[col].nunique()),"nonempty":int(d[col].map(valid).sum()),"top_values":dict(Counter(d[col]).most_common(15)) if d[col].nunique()<=30 else None}
    out["direct_evidence_classification"]="not established; inspect actual publisher README and columns"
    for a,b in nhm.items():
        ids=t.keys()&b.keys()
        out["joint_coverage"][a]={"species":len(ids),"radial":sum(t[x]=="radial" for x in ids),"bilateral":sum(t[x]=="bilateral" for x in ids)}
    ids=t.keys()&nhm["flower colour"].keys()&nhm["flower shape"].keys()
    out["color_shape_symmetry_common_species"]=len(ids)
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=pathlib.Path);a=p.parse_args()
    out={"publisher":"Wang et al., Zenodo 6547106, Dryad DOI 10.5061/dryad.ghx3ffbrh","source_files":{},"tables":{}}
    raw=get(NHM_ROOT+"combined.csv",maxbytes=25_000_000);assert blobsha(raw)==NHM["combined.csv"]
    nhm,audit=strict_nhm(raw);out["nhm_strict"]=audit
    out["nhm_joint_coverage"]={x+"__"+y:len(nhm[x].keys()&nhm[y].keys()) for i,x in enumerate(nhm) for y in list(nhm)[i+1:]}
    print("FD_STRICT_NHM",json.dumps({"traits":audit,"joint":out["nhm_joint_coverage"]},allow_nan=False),flush=True)
    with tempfile.TemporaryDirectory() as td:
        path=pathlib.Path(td)
        for n,h in FILES.items():
            b=get(ROOT+n+"?download=1",maxbytes=15_000_000);assert hashlib.md5(b).hexdigest()==h
            (path/n).write_bytes(b);out["source_files"][n]={"md5":h,"sha256":hashlib.sha256(b).hexdigest(),"bytes":len(b)}
        out["publisher_readme"]=(path/"README.txt").read_text(errors="replace")
        print("FD_WANG_README",json.dumps(out["publisher_readme"]),flush=True)
        # Base R load into an isolated environment; no execution of publisher scripts.
        expr='args<-commandArgs(TRUE); e<-new.env(parent=emptyenv()); nm<-load(args[1],envir=e); for(i in seq_along(nm)){x<-get(nm[i],envir=e,inherits=FALSE);cat(nm[i],class(x),dim(x),"\\n");if(is.data.frame(x)||is.matrix(x)){write.csv(as.data.frame(x),file=file.path(args[2],paste0("table_",i,".csv")),row.names=FALSE,na="")}}'
        proc=subprocess.run(["Rscript","--vanilla","-e",expr,str(path/next(iter(FILES))),str(path)],check=True,capture_output=True,text=True)
        out["r_objects"]=proc.stdout;print("FD_WANG_R_OBJECTS",json.dumps(proc.stdout),flush=True)
        for f in sorted(path.glob("table_*.csv")):
            d=pd.read_csv(f,dtype=str,keep_default_na=False)
            out["tables"][f.name]=audit_table(d,nhm)
            print("FD_WANG_TABLE",f.name,json.dumps(out["tables"][f.name],allow_nan=False),flush=True)
    out["warning"]="Publisher-only source coverage. Strict color lexical coding is intentionally conservative, no imputation. Symmetry numeric/abbreviated codes not interpreted without documentation; no direct evidence claim. Exact names do not establish taxonomic identity or independent observations. Tube length is not tube depth. No island results, no raw taxon rows exported."
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    assert out["tables"],"no tabular publisher objects"
if __name__=="__main__":main()
