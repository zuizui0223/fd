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
CONTROVERSIAL={"asteraceae","poaceae","apiaceae","brassicaceae"}
WANG_COLUMNS=["Order","Family_E","Genus_E","Species_E1","Species_E2","Symmtry","Note1","Note2"]
WANG_LABELS={"RadSymmetry":"radial","BiSymmetry":"bilateral"}
def wang_candidate(raw,genus):
    # Accept an exact binomial consistent with the published genus.
    # Epithet-only fields may be combined with that genus, never with
    # a guessed genus or a synonym. Conflicting valid names are excluded.
    raw=str(raw).strip().replace("_"," ")
    gen=str(genus).strip()
    direct=name(raw)
    if direct is not None:
        return direct if direct.split()[0]==gen.casefold() else None
    if re.fullmatch(r"[a-z][a-z-]+",raw):
        return name(gen+" "+raw)
    return None
def wang_keys(d):
    e1=[wang_candidate(x,g) for x,g in zip(d.Species_E1,d.Genus_E)]
    e2=[wang_candidate(x,g) for x,g in zip(d.Species_E2,d.Genus_E)]
    merged=[];counts=Counter()
    for x,y in zip(e1,e2):
        if x is not None and y is not None and x!=y:
            counts["conflicting_candidate_names"]+=1;merged.append(None)
        elif x is not None or y is not None:
            counts["both_agree" if x==y else "e1_only" if x is not None else "e2_only"]+=1
            merged.append(x if x is not None else y)
        else:
            counts["no_genus_consistent_binomial"]+=1;merged.append(None)
    assert sum(counts.values())==len(d)
    return pd.Series(merged,index=d.index,dtype=object),dict(counts)
def wang_tables(d):
    assert list(d.columns)==WANG_COLUMNS and len(d)==279877
    assert set(d.Symmtry)==set(WANG_LABELS)
    keys,census=wang_keys(d)
    v=d.Symmtry.map(WANG_LABELS)
    notes=d.Note1.str.strip()
    note_clear=notes.eq("")
    family=d.Family_E.str.strip().str.casefold()
    known_family=family.ne("")&~family.isin(["na","nan","unknown","none"])
    # Exclude the whole species across duplicate records when there is
    # any annotation, label conflict or family disagreement.
    bad_notes=set(keys[~note_clear].dropna())
    fam_conf=Counter()
    families={}
    for k,f in zip(keys,family):
        if k is not None:families.setdefault(k,set()).add(f)
    fam_bad={k for k,z in families.items() if len(z)!=1}
    tab,conf=consensus(zip(keys,v))
    strict={k:x for k,x in tab.items() if k not in bad_notes and k not in fam_bad}
    nofam={k:x for k,x in strict.items() if len(families[k])==1 and next(iter(families[k])) not in CONTROVERSIAL|{"","na","nan","unknown","none"}}
    out={"key_census":census,"exact_binomial_rows":int(keys.notna().sum()),"conflicting_symmetry_names":conf,"family_conflicts":len(fam_bad),"note1_counts":dict(Counter(notes).most_common(30)),"note1_annotated_distinct_names":len(bad_notes),"explicit_reference_note_rows":int(d.Note2.map(valid).sum()),"family_known_rows":int(known_family.sum())}
    return {"all_publisher_binary":tab,"unannotated":strict,"unannotated_without_four_families":nofam},out
def audit_table(d,nhm):
    tabs,info=wang_tables(d)
    out={"status":"audited_publisher_binary","rows":len(d),"columns":list(d.columns),"raw_state_counts":dict(Counter(d.Symmtry)),**info,"scopes":{}}
    for scope,t in tabs.items():
        z={"distinct_exact_names":len(t),"states":dict(Counter(t.values())),"joint_coverage":{}}
        for a,b in nhm.items():
            ids=t.keys()&b.keys()
            z["joint_coverage"][a]={"species":len(ids),"radial":sum(t[x]=="radial" for x in ids),"bilateral":sum(t[x]=="bilateral" for x in ids)}
        for label,aa in [("flower_color_form",["flower colour","flower shape"]),("corolla_color_tube",["corolla colour","corolla_tube_range_cm"]),("flower_color_tube",["flower colour","corolla_tube_range_cm"])]:
            ids=set(t).intersection(*(nhm[a].keys() for a in aa))
            z[label+"_symmetry_common_species"]=len(ids)
        out["scopes"][scope]=z
    out["direct_evidence_classification"]="Unavailable for most rows. Binary publisher-coded states, not certified direct observations. Unannotated is a Note1 filter only; not a direct-evidence flag."
    return out
def self_test():
    assert wang_candidate("Abelia biflora","Abelia")=="abelia biflora"
    assert wang_candidate("Abelia_biflora","Abelia")=="abelia biflora"
    assert wang_candidate("biflora","Abelia")=="abelia biflora"
    assert wang_candidate("Other biflora","Abelia") is None
    d=pd.DataFrame({"Species_E1":["Abelia biflora","Abelia biflora","biflora"],"Species_E2":["Abelia_biflora","Abelia other",""],"Genus_E":["Abelia"]*3})
    k,c=wang_keys(d);assert list(k)==["abelia biflora",None,"abelia biflora"] and c["conflicting_candidate_names"]==1
    print("FD_WANG_SCHEMA_SELF_TEST_PASS",flush=True)
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=pathlib.Path);a=p.parse_args();self_test()
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
    out["warning"]="Publisher-only source coverage. Strict color lexical coding is intentionally conservative, no imputation. Publisher RadSymmetry and BiSymmetry labels mapped explicitly according to the study binary symmetry definitions; no direct evidence claim. Note1 exclusions are polymorphism/annotation sensitivity, not proof of species-level measurement. Exact names do not establish taxonomic identity or independent observations. Tube length is not tube depth. No island results, no raw taxon rows exported."
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    assert len(out["tables"])==1 and next(iter(out["tables"].values()))["status"]=="audited_publisher_binary","publisher schema not audited"
if __name__=="__main__":main()
