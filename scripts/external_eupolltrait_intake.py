#!/usr/bin/env python3
"""Publisher-pinned EuPollTrait schema and exact-name EuPPollNet coverage only."""
import argparse,collections,hashlib,io,json,pathlib,re,tempfile,urllib.request,zipfile
import pandas as pd
import numpy as np
from external_euppollnet_bee_loss import BASE,SOURCES,VISITOR,PLANT,download
from external_euppollnet_channel_retention import prepare
from external_floral_replication_intake import name

URL="https://zenodo.org/records/18032357/files/pollinator_traits_bees_hoverflies_repository.zip?download=1"
MD5="7a50cc03e2c5a287e663dc5d029f8026"
BEE={"andrenidae","apidae","colletidae","halictidae","megachilidae","melittidae","stenotritidae"}
NAMES={"dwc_taxon_core_bees_hoverflies.csv","dwc_mof_bees_hoverflies.csv","trait_metadata.csv","metadata_traits_mof.xlsx","traits_raw_bees_hoverflies.xlsx","readme.md","readme.txt"}

def norm(s):return re.sub("[^a-z0-9]","",str(s).casefold())
def valid(v):return str(v).strip().casefold() not in {"","na","nan","none","null","unknown","[missing]"}
def col(d,k):
    c=[v for v in d if norm(v)==norm(k)]
    if len(c)>1:raise ValueError("ambiguous column "+k)
    return c[0] if c else None
def binomial(s):
    if not valid(s):return None
    x=" ".join(str(s).replace("_"," ").strip().split())
    return name(x[0].upper()+x[1:]) if x else None

def get():
    with urllib.request.urlopen(urllib.request.Request(URL,headers={"User-Agent":"FD publisher source audit"}),timeout=90) as r:
        raw=r.read(12_000_001)
    assert len(raw)<=12_000_000 and hashlib.md5(raw).hexdigest()==MD5
    return raw

def csv(z,f):
    b=z.read(f)
    for encoding in ("utf-8-sig","cp1252"):
        try:
            s=b.decode(encoding,errors="strict")
            d=pd.read_csv(io.StringIO(s),sep=None,engine="python",dtype=str,keep_default_na=False)
            return d,encoding
        except UnicodeDecodeError:continue
    raise ValueError("publisher encoding unresolved")

def taxon_table(t):
    idc=col(t,"taxonID");assert idc,"taxonID missing"
    g,e=col(t,"genus"),col(t,"specificEpithet")
    sc=col(t,"scientificName");rank=col(t,"taxonRank");fam=col(t,"family")
    assert fam and (sc or (g and e)),"source name/family schema unresolved"
    rows=[];bad=collections.Counter()
    for _,x in t.iterrows():
        if rank and x[rank].strip().casefold()!="species":
            bad["non_species_rank"]+=1;continue
        if g and e and valid(x[g]) and valid(x[e]):
            nm=binomial(x[g]+" "+x[e])
            if sc:
                sn=binomial(x[sc])
                if sn and nm!=sn:
                    bad["exact_name_conflict"]+=1;continue
        else:nm=binomial(x[sc]) if sc else None
        if not nm or not valid(x[idc]):
            bad["unresolved_name_or_id"]+=1;continue
        family=x[fam].strip().casefold()
        group="bee" if family in BEE else ("hoverfly" if family=="syrphidae" else None)
        if group is None:
            bad["outside_expected_families"]+=1;continue
        rows.append((x[idc],nm,group,family))
    q=pd.DataFrame(rows,columns=["id","name","group","family"]).drop_duplicates()
    if q.empty:return q,dict(bad)
    conflict=set(q.loc[q.id.duplicated(keep=False),"id"]);bad["conflicting_ids"]=len(conflict)
    q=q.loc[~q.id.isin(conflict)]
    conflict=set(q.loc[q.name.duplicated(keep=False),"name"]);bad["duplicate_name_ids_quarantined"]=len(conflict)
    return q.loc[~q.name.isin(conflict)].copy(),dict(bad)

def prepare_network(raw,meta):
    d=prepare(raw,meta)
    d["name"]=d[VISITOR].map(binomial)
    # Network and trait-source group/family must agree; no synonym guessing.
    d["family"]=d.Pollinator_family.str.strip().str.casefold()
    d["group"]=np.where(d.family.isin(BEE)&d.Pollinator_order.eq("hymenoptera"),"bee",
                       np.where(d.family.eq("syrphidae")&d.Pollinator_order.eq("diptera"),"hoverfly","other"))
    d=d.loc[d.name.notna()&d.group.ne("other")].copy()
    # Quarantine normalized name collisions across source visitor labels.
    bad=set(d.groupby("name")[VISITOR].nunique().loc[lambda x:x>1].index)
    conflicts=set(d.groupby("name").family.nunique().loc[lambda x:x>1].index)
    return d.loc[~d.name.isin(bad|conflicts)],{"name_collisions":len(bad),"family_conflicts":len(conflicts)}

def network_support(d,names):
    x=d.loc[d.name.isin(names)].copy()
    if x.empty:return {"visit_rows":0,"visitor_species":0,"studies":0,"eligible_sites":0,"candidate_same_plant_pairs":0,"candidate_two_sided_turnover_pairs":0}
    counts=collections.Counter();focal=set();pair_studies=set()
    for (study,year),g in x.groupby(["Study_id","year"],sort=True):
        sites=[]
        for _,h in g.groupby(["Latitude","Longitude"],sort=True):
            p=h.groupby(PLANT).name.agg(set).to_dict()
            if len(p)>=5 and h.name.nunique()>=5:sites.append(p)
        counts["eligible_sites"]+=len(sites)
        for i,A in enumerate(sites):
            for B in sites[i+1:]:
                for p in A.keys()&B.keys():
                    counts["candidate_same_plant_pairs"]+=1
                    if A[p]-B[p] and B[p]-A[p]:
                        counts["candidate_two_sided_turnover_pairs"]+=1
                        focal.add(p);pair_studies.add(study)
    return {"visit_rows":len(x),"visitor_species":int(x.name.nunique()),"studies":int(x.Study_id.nunique()),
            "eligible_sites":counts["eligible_sites"],"candidate_same_plant_pairs":counts["candidate_same_plant_pairs"],
            "candidate_two_sided_turnover_pairs":counts["candidate_two_sided_turnover_pairs"],
            "turnover_focal_plant_names":len(focal),"turnover_studies":len(pair_studies)}

def measurements(m,t,d):
    idc,tc,vc=map(lambda k:col(m,k),["taxonID","measurementType","measurementValue"])
    assert idc and tc and vc,"MeasurementOrFact schema unresolved"
    x=m.rename(columns={idc:"id",tc:"kind",vc:"value"}).merge(t,on="id",how="inner",validate="many_to_one")
    x=x.loc[x.value.map(valid)].copy()
    unit=col(m,"measurementUnit");method=col(m,"measurementMethod");remarks=col(m,"measurementRemarks")
    # Full measurement-type inventory; no values interpreted as actual
    # morphology without reading original metadata (including sex/caste).
    out={}
    for (group,k),q in x.groupby(["group","kind"],sort=True):
        keys=set(q.name)
        linked=d.loc[d.group.eq(group)]
        overlap=keys&set(linked.name)
        family_map=t.set_index("name").family.to_dict()
        good=linked.loc[linked.name.map(family_map).eq(linked.family)]
        matched=keys&set(good.name)
        v=pd.to_numeric(q.value,errors="coerce")
        ent={"records":len(q),"unique_species":len(keys),"exact_overlap_species":len(overlap),
             "family_agree_overlap_species":len(matched),
             "positive_numeric_records":int((v.gt(0)&np.isfinite(v)).sum()),
             "multiple_rows_species":int(q.groupby("name").size().gt(1).sum())}
        if unit:ent["raw_units"]={str(k):int(v) for k,v in q[unit].value_counts().head(30).items()}
        if method:ent["method_nonempty"]=int(q[method].map(valid).sum())
        if remarks:ent["remarks_nonempty"]=int(q[remarks].map(valid).sum())
        if any(s in norm(k) for s in ("bodylength","intertegular","itd","hairlength","hairdensity","tonguelength","hairiness")):
            ent["candidate_network_support"]=network_support(good,matched)
        out.setdefault(group,{})[k]=ent
    return out

def self_test():
    t=pd.DataFrame({"taxonID":["a","b","c"],"genus":["Apis","Eristalis","Apis"],"specificEpithet":["mellifera","tenax","mellifera"],"family":["Apidae","Syrphidae","Apidae"],"taxonRank":["species"]*3})
    q,a=taxon_table(t)
    assert list(q.name)==["eristalis tenax"] and a["duplicate_name_ids_quarantined"]==1
    assert binomial("Bombus terrestris L.") is None
    print("FD_EUPOLLTRAIT_AUDIT_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    raw=get();info={"doi":"10.5281/zenodo.18032357","md5":MD5,"sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw)}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        files=[f for f in z.namelist() if not f.endswith("/")];info["members"]=files
        tables={}
        for target in ("dwc_taxon_core_bees_hoverflies.csv","dwc_mof_bees_hoverflies.csv","trait_metadata.csv"):
            ff=[f for f in files if pathlib.PurePosixPath(f).name.casefold()==target]
            if len(ff)!=1:
                if target=="trait_metadata.csv":continue
                raise ValueError("publisher member unresolved: "+target)
            d,enc=csv(z,ff[0]);tables[target]=d
            info[target]={"member":ff[0],"rows":len(d),"columns":list(d.columns),"encoding":enc,
                          "sha256":hashlib.sha256(z.read(ff[0])).hexdigest()}
        print("FD_EUPOLLTRAIT_SCHEMA",json.dumps(info,allow_nan=False),flush=True)
        t,ta=taxon_table(tables["dwc_taxon_core_bees_hoverflies.csv"])
        assert len(t)>0,"no exact publisher taxon keys"
        with tempfile.TemporaryDirectory() as td:
            w=pathlib.Path(td);pins={}
            for f,(sha,blob) in SOURCES.items():pins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
            n=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False);assert len(n)==623476
            d,na=prepare_network(n,pd.read_csv(w/"Metadata.csv",dtype=str))
        out={"source":info,"network_source_pins":pins,"taxonomy_audit":ta,"network_name_audit":na,
             "group_counts":{str(k):int(v) for k,v in t.group.value_counts().items()},
             "coverage":measurements(tables["dwc_mof_bees_hoverflies.csv"],t,d),
             "scope":"Publisher source and coverage audit only; no functional effect estimated. Bee ITD and hoverfly body length are not the same measure. Tongue length in this source is ITD-derived, not independent measurement. Exact-name and family concordance are not full taxonomic harmonization. Same EuPPollNet observations; no independent interaction replication. No raw data redistribution."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
        for g,tt in out["coverage"].items():print("FD_EUPOLLTRAIT_COVERAGE",g,json.dumps(tt,allow_nan=False),flush=True)
        print("FD_EUPOLLTRAIT_COMPLETE",json.dumps({"group_counts":out["group_counts"],"taxonomy_audit":ta,"network_name_audit":na}),flush=True)
if __name__=="__main__":main()
