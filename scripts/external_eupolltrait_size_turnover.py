#!/usr/bin/env python3
"""External-only size similarity of exclusive flower visitors on matched carrier subwebs."""
import argparse, collections, hashlib, io, itertools, json, math, pathlib, re, tempfile, zipfile
import numpy as np
import pandas as pd
from external_eupolltrait_intake import get,csv,taxon_table,prepare_network,binomial,norm
from external_euppollnet_bee_loss import BASE,SOURCES,PLANT,download
from external_euppollnet_curveball import simulate
from external_euppollnet_spatial_scaling import infer,assign_genus,fixed_weights

SEED=20261128
DRAWS=399
MIN_ACCEPT=30
BEE_TRAIT="intertegular_distance"
FLY_TRAIT="Body_length_(mm,_raw_data)"
NUM=r"[0-9]+(?:[.,][0-9]+)?"

def parse_range(v):
    s=str(v).strip()
    # A decimal comma is accepted only within a single numeric token.
    m=re.fullmatch(r"("+NUM+r")(?:\s*(?:[-–—]|to)\s*("+NUM+r"))?\s*(?:mm)?",s,re.I)
    if not m:return None
    a=float(m[1].replace(",","."));b=float((m[2] or m[1]).replace(",","."))
    if not (math.isfinite(a) and math.isfinite(b) and 0<a<=b):return None
    return a,b

def caste(sex,column):
    s=str(sex).strip().casefold();c=str(column).strip().casefold()
    # Source sex/caste labels and original column, not an inference from species.
    # No unsexed values or queen values in primary.
    tok=set(re.findall(r"[a-z]+",s+" "+c))
    if "queen" in tok or "queens" in tok or s in {"male","m","males"} or "male" in tok or "males" in tok:return None
    if tok&{"worker","workers"} or s in {"w"}:return "worker"
    if s in {"female","females","f"} or tok&{"female","females"}:return "female"
    return None

def load_traits():
    raw=get()
    pins={"archive_sha256":hashlib.sha256(raw).hexdigest()}
    assert pins["archive_sha256"]=="b6e6b2c0b9f4c348e11185382d3f604610d1aad4181b0088b18a146f6180f782"
    tables=[]
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for target,expected in [
            ("dwc_taxon_core_bees_hoverflies.csv","17e487ecfe775858b786ca0e4e0e60529d5b7adfc670d18bde00d8f2ee193a2c"),
            ("dwc_mof_bees_hoverflies.csv","4c8f13913775483d458ed5e4f7946e6518873ed3b455a209b775bfa9c1cc310d"),
            ("trait_metadata.csv","fae52c364472e36958e49864d509fb3e8c7828aee8b795ae6a4c7ce470292378")]:
            p=[f for f in z.namelist() if pathlib.PurePosixPath(f).name.casefold()==target]
            assert len(p)==1
            digest=hashlib.sha256(z.read(p[0])).hexdigest();assert digest==expected
            t,enc=csv(z,p[0]);tables.append(t);pins[target]={"sha256":digest,"encoding":enc}
    tax,_=taxon_table(tables[0]);m=tables[1].rename(columns={"taxonID":"id"})
    m=m.merge(tax,on="id",validate="many_to_one")
    return tax,m,tables[2],pins

def trait_table(m,group):
    kind=BEE_TRAIT if group=="bee" else FLY_TRAIT
    x=m.loc[m.group.eq(group)&m.measurementType.eq(kind)].copy()
    assert len(x)>0
    census={"raw_rows":len(x),"sex":x.sex.value_counts(dropna=False).to_dict(),
       "original_column":x.originalColumn.value_counts(dropna=False).to_dict(),
       "life_stage":x.lifeStage.value_counts(dropna=False).to_dict(),
       "units":x.measurementUnit.value_counts(dropna=False).to_dict(),
       "measurement_source_nonempty":int(x.measurementSource.str.strip().ne("").sum())}
    count=collections.Counter();tab={};bounds={}
    social=m.loc[m.group.eq("bee")&m.measurementType.eq("Sociality")]
    eu=set(social.loc[social.measurementValue.str.strip().str.casefold().eq("eusocial"),"name"])
    for nm,q in x.groupby("name",sort=True):
        candidates=[]
        for _,v in q.iterrows():
            stage=str(v.lifeStage).strip().casefold()
            if stage not in {"adult","adults","imago","na","nan","none",""}:
                count["nonadult_or_unknown_explicit_stage"]+=1;continue
            unit=str(v.measurementUnit).strip().casefold()
            # ITD NA unit is documented mm in paper; hoverfly requires mm.
            if unit not in ({"mm","na","nan",""} if group=="bee" else {"mm"}):
                count["unsupported_unit"]+=1;continue
            y=parse_range(v.measurementValue)
            if y is None:
                count["unparsed_value"]+=1;continue
            if group=="bee":
                ca=caste(v.sex,v.originalColumn)
                if ca is None:
                    count["excluded_sex_caste"]+=1;continue
                if y[0]!=y[1]:
                    count["nonscalar_itd"]+=1;continue
            else:ca="literature"
            candidates.append((ca,y))
        if not candidates:
            count["no_candidate_species"]+=1;continue
        if group=="bee":
            workers=[y for c,y in candidates if c=="worker"]
            if workers:vals=workers;count["worker_species"]+=1
            elif nm in eu:
                count["eusocial_without_worker"]+=1;continue
            else:vals=[y for c,y in candidates if c=="female"];count["female_species"]+=1
        else:vals=[y for c,y in candidates]
        vals=set(vals)
        if len(vals)!=1:
            count["conflicting_or_missing_preferred_values"]+=1;continue
        lo,hi=next(iter(vals))
        bounds[nm]=(lo,hi)
        tab[nm]=(math.log2(lo)+math.log2(hi))/2
    census["selection"]=dict(count);census["retained_species"]=len(tab)
    census["value_range_mm"]=[min(v[0] for v in bounds.values()),max(v[1] for v in bounds.values())] if bounds else None
    return tab,bounds,census

def metrics(a,b,values,fams):
    u=a-b;v=b-a
    if not u or not v:return None
    i=np.array(sorted(u));j=np.array(sorted(v))
    return np.array([np.abs(values[i,None]-values[j][None,:]).mean(),(fams[i,None]==fams[j][None,:]).mean()])

def conditional(a,b,sa,sb,values,fams):
    obs=metrics(a,b,values,fams)
    if obs is None:return None
    h=len(a&b);accepted=[]
    for k,(x,y) in enumerate(zip(sa,sb)):
        assert len(x)==len(a) and len(y)==len(b)
        if len(x&y)==h:
            z=metrics(x,y,values,fams);assert z is not None
            accepted.append((k,z))
    if len(accepted)<MIN_ACCEPT or sum(k<len(sa)//2 for k,_ in accepted)<10 or sum(k>=len(sa)//2 for k,_ in accepted)<10:return {"accepted":len(accepted),"eligible":False}
    z=np.array([v for _,v in accepted]);null=z.mean(0)
    return {"accepted":len(z),"eligible":True,"size_obs":float(obs[0]),"size_null":float(null[0]),"size_gain":float(null[0]-obs[0]),
            "family_obs":float(obs[1]),"family_null":float(null[1]),"family_gain":float(obs[1]-null[1])}

def network_rows(d,tab,group,tax):
    fm=tax.set_index("name").family.to_dict()
    q=d.loc[d.group.eq(group)&d.name.isin(tab)&d.name.map(fm).eq(d.family)].copy()
    ids={k:i for i,k in enumerate(sorted(q.name.unique()))}
    values=np.array([tab[k] for k in ids]);fams=np.array([fm[k] for k in ids])
    rng=np.random.default_rng(SEED);records=[];diag=collections.Counter();accepted=[]
    for (study,year),g in q.groupby(["Study_id","year"],sort=True):
        sites={}
        for coord,s in g.groupby(["Latitude","Longitude"],sort=True):
            edges=s[[PLANT,"name"]].drop_duplicates()
            if edges[PLANT].nunique()<5 or edges.name.nunique()<5:continue
            p=edges.groupby(PLANT).name.agg(set).to_dict();names=sorted(p)
            rows=[{ids[x] for x in p[k]} for k in names]
            sims,info=simulate(rows,DRAWS,rng)
            sites[coord]=({k:i for i,k in enumerate(names)},rows,list(zip(*sims)))
            diag["eligible_sites"]+=1;diag["frozen_sites"]+=int(info["sample_changes"]==0)
        for (ca,(ia,ra,sa)),(cb,(ib,rb,sb)) in itertools.combinations(sites.items(),2):
            for p in sorted(ia.keys()&ib.keys()):
                a,b=ra[ia[p]],rb[ib[p]]
                if not a-b or not b-a:continue
                diag["candidate_pairs"]+=1
                z=conditional(a,b,sa[ia[p]],sb[ib[p]],values,fams)
                accepted.append(z["accepted"])
                if not z["eligible"]:
                    diag["rejected_conditional_support"]+=1;continue
                records.append(dict(z,study=study,year=int(year),plant=p,site_a=ca,site_b=cb,min_exclusive=min(len(a-b),len(b-a))))
    return pd.DataFrame(records),{"network_visit_rows":len(q),"visitor_species":len(ids),"studies_with_visits":int(q.Study_id.nunique()),**dict(diag),
        "accepted_draw_quantiles":np.quantile(accepted,[0,.25,.5,.75,1]).tolist() if accepted else []}

def report(rows):
    if rows.empty:return {"status":"no_eligible_pairs","supported":False}
    r,audit=assign_genus(rows)
    if r.empty:return {"status":"no_exact_plant_names","name_audit":audit,"supported":False}
    cols=["size_obs","size_null","size_gain","family_obs","family_null","family_gain"]
    q=r.groupby(["study","year","plant","genus"],as_index=False)[cols].mean()
    w=fixed_weights(q)
    size=infer(q,"size_gain");family=infer(q,"family_gain")
    pairs=int(r[["study","year","site_a","site_b"]].drop_duplicates().shape[0])
    # Original gate plus 30 study-year site pairs; direction predicted a priori.
    return {"name_audit":audit,"comparisons":len(r),"site_pairs":pairs,"study_equal":{k:float(np.average(q[k],weights=w)) for k in cols},
      "size":size,"family_descriptive":family,"supported":bool(size["supported"] and pairs>=30)}

def self_test():
    assert parse_range("3–5") == (3.,5.) and parse_range("2,5")== (2.5,2.5) and parse_range("5-3") is None
    assert caste("female","ITD_worker")=="worker" and caste("female","ITD_queen") is None
    assert caste("male","ITD_female") is None and caste("female","ITD")=="female"
    T=np.array([0.,1.,3.]);F=np.array(["a","a","b"])
    assert np.allclose(metrics({0},{1,2},T,F),[2.,.5])
    print("FD_EUPOLLTRAIT_SIZE_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    tax,m,metadata,pins=load_traits();print("FD_EUPOLLTRAIT_TRAIT_METADATA",json.dumps(metadata.loc[metadata.trait.astype(str).map(norm).str.contains("intertegular|bodylength"),:].to_dict("records"),allow_nan=False),flush=True)
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);npins={}
        for f,(sha,blob) in SOURCES.items():npins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
        raw=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False);assert len(raw)==623476
        d,na=prepare_network(raw,pd.read_csv(w/"Metadata.csv",dtype=str))
        results={}
        for g in ("bee","hoverfly"):
            tab,bounds,audit=trait_table(m,g)
            print("FD_EUPOLLTRAIT_SIZE_SELECTION",g,json.dumps(audit,allow_nan=False),flush=True)
            rows,diag=network_rows(d,tab,g,tax)
            result={"selection":audit,"network":diag,"primary":report(rows),
                    "both_exclusive_ge2":report(rows.loc[rows.min_exclusive.ge(2)] if not rows.empty else rows)}
            results[g]=result
            print("FD_EUPOLLTRAIT_SIZE_PANEL",g,json.dumps(result,allow_nan=False),flush=True)
        out={"source_pins":pins,"network_source_pins":npins,"network_name_audit":na,"panels":results,
             "scope":"Bee ITD primary, hoverfly literature size descriptive. Trait/group-specific complete-case local binary networks. Positive gain means more similar log2 sizes than fixed degree and shared-species conditional null; family concordance is descriptive on identical pairs (constant for Syrphidae). Crossed study/plant genus uncertainty, no independent interaction replication, effective pollination or island effect."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
        print("FD_EUPOLLTRAIT_SIZE_COMPLETE",json.dumps({"bee_supported":results["bee"]["primary"].get("supported",False),"hoverfly_descriptive":True}),flush=True)
if __name__=="__main__":main()
