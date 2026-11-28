#!/usr/bin/env python3
"""Cross-study transport of plant-associated visitor size relative to local availability."""
import argparse, collections, json, math, pathlib, tempfile
import numpy as np
import pandas as pd
from external_eupolltrait_size_turnover import load_traits,report
from external_eupolltrait_size_qc_family import corrected_traits,network
from external_eupolltrait_intake import prepare_network,binomial
from external_euppollnet_bee_loss import BASE,SOURCES,PLANT,download
from external_euppollnet_spatial_scaling import infer,fixed_weights

MIN_TRAIN_STUDIES=2
MIN_TRAIN_OTHER_SPECIES=2

def site_rows(d,tab,tax):
    fm=tax.set_index("name").family.to_dict()
    q=d.loc[d.group.eq("bee")&d.name.isin(tab)&d.name.map(fm).eq(d.family)].copy()
    assert q.name.nunique()==484 and len(q)==45164
    q["plant_key"]=q[PLANT].map(binomial)
    bad=set(q.dropna(subset=["plant_key"]).groupby("plant_key")[PLANT].nunique().loc[lambda x:x>1].index)
    # Keep ALL plants for site background, even unresolved focal names.
    rows=[];c=collections.Counter()
    for (study,year,lat,lon),g in q.groupby(["Study_id","year","Latitude","Longitude"],sort=True):
        e=g[[PLANT,"name","plant_key"]].drop_duplicates()
        if e[PLANT].nunique()<5 or e.name.nunique()<5:
            c["small_site_year"]+=1;continue
        c["eligible_site_year"]+=1
        pp=e.groupby(PLANT).name.agg(set).to_dict()
        degree=e.groupby("name")[PLANT].nunique().to_dict()
        keymap=e.drop_duplicates(PLANT).set_index(PLANT).plant_key.to_dict()
        for p,partners in pp.items():
            c["focal_site_records"]+=1
            k=keymap[p]
            if k is None or pd.isna(k) or k in bad:
                c["unresolved_focal"]+=1;continue
            # Exact leave-focal degree-weighted observed visitor pool.
            weights={v:n-int(v in partners) for v,n in degree.items() if n>int(v in partners)}
            if len(weights)<2:
                c["small_background"]+=1;continue
            baseline=sum(w*tab[v] for v,w in weights.items())/sum(weights.values())
            obs=float(np.mean([tab[v] for v in sorted(partners)]))
            rows.append({"study":study,"year":int(year),"plant":k,"genus":k.split()[0],"site":(float(lat),float(lon)),
                         "degree":len(partners),"observed":obs,"local":baseline,"relative":obs-baseline})
    return pd.DataFrame(rows),{**dict(c),"normalized_plant_collision_count":len(bad),"linked_species":int(q.name.nunique())}

def annual_cells(rows):
    if rows.empty:return rows
    keys=["study","year","plant","genus"]
    return rows.groupby(keys,as_index=False)[["observed","local","relative"]].mean()

def prediction(cells):
    """One held-out study at a time; never use a target study to fit predictions."""
    if cells.empty:return pd.DataFrame(),pd.DataFrame(),{}
    train_unit=cells.groupby(["study","plant","genus"],as_index=False).relative.mean()
    all_pred=[];matched=[];audit=collections.Counter()
    for held in sorted(cells.study.unique()):
        train=train_unit.loc[train_unit.study.ne(held)]
        target=cells.loc[cells.study.eq(held)]
        assert not train.study.eq(held).any()
        sums=train.groupby("plant").relative.agg(["mean","size"])
        genus_groups={g:x for g,x in train.groupby("genus",sort=True)}
        for row in target.itertuples(index=False):
            audit["target_plant_year_cells"]+=1
            if row.plant not in sums.index or sums.at[row.plant,"size"]<MIN_TRAIN_STUDIES:
                audit["insufficient_species_training_studies"]+=1;continue
            y=float(row.relative);pred=float(sums.at[row.plant,"mean"])
            a={"study":held,"year":row.year,"plant":row.plant,"genus":row.genus,
               "observed_relative":y,"species_prediction":pred,
               "local_only_squared_error":y*y,"species_squared_error":(y-pred)**2,
               "local_minus_species_mse":y*y-(y-pred)**2,
               "training_studies":int(sums.at[row.plant,"size"])}
            all_pred.append(a)
            # Genus comparator is fit from OTHER species only and OUTSIDE target study.
            gg=genus_groups[row.genus]
            gg=gg.loc[gg.plant.ne(row.plant)]
            if gg.plant.nunique()<MIN_TRAIN_OTHER_SPECIES or gg.study.nunique()<MIN_TRAIN_STUDIES:
                audit["insufficient_other_congeneric_training"]+=1;continue
            gp=float(gg.groupby("plant").relative.mean().mean())
            matched.append(dict(a,genus_prediction=gp,genus_squared_error=(y-gp)**2,
                   local_minus_genus_mse=y*y-(y-gp)**2,genus_minus_species_mse=(y-gp)**2-(y-pred)**2))
    return pd.DataFrame(all_pred),pd.DataFrame(matched),dict(audit)

def summarize(cells,value,extra):
    if cells.empty:return {"status":"no_eligible_cells","supported":False}
    w=fixed_weights(cells)
    metrics={k:float(np.average(cells[k],weights=w)) for k in extra}
    test=infer(cells,value)
    # This is a conditional bootstrap of frozen out-of-study predictions,
    # not an independent-study bootstrap with refitting inside each draw.
    return {"metrics":metrics,"inference":test,
            "training_study_count_quantiles":np.quantile(cells.training_studies,[0,.25,.5,.75,1]).tolist(),
            "supported":bool(test["supported"]),
            "conditional_prediction_bootstrap_only":True}

def analyze(rows):
    cells=annual_cells(rows)
    pred,matched,audit=prediction(cells)
    base=["observed_relative","species_prediction","local_only_squared_error","species_squared_error","local_minus_species_mse"]
    primary=summarize(pred,"local_minus_species_mse",base)
    secondary=summarize(matched,"genus_minus_species_mse",base+["genus_prediction","genus_squared_error","local_minus_genus_mse","genus_minus_species_mse"])
    # Check arithmetic and report coverage; secondary never rescues primary.
    if not pred.empty:
        assert np.allclose(pred.local_only_squared_error-pred.species_squared_error,pred.local_minus_species_mse)
    return {"source_plant_study_year_cells":len(cells),"source_focal_species":int(cells.plant.nunique()) if not cells.empty else 0,
            "source_studies":int(cells.study.nunique()) if not cells.empty else 0,
            "prediction_audit":audit,"primary":primary,"descriptive_other_species_genus_comparator":secondary}

def self_test():
    rows=[]
    for s,v in [("a",1.),("b",2.),("c",3.)]:
        for p,delta in [("g a",0.),("g b",1.),("g c",2.)]:
            rows.append({"study":s,"year":2020,"plant":p,"genus":"g","relative":v+delta})
    r=pd.DataFrame(rows);pred,matched,a=prediction(r)
    z=pred.loc[pred.study.eq("a")&pred.plant.eq("g a")].iloc[0]
    assert z.species_prediction==2.5 and z.training_studies==2
    assert len(pred)==len(matched)==9
    q=matched.loc[matched.study.eq("a")&matched.plant.eq("g a")].iloc[0]
    assert q.genus_prediction==4.
    print("FD_EUPOLLTRAIT_TRANSPORT_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    tax,m,meta,pins=load_traits();tab,selection=corrected_traits(m,"bee")
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);npins={}
        for f,(sha,blob) in SOURCES.items():npins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
        raw=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False);assert len(raw)==623476
        d,na=prepare_network(raw,pd.read_csv(w/"Metadata.csv",dtype=str))
        annual,diag=network(d,tab,"bee",tax)
        ar=report(annual)
        assert len(annual)==1261 and abs(ar["study_equal"]["size_gain"]-.11064886953661242)<1e-10
        rows,site_audit=site_rows(d,tab,tax)
        out={"source_pins":pins,"network_source_pins":npins,"selection":selection,"network_name_audit":na,
            "annual_replay":ar["study_equal"]["size_gain"],"site_audit":site_audit,"prediction":analyze(rows),
            "scope":"Post-discovery leave-one-study-out transport of site-relative mean log2 ITD. Same public sources, distinct estimand and broader same-study-year site cohort; target study excluded from fitted plant/genus means. Observed local background is leave-focal edge degree-weighted, not independently sampled availability. Frozen held-out predictions share training studies; crossed conditional uncertainty is not fully calibrated external validation. No independent dataset, causal floral filtering, actual pollination or island claim."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
        print("FD_EUPOLLTRAIT_TRANSPORT_RESULT",json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()
