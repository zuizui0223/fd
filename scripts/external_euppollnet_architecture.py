#!/usr/bin/env python3
"""External-only cross-site visitor-channel concordance and floral architecture."""
import argparse, collections, hashlib, json, math, pathlib, tempfile
import numpy as np
import pandas as pd
from external_euppollnet_bee_loss import BASE,SOURCES,download
from external_euppollnet_channel_retention import prepare,collect,summarize
from external_floral_replication_intake import name
from external_gift_floral_replication_coverage import source_tables

SEED=20261128
PRIMARY=("symmetry","tube_midpoint")
N_PERM=999
MIN_STUDIES=10
MIN_SPECIES=50
MIN_MOVABLE=30
MIN_GENERA=10

def stratum_table(records):
    r=records.copy()
    # No author truncation, synonym guessing or trait interpolation.
    r["name"]=r.plant.map(name)
    # Quarantine normalized-name collisions rather than treating aliases
    # as independent plants. This also keeps the original background intact.
    counts=r.dropna(subset=["name"]).groupby("name").plant.nunique()
    bad=set(counts[counts>1].index)
    r=r.loc[r.name.notna()&~r.name.isin(bad)].copy()
    r["turnover_class"]=np.where(r.min_exclusive.eq(1),"one_on_either_side","both_at_least_two")
    q=r.groupby(["study","year","turnover_class","name"],as_index=False).channel_residual.mean()
    return q,{"valid_exact_names":int(r.name.nunique()),"normalized_name_collisions":len(bad),
              "unmatched_record_names_or_collision":len(records)-len(r)}

def center(x,groups,w):
    out=np.empty_like(x,dtype=float)
    for idx in groups:
        out[idx]=x[idx]-np.average(x[idx],weights=w[idx])
    return out

def setup(q,traits):
    t=q.loc[q.name.isin(traits)].copy()
    if t.empty:return None
    # Fixed cell selection only by number of matched species, no outcome filtering.
    n=t.groupby(["study","year","turnover_class"]).name.transform("size")
    t=t.loc[n>=4].reset_index(drop=True)
    if t.empty:return None
    species=sorted(t.name.unique())
    labels=np.array([traits[x] for x in species],float)
    ix={x:i for i,x in enumerate(species)}
    codes=t.name.map(ix).to_numpy()
    groups=[np.array(a,dtype=int) for a in t.groupby(["study","year","turnover_class"],sort=True).indices.values()]
    studies=sorted(t.study.unique())
    # Equal cell weights within year, equal year weights within study,
    # then equal species weights within cell. Each study sums to one.
    w=np.empty(len(t),float)
    for idx in groups:
        study=t.study.iloc[idx[0]];year=t.year.iloc[idx[0]]
        ny=t.loc[t.study.eq(study),"year"].nunique()
        nc=t.loc[t.study.eq(study)&t.year.eq(year),"turnover_class"].nunique()
        w[idx]=1/(ny*nc*len(idx))
    sc=t.study.map({x:i for i,x in enumerate(studies)}).to_numpy()
    yc=center(t.channel_residual.to_numpy(float),groups,w)
    gc=collections.defaultdict(list)
    for i,k in enumerate(species):gc[k.split()[0]].append(i)
    movable=[np.array(a,dtype=int) for a in gc.values() if len(a)>1 and np.ptp(labels[a])>1e-12]
    return dict(t=t,labels=labels,codes=codes,groups=groups,studies=studies,w=w,sc=sc,yc=yc,
                species=species,movable=movable)

def crossproducts(z,labels):
    x=center(labels[z["codes"]],z["groups"],z["w"])
    n=len(z["studies"])
    a=np.bincount(z["sc"],weights=z["w"]*x*z["yc"],minlength=n)
    b=np.bincount(z["sc"],weights=z["w"]*x*x,minlength=n)
    return a,b

def effect(a,b):
    d=float(np.sum(b))
    return float(np.sum(a)/d) if d>1e-12 else None

def fit(q,traits,label,draws=N_PERM):
    z=setup(q,traits)
    if z is None:return {"status":"no_eligible_trait_cells","p_for_family":1.,"gate":False}
    a,b=crossproducts(z,z["labels"]);obs=effect(a,b)
    n=len(a);informative=int((b>1e-12).sum())
    k=len(z["species"]);mov=z["movable"];nm=sum(map(len,mov))
    support={"trait_species_before_cell_gate":int(q.loc[q.name.isin(traits),"name"].nunique()),
       "species":k,"studies":n,"informative_studies":informative,"cells":len(z["groups"]),
       "genus_count":len({s.split()[0] for s in z["species"]}),
       "movable_genera":len(mov),"movable_species":nm,"plant_study_year_class_rows":len(z["t"]),
       "cell_class_counts":{str(k):int(v) for k,v in z["t"].groupby(["study","year","turnover_class"]).size().reset_index().turnover_class.value_counts().items()}}
    gate=informative>=MIN_STUDIES and k>=MIN_SPECIES and nm>=MIN_MOVABLE and len(mov)>=MIN_GENERA and obs is not None
    out={"support":support,"gate":bool(gate),"coefficient":obs,"p_for_family":1.}
    if label=="symmetry":
        out["species_state_counts"]={"radial":int(np.sum(z["labels"]==0)),"bilateral":int(np.sum(z["labels"]==1))}
        gate=gate and min(out["species_state_counts"].values())>=15
        out["gate"]=bool(gate)
    if obs is None:return out
    # Study cluster bootstrap on fixed carrier and cell design.
    rng=np.random.default_rng(SEED)
    ix=rng.integers(n,size=(9999,n))
    aa=a[ix].sum(1);bb=b[ix].sum(1)
    good=bb>1e-12;bo=aa[good]/bb[good]
    loo=[effect(np.delete(a,i),np.delete(b,i)) for i in range(n)]
    loo=[v for v in loo if v is not None]
    out["study_bootstrap_ci"]=np.quantile(bo,[.025,.975]).tolist() if len(bo)>=9500 else None
    out["leave_study_range"]=[min(loo),max(loo)] if len(loo)==n and n>=2 else None
    # Coherent genus-restricted permutation: a species keeps one trait
    # value across all studies. Cell composition and weights never reselected.
    if not mov:return out
    vals=[];rng=np.random.default_rng(SEED+1)
    for _ in range(draws):
        x=z["labels"].copy()
        for idx in mov:x[idx]=rng.permutation(x[idx])
        ca,cb=crossproducts(z,x);v=effect(ca,cb)
        if v is not None:vals.append(v)
    if len(vals)<.95*draws:return out
    v=np.array(vals);mu=float(v.mean())
    p=min(1.,2*min((1+int(np.sum(v>=obs-1e-12)))/(len(v)+1),(1+int(np.sum(v<=obs+1e-12)))/(len(v)+1)))
    out.update(genus_null_mean=mu,genus_null_interval=np.quantile(v,[.025,.975]).tolist(),genus_residual=obs-mu,
               genus_permutation_p=p,valid_permutations=len(v),p_for_family=p if gate else 1.)
    # Conditional residual stability: subtract the same coherent-null mean
    # from each bootstrap and leave-study result; not a joint null bootstrap.
    return out

def holm(p):
    order=sorted(range(len(p)),key=lambda i:p[i]);r=[0.]*len(p);prev=0.
    for j,i in enumerate(order):
        prev=max(prev,min(1.,(len(p)-j)*p[i]));r[i]=prev
    return r

def finalize(out):
    p=[out[k]["p_for_family"] if "p_for_family" in out[k] else 1. for k in PRIMARY]
    for k,v in zip(PRIMARY,holm(p)):
        r=out[k];r["holm_p"]=v;ci=r.get("study_bootstrap_ci");loo=r.get("leave_study_range")
        # A positive architecture slope plus genus-restricted two-sided test
        # and sign stability on raw effect. Also require obs above null mean.
        r["primary_supported"]=bool(r.get("gate") and v<.05 and ci and loo and ci[0]>0 and loo[0]>0 and r.get("genus_residual",0)>0)
    return out

def self_test():
    rows=[]
    for i in range(12):
        for j in range(6):
            for status in ("one_on_either_side","both_at_least_two"):
                rows.append({"study":str(i),"year":2020,"turnover_class":status,"name":f"genus{j//2} sp{j}","channel_residual":.1*(j%2)+.2*i+(.3 if status=="one_on_either_side" else 0)})
    q=pd.DataFrame(rows);t={x:float(int(x[-1])%2) for x in q.name.unique()}
    z=setup(q,t);a,b=crossproducts(z,z["labels"])
    assert abs(effect(a,b)-.1)<1e-12
    assert np.allclose(np.bincount(z["sc"],weights=z["w"]),1)
    r=fit(q,t,"symmetry",99);assert abs(r["coefficient"]-.1)<1e-12 and not r["gate"] and r["p_for_family"]==1.
    assert holm([.04,.01])==[.04,.02]
    print("FD_EXTERNAL_ARCHITECTURE_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true")
    a=p.parse_args();self_test()
    if a.self_test:return
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);pins={}
        for f,(sha,blob) in SOURCES.items():pins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
        d=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False)
        assert len(d)==623476
        meta=pd.read_csv(w/"Metadata.csv",dtype=str);records,census=collect(prepare(d,meta))
        original=summarize(records)
        assert len(records)==3043 and abs(original["study_equal"]["channel_residual"]-.039527962549189964)<1e-10
        q,qa=stratum_table(records)
        _,source,nhm,wang=source_tables(include_values=True)
        symmetry={k:float(v=="bilateral") for k,v in wang["unannotated_without_four_families"].items() if v in ("radial","bilateral")}
        symmetry_all={k:float(v=="bilateral") for k,v in wang["all_publisher_binary"].items() if v in ("radial","bilateral")}
        ll=nhm["corolla_tube_range_cm"]
        midpoint={k:(math.log2(lo)+math.log2(hi))/2 for k,(lo,hi) in ll.items()}
        maximum={k:math.log2(hi) for k,(lo,hi) in ll.items()}
        main=finalize({"symmetry":fit(q,symmetry,"symmetry"),"tube_midpoint":fit(q,midpoint,"tube_midpoint")})
        sens={"symmetry_all_publisher":fit(q,symmetry_all,"symmetry",199),
              "tube_max":fit(q,maximum,"tube_max",199),
              "symmetry_both_at_least_two":fit(q.loc[q.turnover_class.eq("both_at_least_two")],symmetry,"symmetry",199),
              "tube_both_at_least_two":fit(q.loc[q.turnover_class.eq("both_at_least_two")],midpoint,"tube_midpoint",199)}
        result={"source_pins":pins,"plant_source_pins":source["source_pins"],"parent_replay_residual":original["study_equal"]["channel_residual"],"parent_census":census,
                "name_audit":qa,"primary":main,"descriptive_sensitivity":sens,
                "scope":"Post-discovery floral moderation, external-only. Trait effects conditional on observed turnover and carrier support. Fixed study-year-turnover-class cells; equal cell/year/study weights, within-cell regression and coherent genus trait shuffles; study bootstrap conditional on observed species. No pollination effectiveness, evolutionary or island effect claim. Publisher symmetry not direct measurement; tube length not accessibility."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False))
        print("FD_EXTERNAL_ARCHITECTURE_RESULT",json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":main()
