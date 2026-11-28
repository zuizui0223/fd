#!/usr/bin/env python3
"""External three-dimensional bee morphology, common-carrier fixed-margin turnover."""
import argparse, collections, json, math, pathlib, tempfile
import numpy as np
import pandas as pd
from external_eupolltrait_size_turnover import load_traits,parse_range
from external_eupolltrait_size_qc_family import corrected_traits,network
from external_eupolltrait_intake import prepare_network,norm
from external_euppollnet_bee_loss import BASE,SOURCES,download
from external_euppollnet_spatial_scaling import assign_genus,infer,fixed_weights
from external_eupolltrait_energy import components
from external_eupolltrait_visitor_genus import genus_design

AXES=("itd","hair_length","hair_density")
SEED=20261128
MIN_ACCEPT=30

def hair_table(m,kind):
    x=m.loc[m.group.eq("bee")&m.measurementType.eq(kind)].copy()
    assert len(x)>0
    census={k:{str(a):int(b) for a,b in x[k].value_counts(dropna=False).items()} for k in ("sex","originalColumn","lifeStage","measurementUnit")}
    census["source_nonempty"]=int(x.measurementSource.str.strip().ne("").sum())
    # Same source's female category; require explicit F and a single
    # scalar, no inference of sex from taxonomic identity. NA unit is
    # documented source export convention; uniform raw unit per trait.
    f=x.loc[x.sex.str.strip().str.casefold().eq("f")].copy()
    if f.empty:return {},{"raw_rows":len(x),"schema":census,"status":"no_exact_f"}
    unit=set(f.measurementUnit.str.strip().str.casefold())
    known={"na","nan","","mm","µm","μm","um"} if kind=="hair_length" else {"na","nan","","hairs/0.1 mm²","hairs/0.1 mm2"}
    if len(unit)!=1 or not unit.issubset(known):
        return {},{"raw_rows":len(x),"schema":census,"status":"unresolved_or_mixed_units"}
    bad=collections.Counter();tab={}
    for k,q in f.groupby("name",sort=True):
        vals=set()
        for row in q.itertuples():
            if row.lifeStage.strip().casefold() not in {"adult","adults","imago","na","nan","none",""}:
                bad["stage"]+=1;continue
            v=parse_range(row.measurementValue)
            if v is None or v[0]!=v[1]:
                bad["nonpositive_or_nonscalar"]+=1;continue
            vals.add(v[0])
        if len(vals)!=1:
            bad["no_unique_value"]+=1;continue
        tab[k]=math.log2(next(iter(vals)))
    return tab,{"raw_rows":len(x),"exact_f_rows":len(f),"schema":census,"exclusions":dict(bad),"retained_species":len(tab),
                "raw_positive_range":[2**min(tab.values()),2**max(tab.values())] if tab else None}

def matrix(tables):
    keys=sorted(set.intersection(*(set(t) for t in tables)))
    if not keys:return {},None,{}
    X=np.array([[t[k] for t in tables] for k in keys])
    sd=X.std(0);mu=X.mean(0)
    if len(keys)<2 or (sd<=1e-12).any():return {},None,{"status":"constant_or_insufficient_common_carriers","species":len(keys)}
    Z=(X-mu)/sd
    return {k:Z[i] for i,k in enumerate(keys)},Z,{"common_source_species":len(keys),"log2_means":mu.tolist(),"log2_sd":sd.tolist(),
                    "source_axis_correlation":np.corrcoef(Z.T).tolist()}

def metric_bank(Z):
    dif=Z[:,None,:]-Z[None,:,:]
    # Equal standardized axes, root mean squared Euclidean distance.
    full=np.sqrt(np.mean(dif*dif,axis=2))
    out={"joint":full}
    for j,a in enumerate(AXES):out[a]=np.abs(dif[:,:,j])
    out["hair_joint"]=np.sqrt(np.mean(dif[:,:,1:]**2,axis=2))
    return out

def pair(a,b,sa,sb,bank):
    u,v=a-b,b-a
    if not u or not v:return None
    obs={k:components(u,v,D,True) for k,D in bank.items()}
    sims={k:[] for k in bank};halves=[0,0]
    for i,(x,y) in enumerate(zip(sa,sb)):
        assert len(x)==len(a) and len(y)==len(b)
        if len(x&y)!=len(a&b):continue
        for k,D in bank.items():sims[k].append(components(x-y,y-x,D,True))
        halves[int(i>=len(sa)//2)]+=1
    n=sum(halves)
    if n<MIN_ACCEPT or min(halves)<10:return {"accepted":n,"eligible":False}
    out={"accepted":n,"eligible":True}
    for k,x in obs.items():
        y=np.mean(sims[k],axis=0);g=y-x
        assert abs(g[0]-g[1]-g[2])<1e-10
        for j,name in enumerate(("cross","dispersion","energy")):
            for m,val in (("obs",x[j]),("null",y[j]),("gain",g[j])):out[f"{k}_{name}_{m}"]=float(val)
    return out

def panel(rows):
    if rows.empty:return {"status":"no_eligible_pairs","supported":False}
    r,audit=assign_genus(rows)
    if r.empty:return {"status":"no_exact_names","supported":False}
    cols=[c for c in r if c.endswith(("_obs","_null","_gain"))]
    q=r.groupby(["study","year","plant","genus"],as_index=False)[cols].mean();w=fixed_weights(q)
    stats={k:float(np.average(q[k],weights=w)) for k in cols}
    pairs=int(r[["study","year","site_a","site_b"]].drop_duplicates().shape[0])
    tests={k:infer(q,k) for k in ("joint_energy_gain","joint_cross_gain","joint_dispersion_gain",
                                  "itd_energy_gain","hair_length_energy_gain","hair_density_energy_gain","hair_joint_energy_gain")}
    return {"name_audit":audit,"comparisons":len(r),"site_pairs":pairs,"study_equal":stats,
            "inference":tests,"supported":bool(pairs>=30 and tests["joint_energy_gain"]["supported"]),
            "secondary_metrics_descriptive_only":True}

def self_test():
    T={"a":1.,"b":2.,"c":3.};u={"a":0.,"b":1.,"c":2.};v={"a":2.,"b":0.,"c":1.}
    z,X,a=matrix([T,u,v]);assert len(z)==3 and np.allclose(X.std(0),1)
    bank=metric_bank(X);assert np.allclose(bank["joint"].diagonal(),0)
    out=pair({0},{1},[{0}]*40,[{2}]*40,bank)
    assert out["eligible"] and out["joint_dispersion_gain"]==0
    print("FD_EUPOLLTRAIT_MULTIVARIATE_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    tax,m,meta,pins=load_traits()
    # Inspect full definitions; no derived hairiness/tongue index.
    mask=meta.trait.astype(str).map(norm).str.contains("hairlength|hairdensity|itd")
    print("FD_EUPOLLTRAIT_MULTIVARIATE_METADATA",json.dumps(meta.loc[mask].to_dict("records"),allow_nan=False),flush=True)
    itd,ia=corrected_traits(m,"bee");hl,la=hair_table(m,"hair_length");hd,da=hair_table(m,"hair_density")
    z,Z,ma=matrix([itd,hl,hd])
    audit={"itd":ia,"hair_length":la,"hair_density":da,"scaling":ma}
    print("FD_EUPOLLTRAIT_MULTIVARIATE_CARRIERS",json.dumps(audit,allow_nan=False),flush=True)
    if not z:
        out={"source_pins":pins,"selection":audit,"status":"no_resolved_common_trait_schema","primary":{"supported":False}}
    else:
        with tempfile.TemporaryDirectory() as td:
            w=pathlib.Path(td);npins={}
            for f,(sha,blob) in SOURCES.items():npins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
            raw=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False);assert len(raw)==623476
            d,na=prepare_network(raw,pd.read_csv(w/"Metadata.csv",dtype=str))
            fm=tax.set_index("name").family.to_dict()
            q=d.loc[d.group.eq("bee")&d.name.isin(z)&d.name.map(fm).eq(d.family)]
            names=sorted(q.name.unique())
            X=np.array([z[k] for k in names])
            bank=metric_bank(X)
            tab={k:float(z[k][0]) for k in z}
            def pair_fn(a,b,sa,sb,T,F,K):
                assert np.array_equal(T,X[:,0])
                return pair(a,b,sa,sb,bank)
            rows,diag=network(d,tab,"bee",tax,pair_fn=pair_fn)
            out={"source_pins":pins,"network_source_pins":npins,"selection":audit,"network_name_audit":na,
                 "network":diag,"axes":list(AXES),"primary":panel(rows),
                 "descriptive_both_exclusive_ge2":panel(rows.loc[rows.min_exclusive.ge(2)] if not rows.empty else rows),
                 "scope":"Post-discovery common measured trait carrier subweb, adult morphology only. Joint standardized Euclidean half-energy, fixed local degrees/shared counts, crossed study and plant genus uncertainty. No independent interaction replication, visitor-genus-adjusted claim, actual pollen transport or island effect. Correlated axes are not independent functions."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("FD_EUPOLLTRAIT_MULTIVARIATE_RESULT",json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()
