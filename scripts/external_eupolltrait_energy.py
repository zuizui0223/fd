#!/usr/bin/env python3
"""Separate cross-partner size distance into distribution divergence and dispersion."""
import argparse, json, pathlib, tempfile
import numpy as np
import pandas as pd
from external_eupolltrait_size_turnover import load_traits
from external_eupolltrait_size_qc_family import corrected_traits,network,MIN_ACCEPT
from external_eupolltrait_visitor_genus import genus_design
from external_eupolltrait_intake import prepare_network
from external_euppollnet_bee_loss import BASE,SOURCES,download
from external_euppollnet_spatial_scaling import assign_genus,infer,fixed_weights

SEED=20261128
METRICS=("cross","dispersion","energy_half")

def components(u,v,D,check_nonnegative=False):
    i=np.array(sorted(u),dtype=int)
    j=np.array(sorted(v),dtype=int)
    assert len(i)>0 and len(j)>0 and not set(i)&set(j)
    cross=float(D[np.ix_(i,j)].mean())
    within=float((D[np.ix_(i,i)].mean()+D[np.ix_(j,j)].mean())/2)
    energy=cross-within
    if check_nonnegative:assert energy>=-1e-10,"half-energy must be nonnegative for absolute metric"
    return np.array([cross,within,energy])

def compare(a,b,sa,sb,D,E):
    u,v=a-b,b-a
    if not u or not v:return None
    observed=components(u,v,D,True)
    genus_observed=components(u,v,E)
    null=[];genus_null=[];halves=[0,0]
    for k,(x,y) in enumerate(zip(sa,sb)):
        assert len(x)==len(a) and len(y)==len(b)
        if len(x&y)!=len(a&b):continue
        xx,yy=x-y,y-x
        assert len(xx)==len(u) and len(yy)==len(v)
        null.append(components(xx,yy,D,True))
        genus_null.append(components(xx,yy,E))
        halves[int(k>=len(sa)//2)]+=1
    n=sum(halves)
    if n<MIN_ACCEPT or min(halves)<10:return {"accepted":n,"eligible":False}
    expected=np.mean(null,axis=0);genus_expected=np.mean(genus_null,axis=0)
    gain=expected-observed;genus_gain=genus_expected-genus_observed
    residual=gain-genus_gain
    for v in (observed,expected,gain,genus_gain,residual):
        assert abs(v[0]-v[1]-v[2])<1e-10
    out={"accepted":n,"eligible":True}
    for k,arr in [("obs",observed),("null",expected),("gain",gain),("genus_expected_gain",genus_gain),("genus_residual",residual)]:
        for m,value in zip(METRICS,arr):out[m+"_"+k]=float(value)
    return out

def panel(rows):
    if rows.empty:return {"status":"no_rows","supported":False}
    r,audit=assign_genus(rows)
    if r.empty:return {"status":"no_exact_focal_names","supported":False}
    keys=["study","year","plant","genus"]
    cols=[m+"_"+v for m in METRICS for v in ("obs","null","gain","genus_expected_gain","genus_residual")]
    q=r.groupby(keys,as_index=False)[cols].mean();w=fixed_weights(q)
    means={k:float(np.average(q[k],weights=w)) for k in cols}
    for v in ("obs","null","gain","genus_expected_gain","genus_residual"):
        assert abs(means["cross_"+v]-means["dispersion_"+v]-means["energy_half_"+v])<1e-10
    tests={k:infer(q,k) for k in ("cross_gain","dispersion_gain","energy_half_gain","energy_half_genus_residual")}
    pairs=int(r[["study","year","site_a","site_b"]].drop_duplicates().shape[0])
    return {"name_audit":audit,"comparisons":len(r),"site_pairs":pairs,"study_equal":means,
            "inference":tests,"supported":bool(pairs>=30 and tests["energy_half_gain"]["supported"]),
            "genus_residual_descriptive_only":True}

def self_test():
    T=np.array([0.,1.,3.,4.]);D=np.abs(T[:,None]-T[None,:])
    a=components({0,1},{2,3},D,True)
    assert np.allclose(a,[3.,.5,2.5])
    b=components({0,2},{1,3},D,True)
    assert np.allclose(b,[2.,1.5,.5])
    # The cross gain can contain a negative dispersion contribution.
    z=compare({0,1},{2,3},[{0,2}]*40,[{1,3}]*40,D,D)
    assert z["eligible"] and np.isclose(z["cross_gain"],-1.)
    assert np.isclose(z["dispersion_gain"],1.) and np.isclose(z["energy_half_gain"],-2.)
    assert np.isclose(z["energy_half_genus_residual"],0.)
    print("FD_EUPOLLTRAIT_ENERGY_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    tax,m,meta,pins=load_traits();tab,audit=corrected_traits(m,"bee")
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);npins={}
        for f,(sha,blob) in SOURCES.items():npins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
        raw=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False);assert len(raw)==623476
        d,na=prepare_network(raw,pd.read_csv(w/"Metadata.csv",dtype=str))
        fm=tax.set_index("name").family.to_dict()
        q=d.loc[d.group.eq("bee")&d.name.isin(tab)&d.name.map(fm).eq(d.family)]
        names=sorted(q.name.unique());T=np.array([tab[k] for k in names])
        D=np.abs(T[:,None]-T[None,:])
        groups,codes,E0,movable=genus_design(names,fm,T)
        E=E0[codes[:,None],codes[None,:]]
        # Identical species has zero distance even after label permutation.
        np.fill_diagonal(E,0.)
        def pair_fn(a,b,sa,sb,tt,F,K):
            assert np.array_equal(T,tt)
            return compare(a,b,sa,sb,D,E)
        rows,diag=network(d,tab,"bee",tax,pair_fn=pair_fn)
        assert len(names)==484 and len(rows)==1261
        main=panel(rows);sens=panel(rows.loc[rows.min_exclusive.ge(2)])
        assert abs(main["study_equal"]["cross_gain"]-.11064886953661242)<1e-10
        assert len(rows.loc[rows.min_exclusive.ge(2)])==333
        out={"source_pins":pins,"network_source_pins":npins,"selection":audit,"network":diag,
             "primary":main,"descriptive_both_exclusive_ge2":sens,
             "scope":"Post-discovery half-energy distribution-distance test on identical corrected bee pairs; exact cross distance = within-exclusive-set dispersion + half-energy. Genus trait-label expected distances descriptive, not causal. Fixed local degrees and shared counts, crossed study/plant genus uncertainty. Not an independent replication or direct pollination measure."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
        print("FD_EUPOLLTRAIT_ENERGY_RESULT",json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()
