#!/usr/bin/env python3
"""Coherent visitor-genus trait permutation after exact family-mix accounting."""
import argparse, collections, json, pathlib, tempfile
import numpy as np
import pandas as pd
from external_eupolltrait_size_turnover import load_traits
from external_eupolltrait_size_qc_family import corrected_traits,network,conditional,MIN_STRATUM_DRAWS
from external_eupolltrait_intake import prepare_network
from external_euppollnet_bee_loss import BASE,SOURCES,download
from external_euppollnet_spatial_scaling import assign_genus,infer,fixed_weights

SEED=20261128
DRAWS=999

def cross_codes(a,b,F,K):
    i=np.array(sorted(a-b),dtype=int);j=np.array(sorted(b-a),dtype=int)
    x=np.minimum(i[:,None],j[None,:]).ravel()
    y=np.maximum(i[:,None],j[None,:]).ravel()
    n=len(F);species=x*n+y
    family=np.minimum(F[x],F[y])*K+np.maximum(F[x],F[y])
    return species,family

def sparse_weights(a,b,sa,sb,T,F,K):
    z=conditional(a,b,sa,sb,T,F,K)
    if z is None or not z.get("eligible") or not z.get("family_standardization_eligible"):
        return z
    obs,oh=cross_codes(a,b,F,K);N=len(obs);n=len(T)
    p=np.bincount(oh,minlength=K*K)/N
    C=np.zeros(K*K,dtype=np.int64);hits=C.copy();keys=[]
    for x,y in zip(sa,sb):
        if len(x&y)!=len(a&b):continue
        k,h=cross_codes(x,y,F,K)
        assert len(k)==N
        counts=np.bincount(h,minlength=K*K)
        C+=counts;hits+=counts>0;keys.append(k)
    assert len(keys)==z["accepted"] and int(C.sum())==len(keys)*N
    assert np.all(hits[p>0]>=MIN_STRATUM_DRAWS)
    k,c=np.unique(np.concatenate(keys),return_counts=True)
    i,j=k//n,k%n
    h=np.minimum(F[i],F[j])*K+np.maximum(F[i],F[j])
    w=collections.defaultdict(float)
    for kk,v in zip(k,c*p[h]/C[h]):w[int(kk)]+=float(v)
    for kk in obs:w[int(kk)]-=1/N
    w={kk:v for kk,v in w.items() if abs(v)>1e-15}
    assert abs(sum(w.values()))<1e-10
    # Each exact family-pair stratum separately has zero signed mass.
    check=np.zeros(K*K)
    for kk,v in w.items():
        i,j=divmod(kk,n);check[min(F[i],F[j])*K+max(F[i],F[j])]+=v
    assert np.max(np.abs(check))<1e-10
    actual=sum(v*abs(T[k//n]-T[k%n]) for k,v in w.items())
    assert abs(actual-z["family_standardized_size_gain"])<1e-10
    return dict(z,weights=w)

def genus_design(names,families,T):
    groups=collections.defaultdict(list)
    for i,s in enumerate(names):groups[(families[s],s.split()[0])].append(i)
    groups=[np.array(v,dtype=int) for k,v in sorted(groups.items())]
    code=np.empty(len(T),dtype=int)
    for j,g in enumerate(groups):code[g]=j
    E=np.zeros((len(groups),len(groups)))
    for i,a in enumerate(groups):
        for j in range(i,len(groups)):
            b=groups[j];d=np.abs(T[a,None]-T[b][None,:])
            if i==j:
                mu=float(d.sum()/(len(a)*(len(a)-1))) if len(a)>1 else 0.
            else:mu=float(d.mean())
            E[i,j]=E[j,i]=mu
    movable=[g for g in groups if len(g)>1 and np.ptp(T[g])>1e-12]
    return groups,code,E,movable

def panel(rows,T,names,families):
    if rows.empty:return {"status":"no_eligible_rows","supported":False}
    r=rows.loc[rows.family_standardization_eligible].copy()
    r,audit=assign_genus(r)
    if r.empty:return {"status":"no_family_common_support","supported":False}
    n=len(T);groups,codes,E,movable=genus_design(names,families,T)
    r["visitor_genus_null_expectation"]=[sum(w*E[codes[k//n],codes[k%n]] for k,w in z.items()) for z in r.weights]
    r["visitor_genus_adjusted_gain"]=r.family_standardized_size_gain-r.visitor_genus_null_expectation
    keys=["study","year","plant","genus"]
    cols=["family_standardized_size_gain","visitor_genus_null_expectation","visitor_genus_adjusted_gain"]
    q=r.groupby(keys,as_index=False)[cols].mean()
    q["w"]=fixed_weights(q)
    ww={tuple(v[:4]):v[4] for v in q[keys+["w"]].itertuples(index=False,name=None)}
    counts=r.groupby(keys).size().to_dict()
    W=np.zeros(n*n)
    for row in r.itertuples():
        k=(row.study,row.year,row.plant,row.genus);scale=ww[k]/counts[k]/q.w.sum()
        for j,v in row.weights.items():W[j]+=v*scale
    inds=np.flatnonzero(abs(W)>1e-15);ii=inds//n;jj=inds%n;w=W[inds]
    obs=float(w@np.abs(T[ii]-T[jj]));mu=float(w@E[codes[ii],codes[jj]])
    assert abs(obs-np.average(q.family_standardized_size_gain,weights=q.w))<1e-10
    assert abs(mu-np.average(q.visitor_genus_null_expectation,weights=q.w))<1e-10
    active=set(ii)|set(jj);mov=[g for g in movable if any(x in active for x in g)]
    nm=sum(sum(x in active for x in g) for g in mov)
    rng=np.random.default_rng(SEED+51)
    perm=[]
    for _ in range(DRAWS):
        v=T.copy()
        for g in mov:v[g]=rng.permutation(v[g])
        perm.append(float(w@np.abs(v[ii]-v[jj])))
    z=np.array(perm)
    p=min(1.,2*min((1+int(np.sum(z>=obs-1e-12)))/(DRAWS+1),(1+int(np.sum(z<=obs+1e-12)))/(DRAWS+1)))
    ci=infer(q,"visitor_genus_adjusted_gain")
    pairs=int(r[["study","year","site_a","site_b"]].drop_duplicates().shape[0])
    gate=ci["gate"] and pairs>=30 and len(mov)>=10 and nm>=30
    return {"name_audit":audit,"comparisons":len(r),"site_pairs":pairs,
        "visitor_species":n,"visitor_genera_with_family_stratum":len(groups),"active_visitor_species":len(active),
        "movable_visitor_genera":len(mov),"active_movable_visitor_species":int(nm),
        "family_standardized_gain":obs,"visitor_genus_label_null_exact_mean":mu,
        "visitor_genus_label_null_mc_mean":float(z.mean()),"visitor_genus_label_null_interval":np.quantile(z,[.025,.975]).tolist(),
        "visitor_genus_adjusted_gain":obs-mu,"permutation_p":p,"draws":DRAWS,
        "crossed_study_plant_genus_residual":ci,"support_gate":bool(gate),
        "supported":bool(gate and p<.05 and obs>mu and ci["supported"])}

def self_test():
    names=["a aa","a ab","b ba","b bb"];fam={s:"f" for s in names};T=np.array([0.,1.,3.,4.])
    g,c,E,m=genus_design(names,fam,T)
    assert len(m)==2 and E[0,0]==1 and E[0,1]==3
    # Same degrees, overlap and family mass; coefficient replay.
    a={0,1};b={2,3};sa=[{0,2}]*40;sb=[{1,3}]*40
    z=sparse_weights(a,b,sa,sb,T,np.zeros(4,dtype=int),1)
    assert z["family_standardization_eligible"]
    assert abs(sum(z["weights"].values()))<1e-12
    assert abs(z["family_standardized_size_gain"]+1)<1e-12
    print("FD_EUPOLLTRAIT_VISITOR_GENUS_SELF_TEST_PASS",flush=True)

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
        rows,diag=network(d,tab,"bee",tax,pair_fn=sparse_weights)
        assert len(names)==484 and len(rows)==1261
        main=panel(rows,T,names,fm)
        sens=panel(rows.loc[rows.min_exclusive.ge(2)],T,names,fm)
        assert main["comparisons"]==1117 and abs(main["family_standardized_gain"]-.04219365065251101)<1e-10
        assert sens["comparisons"]==314 and abs(sens["family_standardized_gain"]-.055920059187065774)<1e-10
        out={"source_pins":pins,"network_source_pins":npins,"selection":audit,"network":diag,
          "primary":main,"descriptive_both_exclusive_ge2":sens,
          "scope":"Post-discovery visitor-genus trait label audit. Exact parent pair support and family-standardized weights replayed; exchangeability only within visitor family/genus, one size label per species across all sites. Crossed bootstrap by study and focal plant genus, not a calibrated phylogenetic or visitor-degree model. No independent observations, actual pollination, or causal interpretation."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
        print("FD_EUPOLLTRAIT_VISITOR_GENUS_RESULT",json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()
