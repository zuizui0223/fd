#!/usr/bin/env python3
"""Post-discovery EuPPollNet/EuPollTrait tests of genus and plant morphology.

Replay the exact original eligibility/null first. Species-coherent within-genus
size exchange on fixed observed links. Independent NHM recorded tube lengths
against mean visitor body size on fixed-margin local measured subwebs.
No project botanical inputs or inferred historical island response.
"""
import argparse
import collections
import hashlib
import json
import pathlib
import numpy as np
from scipy import sparse
from numba import njit
from external_bee_itd_network_null import traits, source_visits, swaps, breadth, summary, SEED, DRAWS
from tropical_pollinator_source_preflight import get,blobsha
from external_wang_nhm_trait_bridge import NHM_ROOT,NHM,strict_nhm

@njit
def pairweights(A):
    W=np.zeros((A.shape[1],A.shape[1]))
    n=0
    for row in A:
        ii=np.where(row)[0]
        k=len(ii)
        if k<2:
            continue
        n+=1
        w=2./(k*(k-1))
        for i in range(k):
            for j in range(i):
                a,b=ii[i],ii[j]
                W[max(a,b),min(a,b)]+=w
    return W/n

def collect(d,t):
    out=[];reject=collections.Counter()
    for key,df in d.groupby(["Study_id","Network_id","year"],sort=True):
        loc=df[["lat","lon"]].drop_duplicates()
        if len(loc)!=1:
            reject["multi_coord_site_year"]+=1;continue
        df=df[df.visitor.isin(t)]
        if df.empty:
            reject["no_trait_visitor"]+=1;continue
        P=sorted(set(df.plant));B=sorted(set(df.visitor))
        if len(P)<5 or len(B)<5:
            reject["small_subweb"]+=1;continue
        A=np.zeros((len(P),len(B)),dtype=bool)
        pi={x:i for i,x in enumerate(P)};bi={x:i for i,x in enumerate(B)}
        A[[pi[x] for x in df.plant],[bi[x] for x in df.visitor]]=True
        vals=np.array([t[x] for x in B])
        if (A.sum(1)>=2).sum()<3 or len(set(vals))<3:
            reject["insufficient_breadth"]+=1;continue
        orig=A.copy();E=np.argwhere(A).astype(np.int64);n=len(E)
        seed=(SEED+int(hashlib.sha256("|".join(map(str,key)).encode()).hexdigest()[:8],16))%(2**31-1)
        if swaps(E,A,20*n,seed)<max(10,n):
            reject["burnin_low_switch"]+=1;continue
        draws=[];moves=[]
        for i in range(DRAWS):
            moves.append(swaps(E,A,5*n,(seed+i+1)%(2**31-1)))
            draws.append(breadth(A,vals))
        assert np.array_equal(orig.sum(0),A.sum(0)) and np.array_equal(orig.sum(1),A.sum(1))
        if np.mean(moves)<max(2,.1*n) or np.std(draws)<1e-12:
            reject["low_switch_or_degenerate_null"]+=1;continue
        lat,lon=loc.iloc[0]
        out.append({"study":str(key[0]),"block":(int(np.floor(lat/10)),int(np.floor(lon/10))),"key":key,"plants":P,"bees":B,"A":orig,"observed":float(breadth(orig,vals)),"draw":np.array(draws)})
    return out,dict(reject)

def permutations(names,values,draws=999):
    idx=collections.defaultdict(list)
    for i,x in enumerate(names):idx[x.split()[0]].append(i)
    varying=[np.array(v) for v in idx.values() if len(v)>1 and len(set(values[v]))>1]
    rng=np.random.default_rng(SEED+1)
    res=np.tile(values,(draws+1,1))
    for k in range(1,draws+1):
        for v in varying:res[k,v]=values[rng.permutation(v)]
    return res,{"movable_genera":len(varying),"movable_species":sum(map(len,varying)),"species":len(names)}

def coherent_breadth(rr,t):
    ids=sorted(set().union(*(set(z["bees"]) for z in rr)))
    ix={s:i for i,s in enumerate(ids)};n=len(ids)
    values=np.array([t[s] for s in ids])
    sims,census=permutations(ids,values)
    # Use sparse, block/study-equal weights to avoid recalculating local
    # pair loops for every coherent species draw. Each observed plant with
    # >=2 visitors gets equal weight within its network.
    pairs=[];weights=[]
    for z in rr:
        W=pairweights(z["A"]);a,b=np.nonzero(W)
        q=np.array([ix[x] for x in z["bees"]])
        lo=np.minimum(q[a],q[b]);hi=np.maximum(q[a],q[b])
        pairs.append(hi*n+lo);weights.append(W[a,b])
    pair_ids=np.unique(np.concatenate(pairs));row=[];col=[];w=[]
    for i,(p,v) in enumerate(zip(pairs,weights)):
        row.extend([i]*len(p));col.extend(np.searchsorted(pair_ids,p));w.extend(v)
    A=sparse.csr_matrix((w,(row,col)),shape=(len(rr),len(pair_ids)))
    assert np.allclose(np.asarray(A.sum(1)).ravel(),1)
    # A sparse aggregation of within-network coefficients, one group at
    # a time; trait draws coherent over repeated species across all groups.
    a,b=pair_ids//n,pair_ids%n
    observed=np.array([z["observed"] for z in rr])
    assert np.allclose(A@np.abs(values[a]-values[b]),observed)
    output={}
    for k in ("block","study"):
        groups=sorted({z[k] for z in rr});g=np.array([groups.index(z[k]) for z in rr])
        counts=np.bincount(g)
        G=sparse.csr_matrix((1/counts[g],(g,np.arange(len(rr)))),shape=(len(groups),len(rr)))@A
        X=np.stack([G@np.abs(v[a]-v[b]) for v in sims],axis=1)
        fake=[{"observed":x[0],"draw":x[1:],"block":i,"study":i} for i,x in enumerate(X)]
        r=summary(fake,k)
        r.update(networks=len(rr),groups=len(groups),studies=len({z["study"] for z in rr}),support_gate=len(rr)>=20 and len(groups)>=10 and census["movable_species"]>=30 and census["movable_genera"]>=10)
        output[k]=r
    return {"permuted_species":census,"panels":output}

@njit
def slope(A,x,v):
    # plant-equal regression of mean log2 visitor size on log2 tube length.
    y=(A.astype(np.float64)@v)/A.sum(1)
    return np.dot(x,y)/np.dot(x,x)

def plant_matching(rr,t,length,mode):
    out=[];rejection=collections.Counter();pp=set();bb=set()
    for z in rr:
        keep=np.array([s.casefold() in length for s in z["plants"]])
        if keep.sum()<5:
            rejection["fewer_than_5_documented_plants"]+=1;continue
        A=z["A"][keep].copy();cols=A.sum(0)>0;A=A[:,cols]
        if A.shape[1]<5:
            rejection["fewer_than_5_documented_bees"]+=1;continue
        P=np.array(z["plants"])[keep];B=np.array(z["bees"])[cols]
        lengths=np.array([length[s.casefold()] for s in P])
        x=np.log2(lengths[:,1]) if mode=="max" else np.log2(lengths).mean(1)
        x-=x.mean()
        if np.sqrt(np.mean(x*x))<.25:
            rejection["low_flower_length_variation"]+=1;continue
        v=np.array([t[s] for s in B])
        if len(set(v))<3:
            rejection["low_bee_size_variation"]+=1;continue
        obs=slope(A,x,v);E=np.argwhere(A).astype(np.int64);n=len(E)
        d0=A.sum(0).copy();d1=A.sum(1).copy()
        seed=(SEED+3+int(hashlib.sha256("|".join(map(str,z["key"])).encode()).hexdigest()[:8],16))%(2**31-1)
        if swaps(E,A,20*n,seed)<max(10,n):
            rejection["low_burnin_switch"]+=1;continue
        sims=[];moves=[]
        for i in range(DRAWS):
            moves.append(swaps(E,A,5*n,(seed+i+1)%(2**31-1)))
            sims.append(slope(A,x,v))
        assert np.array_equal(A.sum(0),d0) and np.array_equal(A.sum(1),d1)
        if np.mean(moves)<max(2,.1*n) or np.std(sims)<1e-12:
            rejection["low_switch_or_degenerate_null"]+=1;continue
        pp.update(P);bb.update(B)
        out.append({"observed":float(obs),"draw":np.array(sims),"block":z["block"],"study":z["study"]})
    return {"panels":{k:summary(out,k) for k in ("block","study")},"distinct_plants":len(pp),"distinct_bees":len(bb),"rejections":dict(rejection)}

def self_test():
    A=np.array([[1,1,0],[0,1,1],[1,0,1]],dtype=bool);v=np.array([0.,1.,3.])
    assert np.isclose(np.sum(pairweights(A)*np.abs(v[:,None]-v[None,:])),breadth(A,v))
    x=np.array([-1.,0.,1.])
    assert np.isclose(slope(A,x,v),.5)
    names=["A a","A b","B c"];r,c=permutations(names,v,5)
    assert c["movable_species"]==2 and np.all(r[:,2]==3)
    assert all(sorted(x[:2])==[0.,1.] for x in r)
    print("FD_GENUS_MORPHOLOGY_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true")
    a=p.parse_args();self_test()
    if a.self_test:return
    if a.out is None:p.error("--out required")
    t,ti=traits();assert t
    d,di=source_visits();rr,reject=collect(d,t)
    replay=summary(rr,"block")
    assert len(rr)==980 and np.isclose(replay["observed"],.47335833402980687,atol=1e-10,rtol=0)
    assert np.isclose(replay["null_mean"],.5890613488588909,atol=1e-10,rtol=0),"original null replay differs"
    tax=coherent_breadth(rr,t)
    print("FD_GENUS_COHERENT_RESULT",json.dumps(tax,allow_nan=False),flush=True)
    raw=get(NHM_ROOT+"combined.csv",maxbytes=25_000_000)
    assert blobsha(raw)==NHM["combined.csv"]
    table,info=strict_nhm(raw);length=table["corolla_tube_range_cm"]
    match={k:plant_matching(rr,t,length,k) for k in ("midpoint","max")}
    print("FD_PLANT_TUBE_MATCHING_RESULT",json.dumps(match,allow_nan=False),flush=True)
    prim=[tax["panels"]["block"],match["midpoint"]["panels"]["block"]]
    pvals=[z.get("conditional_p",1) if z.get("support_gate",False) else 1 for z in prim]
    order=np.argsort(pvals);adj=[1.,1.];last=0.
    for i,j in enumerate(order):
        last=max(last,min(1.,(2-i)*pvals[j]));adj[j]=last
    result={"source_replay":replay,"rejections":reject,"source_pins":{"bee_archive_sha256":ti["archive_sha256"],"visits":di["source_pins"],"nhm_blob":blobsha(raw),"nhm_sha256":hashlib.sha256(raw).hexdigest()},"genus_coherent":tax,"plant_tube_matching":match,"holm2":{"genus_conditioned_size_breadth":adj[0],"tube_matching_midpoint":adj[1]},"scope":"Post-discovery tests; Europe-only observed bee subwebs. Within-genus coherent species exchange is not a phylogenetic model. Tube length is not depth or measured reach; ITD is not independent tongue length; no causal mechanism or island effect. Study-equal and max-length sensitivities descriptive. Positive matching slope is prediction, two-sided conditional test."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False))
    print("FD_GENUS_MORPHOLOGY_COMPLETE",json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":main()
