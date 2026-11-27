#!/usr/bin/env python3
"""Post-discovery new-target test of the latitude-dependent corolla-length contrast.

Only the metadata-defined targets after the already inspected 288 are downloaded.
Independent target-list extension, not independent trait observations or donors.
"""
import argparse,collections,json,pathlib
import numpy as np
from scipy import sparse
import external_gift_native_flora_roster as roster
from external_gift_floral_replication_coverage import source_tables
from external_gift_continuous_flower_space import groups,REG,REGIONS
SEED=20261128
PREVIOUS_SHA="8945ce01bd1b5a9294fe18ee570ec6ed9d207e93363c70ec3732ffe013f3c935"
PRIMARY="loose__disjoint__midpoint__n5"
SPECS=[("loose","disjoint","midpoint",5),("loose","inclusive","midpoint",5),("loose","disjoint","max",5),("loose","disjoint","midpoint",10),("strict","disjoint","midpoint",5)]

def assemble(rows,tab,matched):
    kept=[];coeff=[];labels=("total","composition","within_radial","within_bilateral") if matched else ("total",)
    names=sorted(set().union(*(z["t"]|z["s"] for z in rows)))
    idx={s:i for i,s in enumerate(names)}
    mats=[([],[],[]) for _ in labels]
    for z in rows:
        a=[sorted(s for s in z["t"] if tab[s][0]==g) for g in (0,1)]
        b=[sorted(s for s in z["s"] if tab[s][0]==g) for g in (0,1)]
        if matched and min(map(len,a+b))<2:continue
        j=len(kept);kept.append(z)
        def add(k,ss,w):
            rr,cc,dd=mats[k];rr.extend([j]*len(ss));cc.extend(idx[s] for s in ss);dd.extend([w]*len(ss))
        add(0,z["t"],1/len(z["t"]));add(0,z["s"],-1/len(z["s"]))
        if not matched:continue
        pt=len(a[1])/len(z["t"]);ps=len(b[1])/len(z["s"]);d=pt-ps;p=(pt+ps)/2
        for g,sign in ((0,-1),(1,1)):
            add(1,a[g],sign*d/(2*len(a[g])));add(1,b[g],sign*d/(2*len(b[g])))
            w=(1-p) if g==0 else p
            add(g+2,a[g],w/len(a[g]));add(g+2,b[g],-w/len(b[g]))
    M=[sparse.csr_matrix((d,(r,c)),shape=(len(kept),len(names))) for r,c,d in mats]
    if matched and len(kept):
        err=M[0]-M[1]-M[2]-M[3]
        assert not err.nnz or max(abs(err.data))<1e-12
    return kept,names,M,labels

def test(rows,tab,matched,perms=999,boots=999):
    if not rows:return {"status":"empty","pairs":0}
    rows,names,M,labels=assemble(rows,tab,matched)
    if not rows:return {"status":"empty","pairs":0}
    x=np.array([tab[s][1] for s in names]);sym=np.array([tab[s][0] for s in names]);g=[s.split()[0] for s in names]
    # Shuffle length values only, within genus AND fixed symmetry.
    cells=collections.defaultdict(list)
    for i,(gen,s) in enumerate(zip(g,sym)):cells[(gen,int(s))].append(i)
    code=np.empty(len(names),int);mov=np.zeros(len(names),bool);expect=np.zeros(len(names))
    for j,ix in enumerate(cells.values()):
        code[ix]=j;expect[ix]=x[ix].mean()
        if len(set(x[ix]))>1:mov[ix]=True
    order=np.argsort(code,kind="stable");idx={s:i for i,s in enumerate(names)}
    support={}
    for r in REGIONS:
        rr=[z for z in rows if z["region"]==r];ids=set().union(*(z["t"]|z["s"] for z in rr))
        mm=[idx[s] for s in ids if mov[idx[s]]]
        support[r]={"pairs":len(rr),"blocks":len({z["block"] for z in rr}),"donors":len({z["donor"] for z in rr}),"movable_species":len(mm),"movable_genera":len({g[i] for i in mm})}
    Y=np.column_stack([m@x for m in M]);E=np.column_stack([m@expect for m in M])
    ds={tag:groups(rows,tag) for tag in ("block","donor")}
    rng=np.random.default_rng(SEED)
    draws={tag:np.empty((perms,4,len(labels))) for tag in ds}
    for j in range(perms):
        u=np.empty_like(x);u[order]=x[np.lexsort((rng.random(len(x)),code))]
        z=np.column_stack([m@u for m in M])
        for tag,(G,agg,N) in ds.items():draws[tag][j]=agg(G@z,np.ones(N))
    out={}
    for tag,(G,agg,N) in ds.items():
        B=G@Y;D=G@E;obs=agg(B,np.ones(N));ex=agg(D,np.ones(N));res=obs-ex
        raw=[];resid=[];leave=[]
        for j in range(boots):
            w=np.bincount(rng.integers(N,size=N),minlength=N)
            raw.append(agg(B,w));resid.append(agg(B-D,w))
        for j in range(N):
            w=np.ones(N);w[j]=0;leave.append(agg(B-D,w))
        raw=np.array(raw);resid=np.array(resid);leave=np.array(leave)
        out[tag]={}
        for i,r in enumerate(REG):
            req=[r] if r in REGIONS else ["north","tropical"]
            if any(support[q]["pairs"]==0 for q in req):continue
            gate=all(support[q]["pairs"]>=50 and support[q]["blocks"]>=30 and support[q]["donors"]>=20 for q in req)
            nullgate=gate and all(support[q]["movable_species"]>=30 and support[q]["movable_genera"]>=10 for q in req)
            out[tag][r]={"geography_gate":bool(gate),"null_gate":bool(nullgate),"metrics":{}}
            for k,m in enumerate(labels):
                v=draws[tag][:,i,k];o=obs[i,k];pv=min(1.,2*min((1+int(np.sum(v<=o+1e-12)))/(perms+1),(1+int(np.sum(v>=o-1e-12)))/(perms+1)))
                def ci(z):
                    v=z[:,i,k];v=v[np.isfinite(v)]
                    return np.quantile(v,[.025,.975]).tolist() if len(v)>=.9*boots else None
                lv=leave[:,i,k]
                out[tag][r]["metrics"][m]={"observed":float(o),"observed_ci":ci(raw),"conditional_null_mean":float(ex[i,k]),"residual":float(res[i,k]),"residual_ci":ci(resid),"coherent_p":float(pv),"residual_leave_range":[float(lv.min()),float(lv.max())] if np.isfinite(lv).all() else None}
    return {"status":"estimated","pairs":len(rows),"species":len(names),"support":support,"group_weightings":out}

def self_test():
    tab={"a a":(0,1.),"a b":(0,2.),"a c":(1,4.),"a d":(1,5.),"a e":(1,6.),"a f":(0,2.),"a g":(0,3.),"a h":(1,5.),"a i":(1,7.)}
    rr=[{"region":"north","block":"b","donor":"d","t":{"a a","a b","a c","a d","a e"},"s":{"a f","a g","a h","a i"}}]
    kept,names,M,labels=assemble(rr,tab,True);x=np.array([tab[s][1] for s in names]);v=[float((m@x)[0]) for m in M]
    assert len(kept)==1 and np.isclose(v[0],sum(v[1:])) and np.isclose(v[0],3.6-4.25)
    print("FD_EXTERNAL_LENGTH_HOLDOUT_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args();self_test()
    if a.self_test:return
    roster.LIMIT=576
    plan,data,prov=roster.cohort(start_after=288,expected_prefix_sha256=PREVIOUS_SHA)
    print("FD_EXTERNAL_LENGTH_HOLDOUT_COHORT",json.dumps({k:v for k,v in prov.items() if k!="checklist_sources"},allow_nan=False),flush=True)
    _,source,nhm,wang=source_tables(include_values=True)
    ll=nhm["corolla_tube_range_cm"];sy=wang["unannotated_without_four_families"];panels={}
    for native,mode,length,minimum in SPECS:
        key="__".join([native,mode,length,"n"+str(minimum)])
        tab={s:(int(sy[s]=="bilateral"),float((np.log2(ll[s][0])+np.log2(ll[s][1]))/2 if length=="midpoint" else np.log2(ll[s][1]))) for s in sy.keys()&ll.keys()}
        sc=roster.SCOPES[0 if native=="loose" else 1];rows=[]
        for z in plan:
            i=z["target"];m=z["donor"]
            if i not in data or m not in data:continue
            t=data[i][sc]&tab.keys();s=data[m][sc]&tab.keys();d=s-data[i][roster.SCOPES[0]]
            if min(len(t),len(d))<minimum:continue
            rows.append(z|{"t":t,"s":d if mode=="disjoint" else s})
        perms=999 if key==PRIMARY else 199
        panels[key]={"all_pairs":test(rows,tab,False,perms),"matched_symmetry_cells":test(rows,tab,True,perms)}
    P=panels[PRIMARY]["all_pairs"].get("group_weightings",{})
    B=P.get("block",{}).get("tropical_minus_north",{});D=P.get("donor",{}).get("tropical_minus_north",{})
    def supported(q,conditional=False):
        m=q.get("metrics",{}).get("total",{});ci=m.get("residual_ci" if conditional else "observed_ci");lr=m.get("residual_leave_range")
        return bool(q.get("null_gate" if conditional else "geography_gate") and ci and ci[0]>0 and (not conditional or (m.get("coherent_p",1)<.05 and lr and lr[0]>0)))
    verdict={"raw_regional_contrast_reproduced":supported(B) and supported(D),"conditional_length_residual_supported":supported(B,True) and supported(D,True)}
    out={"primary":PRIMARY,"verdict":verdict,"rosters":prov,"source":source,"panels":panels,"warning":"After seeing the 288-target discovery set, freeze positive tropical-minus-north source-relative mean log2 tube-length difference on new target lists only. Donors, blocks, species and underlying trait publications may overlap. Same metadata plan prefix required; cap 576, no replacement of failed lists. Primary raw support uses both group bootstrap intervals and minimum geographic support; stronger conditional test fixes genus and symmetry, one endpoint, requires both p<.05, residual intervals and leave-group signs. Matched-cell decomposition exact composition + within-radial + within-bilateral, only when each cell has >=2 species; diagnostic, not causal. Remaining sensitivities descriptive. No confirmed direct pollinator or evolutionary response."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    for k,v in panels.items():print("FD_EXTERNAL_LENGTH_HOLDOUT_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
    print("FD_EXTERNAL_LENGTH_HOLDOUT_COMPLETE",json.dumps(verdict,allow_nan=False),flush=True)
if __name__=="__main__":main()
