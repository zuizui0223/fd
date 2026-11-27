#!/usr/bin/env python3
"""Independent publisher-only native flora symmetry x continuous corolla length.

Post-coverage exploratory, metadata-defined 288 targets. The first 96
must reproduce the earlier geography-only manifest before scoring traits.
No new claim of direct observations, realized pollinators, or isolation.
"""
import argparse,collections,hashlib,json,pathlib
import numpy as np
from scipy import sparse
import external_gift_native_flora_roster as roster
from external_gift_floral_replication_coverage import source_tables
SEED=20261127
PREFIX_SHA="0fd9dd4096fb1ff6512d5c625312c397211c3d3d928549503f46d564413c6a05"
REGIONS=("north","tropical","south")
REG=REGIONS+("tropical_minus_north",)
METRICS=("joint_Q","additive_Q","joint_minus_fixed_marginal_Q","symmetry_Q","length_Q","bilateral_share","mean_log2_length_cm")
PRIMARY="loose__disjoint__unannotated_without_four_families__midpoint__n5"
SPECS=[
("loose","disjoint","unannotated_without_four_families","midpoint",5),
("loose","inclusive","unannotated_without_four_families","midpoint",5),
("loose","disjoint","all_publisher_binary","midpoint",5),
("loose","disjoint","unannotated","midpoint",5),
("strict","disjoint","unannotated_without_four_families","midpoint",5),
("loose","disjoint","unannotated_without_four_families","max",5),
("loose","disjoint","unannotated_without_four_families","midpoint",10),
]

def incidence(rosters,idx):
    rr=[];cc=[]
    for i,s in enumerate(rosters):
        for x in sorted(s):rr.append(i);cc.append(idx[x])
    return sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(rosters),len(idx)))

def pair_kernel(C,x):
    # Ordered distinct pairs, exponential similarity of log2 cm.
    # O(rows * unique lengths) recurrence, no pair-distance matrix.
    n=C.shape[0];prev=np.zeros(n);ans=np.sum(C*(C-1),axis=1)
    for j in range(1,len(x)):
        prev=(prev+C[:,j-1])*np.exp(-(x[j]-x[j-1]))
        ans+=2*C[:,j]*prev
    return ans

def evaluator(A,n,profiles):
    x=np.array(sorted({v[1] for v in profiles}))
    ix={v:i for i,v in enumerate(x)};K=len(x)
    col=np.array([int(s)*K+ix[l] for s,l in profiles]);N=A.shape[1];den=n*(n-1)
    def fn(code):
        Z=sparse.csr_matrix((np.ones(N),(np.arange(N),col[code])),shape=(N,2*K))
        C=(A@Z).toarray();a=C[:,:K];b=C[:,K:];m=a+b
        nb=b.sum(1);na=n-nb
        same=(na*(na-1)+nb*(nb-1))/den
        lk=pair_kernel(m,x)/den
        jk=(pair_kernel(a,x)+pair_kernel(b,x))/den
        qs=1-same;ql=1-lk;jq=1-jk
        return np.column_stack([jq,(qs+ql)/2,jq-(1-same*lk),qs,ql,nb/n,(m@x)/n])
    return fn

def groups(rows,tag):
    # Group equal weights within regions. Shared donor resampling
    # remains coherent when a donor occurs in two regions.
    keys=sorted({z[tag] for z in rows});ix={x:i for i,x in enumerate(keys)}
    pairs=sorted({(z["region"],z[tag]) for z in rows})
    px={p:i for i,p in enumerate(pairs)}
    g=np.array([px[(z["region"],z[tag])] for z in rows]);n=np.bincount(g)
    G=sparse.csr_matrix((1/n[g],(g,np.arange(len(rows)))),shape=(len(pairs),len(rows)))
    region=np.array([REGIONS.index(p[0]) for p in pairs]);key=np.array([ix[p[1]] for p in pairs])
    def agg(B,w):
        out=np.full((3,B.shape[1]),np.nan)
        for i in range(3):
            keep=region==i;ww=w[key[keep]]
            if ww.sum()>0:out[i]=np.average(B[keep],axis=0,weights=ww)
        return np.vstack([out,out[1]-out[0]])
    return G,agg,len(keys)

def evaluate(rows,tab,perms,boots=999):
    if not rows:return {"status":"empty","pairs":0}
    I=len(rows);rr=[z["t"] for z in rows];ss=[z["s"] for z in rows]
    names=sorted(set().union(*rr,*ss));ix={s:i for i,s in enumerate(names)}
    A=incidence(rr+ss,ix);n=np.asarray(A.sum(1)).ravel();assert min(n)>=5
    v=[tab[s] for s in names];profiles=sorted(set(v));pi={v:i for i,v in enumerate(profiles)}
    code=np.array([pi[z] for z in v]);fn=evaluator(A,n,profiles)
    genera=np.array([s.split()[0] for s in names]);_,gc=np.unique(genera,return_inverse=True)
    order=np.argsort(gc,kind="stable");movable=np.zeros(len(names),bool)
    by=collections.defaultdict(list)
    for i,g in enumerate(gc):by[g].append(i)
    for z in by.values():
        if len(set(code[z]))>1:movable[z]=True
    support={}
    for r in REGIONS:
        ids=np.flatnonzero(np.array([z["region"]==r for z in rows]))
        pres=np.asarray(A[np.r_[ids,ids+I]].sum(0)).ravel()>0
        support[r]={"pairs":len(ids),"blocks":len({rows[i]["block"] for i in ids}),"donors":len({rows[i]["donor"] for i in ids}),"movable_species":int((pres&movable).sum()),"movable_genera":int(len(set(genera[pres&movable]))),"median_target_carriers":float(np.median(n[ids])) if len(ids) else None,"median_donor_carriers":float(np.median(n[ids+I])) if len(ids) else None}
    def outcome(c):
        y=fn(c);return y[:I]-y[I:]
    Y=outcome(code);E=np.zeros_like(Y)
    designs={tag:groups(rows,tag) for tag in ("block","donor")}
    rng=np.random.default_rng(SEED)
    draws={tag:np.empty((perms,4,len(METRICS))) for tag in designs}
    for k in range(perms):
        u=np.empty_like(code);u[order]=code[np.lexsort((rng.random(len(code)),gc))]
        z=outcome(u);E+=z/perms
        for tag,(G,agg,N) in designs.items():draws[tag][k]=agg(G@z,np.ones(N))
    result={}
    for tag,(G,agg,N) in designs.items():
        B=G@Y;EN=G@E;obs=agg(B,np.ones(N));ex=agg(EN,np.ones(N));res=obs-ex
        raw=[];resid=[];leave=[]
        for k in range(boots):
            w=np.bincount(rng.integers(N,size=N),minlength=N)
            raw.append(agg(B,w));resid.append(agg(B-EN,w))
        for k in range(N):
            w=np.ones(N);w[k]=0;leave.append(agg(B-EN,w))
        raw=np.array(raw);resid=np.array(resid);leave=np.array(leave);out={}
        for j,r in enumerate(REG):
            req=[r] if r in REGIONS else ["north","tropical"]
            if any(support[q]["pairs"]==0 for q in req):continue
            gate=all(support[q]["pairs"]>=50 and support[q]["blocks"]>=30 and support[q]["donors"]>=20 and support[q]["movable_species"]>=30 and support[q]["movable_genera"]>=10 for q in req)
            out[r]={"gate":bool(gate),"metrics":{}}
            for i,k in enumerate(METRICS):
                x=draws[tag][:,j,i];o=obs[j,i];p=min(1.,2*min((1+int(np.sum(x<=o+1e-12)))/(perms+1),(1+int(np.sum(x>=o-1e-12)))/(perms+1)))
                def ci(z):
                    v=z[:,j,i];v=v[np.isfinite(v)]
                    return np.quantile(v,[.025,.975]).tolist() if len(v)>=.9*boots else None
                lv=leave[:,j,i]
                out[r]["metrics"][k]={"observed_delta":float(o),"observed_ci":ci(raw),"genus_null_mean":float(ex[j,i]),"residual":float(res[j,i]),"residual_ci":ci(resid),"coherent_p":float(p),"residual_leave_range":[float(lv.min()),float(lv.max())] if np.isfinite(lv).all() else None}
        result[tag]=out
    return {"status":"estimated","pairs":I,"species":len(names),"profiles":len(profiles),"permutations":perms,"support":support,"group_weightings":result}

def self_test():
    p=[(0,0.),(0,1.),(1,0.),(1,2.)]
    code=np.array([0,0,1,2,3,3]);A=sparse.csr_matrix(np.ones((1,len(code))));n=np.array([len(code)])
    v=np.array([p[i] for i in code]);ds=(v[:,0,None]!=v[None,:,0]).astype(float);dl=1-np.exp(-np.abs(v[:,1,None]-v[None,:,1]))
    mask=~np.eye(len(code),dtype=bool);a=ds[mask].mean();b=dl[mask].mean();j=(1-(1-ds)*(1-dl))[mask].mean()
    got=evaluator(A,n,p)(code)[0]
    assert np.allclose(got[:5],[j,(a+b)/2,j-(1-(1-a)*(1-b)),a,b],atol=1e-12)
    assert np.isclose(got[5],v[:,0].mean()) and np.isclose(got[6],v[:,1].mean())
    print("FD_EXTERNAL_CONTINUOUS_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true");p.add_argument("--out",type=pathlib.Path);a=p.parse_args()
    self_test()
    if a.self_test:return
    roster.LIMIT=288
    plan,data,provenance=roster.cohort()
    prefix=hashlib.sha256(json.dumps(plan[:96],sort_keys=True,separators=(",",":")).encode()).hexdigest()
    assert len(plan)==288 and prefix==PREFIX_SHA,"metadata plan changed"
    print("FD_EXTERNAL_CONTINUOUS_COHORT",json.dumps({"prefix_sha256":prefix,"pairs":len(plan),"provenance":{k:v for k,v in provenance.items() if k!="checklist_sources"}},allow_nan=False),flush=True)
    carriers,source,nhm,wang=source_tables(include_values=True)
    panels={}
    for native,mode,scope,length,minimum in SPECS:
        key="__".join([native,mode,scope,length,"n"+str(minimum)])
        sc=roster.SCOPES[0 if native=="loose" else 1];sy=wang[scope];ll=nhm["corolla_tube_range_cm"]
        tab={s:(int(sy[s]=="bilateral"),float((np.log2(ll[s][0])+np.log2(ll[s][1]))/2 if length=="midpoint" else np.log2(ll[s][1]))) for s in sy.keys()&ll.keys()}
        rows=[]
        for z in plan:
            i=z["target"];m=z["donor"]
            if i not in data or m not in data:continue
            t=data[i][sc]&tab.keys();s=data[m][sc]&tab.keys();disj=s-data[i][roster.SCOPES[0]]
            if min(len(t),len(disj))<minimum:continue
            rows.append(z|{"t":t,"s":s if mode=="inclusive" else disj})
        panels[key]=evaluate(rows,tab,999 if key==PRIMARY else 199)
        print("FD_EXTERNAL_CONTINUOUS_PROGRESS",key,json.dumps({"pairs":len(rows),"global_carriers":len(tab)}),flush=True)
    # Fixed 3-endpoint family for regional block-weighted residuals.
    primary=panels[PRIMARY].get("group_weightings",{}).get("block",{}).get("tropical_minus_north",{})
    donor=panels[PRIMARY].get("group_weightings",{}).get("donor",{}).get("tropical_minus_north",{})
    tests=[]
    for metric in METRICS[:3]:
        x=primary.get("metrics",{}).get(metric,{})
        tests.append((x.get("coherent_p",1.) if primary.get("gate") else 1.,metric,x))
    last=0
    for j,(pv,k,x) in enumerate(sorted(tests)):
        last=max(last,min(1.,pv*(len(tests)-j)))
        y=donor.get("metrics",{}).get(k,{})
        x["primary_holm_p"]=last
        def sign_ci(z):
            c=z.get("residual_ci");l=z.get("residual_leave_range")
            if not c or not l:return 0
            return 1 if c[0]>0 and l[0]>0 else -1 if c[1]<0 and l[1]<0 else 0
        x["primary_supported"]=bool(primary.get("gate") and donor.get("gate") and last<.05 and y.get("coherent_p",1)<.05 and sign_ci(x)!=0 and sign_ci(x)==sign_ci(y))
    out={"primary":PRIMARY,"primary_family_size":3,"metadata_prefix_sha256":prefix,"rosters":provenance,"source":source,"panels":panels,"warning":"Post-coverage exploratory external publisher-only analysis. Pairwise kernel distance 1-exp(-abs(log2 length ratio)); symmetry mismatch binary, joint product-complement, exact independent-marginal departure, finite-sample distinct pairs. Midpoint is geometric midpoint of published min/max, not measured mean. Unannotated and four-family exclusion do not establish direct evidence. Whole-pair coherent genus permutations conditional on observed carrier names and rosters. Block and donor bootstrap separate, not independent replication; require both residual interval and leave-group signs. Loose native flag missingness explicitly allowed; strict sensitivity. Not geographic isolation or realized pollination; no raw plant rows exported."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    for k,v in panels.items():print("FD_EXTERNAL_CONTINUOUS_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
    print("FD_EXTERNAL_CONTINUOUS_COMPLETE",json.dumps({"panels":len(panels),"primary_family_size":3}),flush=True)
if __name__=="__main__":main()
