#!/usr/bin/env python3
"""One-study external measured morphology comparison, no raw export."""
import argparse, collections, itertools, json, math, pathlib
import numpy as np
from external_emln_spatial_audit import M,source_get,pivot,scalar
from external_euppollnet_curveball import simulate
DIR="26_Spatial_Chacoff_2020"
DRAWS=999
SEED=20261128

def positive_log(node,key):
    try:x=float(scalar(node,key))
    except (TypeError,ValueError):return None
    return math.log2(x) if math.isfinite(x) and x>0 else None

def source(n,l,e):
    N=pivot(n,"node_id");L=pivot(l,"layer");E=pivot(e,"interaction_id")
    assert (len(N),len(L),len(E))==(180,4,219)
    assert all(scalar(v,"type")=="space" and scalar(v,"year")=="2006" for v in L.values())
    def lookup(tab,key):
        out=collections.defaultdict(set)
        for k,d in tab.items():
            out[k].add(k)
            if scalar(d,key):out[scalar(d,key)].add(k)
        return out
    ni,li=lookup(N,"node_name"),lookup(L,"name")
    def resolve(key,mp):
        s=mp[key];assert len(s)==1,"ambiguous source key";return next(iter(s))
    roles={k:scalar(v,"type") for k,v in N.items()}
    traits={}
    for k,v in N.items():
        if roles[k]!="pollinator":continue
        x=[positive_log(v,t) for t in ("proboscis_length","body_length")]
        if None not in x:traits[k]=np.array(x)
    ids={v:i for i,v in enumerate(sorted(traits))}
    edges=collections.defaultdict(set);counts=collections.Counter()
    for d in E.values():
        a,b=resolve(scalar(d,"node_from"),ni),resolve(scalar(d,"node_to"),ni)
        s,t=resolve(scalar(d,"layer_from"),li),resolve(scalar(d,"layer_to"),li)
        assert s==t and {roles[a],roles[b]}=={"plant","pollinator"}
        assert math.isfinite(float(scalar(d,"weight"))) and float(scalar(d,"weight"))>0
        p,v=(a,b) if roles[a]=="plant" else (b,a)
        counts["positive_local_records"]+=1
        if v in ids:
            edges[s].add((p,ids[v]));counts["complete_trait_records"]+=1
    nets={};census=[]
    for s in sorted(L):
        ee=edges[s];p=collections.defaultdict(set)
        for a,b in ee:p[a].add(b)
        nvis=len({v for _,v in ee});census.append({"plant_taxa":len(p),"visitor_taxa":nvis,"unique_links":len(ee)})
        if len(p)>=5 and nvis>=5:
            names=sorted(p);nets[s]={"ix":{p:i for i,p in enumerate(names)},"rows":[p[k] for k in names]}
    return np.array([traits[k] for k in sorted(ids)]),nets,dict(counts),census

def distance(a,b,T):
    u=a-b;v=b-a
    if not u or not v:return None
    x=T[sorted(u)];y=T[sorted(v)]
    return np.abs(x[:,None,:]-y[None,:,:]).mean((0,1))

def conditional(a,b,sa,sb,T):
    ob=distance(a,b,T)
    if ob is None:return None
    shared=len(a&b);z=[];halves=[0,0]
    for i,(x,y) in enumerate(zip(sa,sb)):
        assert len(x)==len(a) and len(y)==len(b)
        if len(x&y)!=shared:continue
        d=distance(x,y,T);assert d is not None
        z.append(d);halves[int(i>=len(sa)//2)]+=1
    if len(z)<50 or min(halves)<20:return {"accepted":len(z),"eligible":False}
    nul=np.mean(z,axis=0)
    return {"accepted":len(z),"eligible":True,"obs":ob,"null":nul,"gain":nul-ob}

def aggregate(rows):
    if not rows:return None
    by=collections.defaultdict(list)
    for q in rows:by[q["plant"]].append(q)
    def avg(k):return np.mean([np.mean([x[k] for x in v],axis=0) for v in by.values()],axis=0)
    gain=avg("gain")
    return {"comparisons":len(rows),"focal_plant_taxa":len(by),"site_pairs":len({(q["a"],q["b"]) for q in rows}),
        "sites":len({k for q in rows for k in (q["a"],q["b"])}),
        "observed_mean_abs_log2_distance":avg("obs").tolist(),"null_mean_abs_log2_distance":avg("null").tolist(),
        "null_minus_observed_by_trait":gain.tolist(),"joint_gain":float(gain.mean())}

def analyze(T,nets):
    rng=np.random.default_rng(SEED);bank={};changes=[]
    for k,n in sorted(nets.items()):
        sim,diag=simulate(n["rows"],DRAWS,rng);bank[k]=list(zip(*sim));changes.append(diag["sample_changes"])
    rows=[];rejected=0;candidate=0
    for (a,A),(b,B) in itertools.combinations(sorted(nets.items()),2):
        for p in sorted(A["ix"].keys()&B["ix"].keys()):
            ia,ib=A["ix"][p],B["ix"][p];u,v=A["rows"][ia],B["rows"][ib]
            q=conditional(u,v,bank[a][ia],bank[b][ib],T)
            if q is None:continue
            candidate+=1
            if not q["eligible"]:rejected+=1;continue
            rows.append(dict(q,a=a,b=b,plant=p,min_exclusive=min(len(u-v),len(v-u))))
    main=aggregate(rows)
    return {"traits":["log2_mouthpart_length","log2_body_length"],"candidate_pairs":candidate,"rejected_pairs":rejected,
        "complete_trait_visitor_nodes":len(T),"eligible_sites":len(nets),"site_changed_trades":changes,
        "primary":main,"descriptive_support_gate":bool(main and main["focal_plant_taxa"]>=10 and main["sites"]>=3 and main["site_pairs"]>=3),
        "leave_one_site":[aggregate([q for q in rows if k not in (q["a"],q["b"])]) for k in sorted(nets)],
        "both_exclusive_ge2":aggregate([q for q in rows if q["min_exclusive"]>=2]),
        "scope":"One study, original morphotaxa, no cross-study inference. Positive gain means smaller morphological distance than conditional null. Mouthparts not homologous reach across orders; no pollination equivalence or geographical slope."}

def self_test():
    T=np.array([[0.,0.],[1.,2.],[2.,4.]])
    assert np.allclose(distance({0},{1,2},T),[1.5,3.])
    assert distance({0},{0,1},T) is None
    print("FD_EMLN_TRAIT_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    m=json.loads(M.read_text());tables=[];pins={}
    for f in ("nodes.csv","layers.csv","interactions.csv"):
        d,h=source_get(m,DIR,f);tables.append(d);pins[f]=h
    T,nets,c,census=source(*tables);r=analyze(T,nets);r.update(source_pins=pins,source_counts=c,layer_census=census)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2,allow_nan=False))
    print("FD_EMLN_TRAIT_RESULT",json.dumps(r,allow_nan=False),flush=True)
if __name__=="__main__":main()
