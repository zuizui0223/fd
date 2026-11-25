#!/usr/bin/env python3
"""Post-pilot, exact-cohort source-composition and guild-specific genus audit.

Reproduce the 96-pair official-source intake using signed prior artifact
metadata and byte hashes; fail closed on any source drift. No new geographic
selection or lowered primary support gate. GIFT self fertilization is not
autonomous selfing; publisher-coded syndromes are not necessarily observed
effective pollinators. All output is aggregated.
"""
import argparse,collections,hashlib,json,math,pathlib,time,zipfile
import numpy as np
from external_gift_reproductive_assurance_pilot import (
    fetch,eligible_lists,geo_plan,roster,carriers,contrast,summarize,
    SCOPES,ENDPOINTS,SEED,PIN,trait_rows)
PRIOR_SHA="1ee08e7e7f7f5326dc85223310d8fc85a91b248a9c56ae2216fa008c61148812"
GUILDS=("insect","wind")
METRICS=("raw_full_carriers","raw_common_reference_carriers","reference_equal_strata",
         "carrier_restriction_change","reference_reweighting_change","total_change")

def prior(path):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==PRIOR_SHA
    with zipfile.ZipFile(path) as z:
        names=[n for n in z.namelist() if n.endswith(".json")]
        assert len(names)==1
        a=json.loads(z.read(names[0]))
    assert a["geography_census"]["selected_island_entities"]==96
    assert a["geography_census"]["selected_mainland_comparators"]==59
    assert len(a["panels"])==6
    return a

def partition(rows,by):
    out=collections.defaultdict(list)
    for row in rows:
        i,g,sig,guild,y=row
        if by=="reference":k=(sig,guild)
        elif by=="guild_genus":k=(sig,genus_or_raise(g),guild)
        else:raise ValueError("unknown grouping")
        out[k].append(row)
    return out

def genus_or_raise(g):
    assert g not in ("","None")
    return g

def rate(rows):
    return float(np.mean([r[-1] for r in rows]))

def raw(t,s):
    z=[]
    for rows in (t,s):
        a=[[r for r in rows if r[3]==g] for g in GUILDS]
        if min(map(len,a))<5:return None
        z.append(rate(a[0])-rate(a[1]))
    return z[0]-z[1]

def ref_decomp(t,s):
    old=contrast(t,s,"common_reference_pattern")
    if old is None:return None
    a=partition(t,"reference");b=partition(s,"reference")
    sigs=sorted({r[2] for r in t}&{r[2] for r in s})
    keep=[sig for sig in sigs if all(len(x[(sig,g)])>=2 for x in (a,b) for g in GUILDS)]
    tt=[r for r in t if r[2] in keep];ss=[r for r in s if r[2] in keep]
    r0=raw(t,s);r1=raw(tt,ss)
    assert r0 is not None and r1 is not None
    v=[]
    for sig in keep:
        v.append((rate(a[(sig,"insect")])-rate(a[(sig,"wind")]))-
                 (rate(b[(sig,"insect")])-rate(b[(sig,"wind")])))
    r2=float(np.mean(v))
    assert np.isclose(r2,old["delta"],atol=1e-12)
    out=np.array([r0,r1,r2,r1-r0,r2-r1,r2-r0])
    assert np.isclose(out[3]+out[4],out[5],atol=1e-12)
    return {"values":out,"matched_target_taxa":len(tt),"matched_source_taxa":len(ss),
            "target_carriers":len(t),"source_carriers":len(s),"matched_reference_strata":len(keep)}

def guild_genus(t,s):
    """Match genus+reference within each guild across island/mainland.
    Distinct guilds need not occur within the same genus. This is not the
    prior within-genus cross-guild effect, nor an independent phylogenetic test.
    """
    a=partition(t,"guild_genus");b=partition(s,"guild_genus")
    vals=[];counts=[];unique=set();ng=set();n_info=0;n_strata=0;shared=0
    for guild in GUILDS:
        rows=[]
        for k in sorted(a.keys()&b.keys()):
            if k[-1]!=guild or min(len(a[k]),len(b[k]))<2:continue
            i={r[0] for r in a[k]};j={r[0] for r in b[k]}
            assert len(i)==len(a[k]) and len(j)==len(b[k])
            informative=i!=j and len({r[-1] for r in a[k]+b[k]})>=2
            # Count every matched stratum, but require at least one stratum
            # per guild with actual roster and outcome variation.
            rows.append((rate(a[k])-rate(b[k]),len(i),len(j),informative))
            unique|=i|j;ng.add(k[1]);shared+=len(i&j)
        if not rows:return None
        n=np.sum([[r[1],r[2]] for r in rows],axis=0)
        if min(n)<5 or not any(r[3] for r in rows):return None
        vals.append([float(np.mean([r[0] for r in rows])),
                     float(np.average([r[0] for r in rows],weights=[min(r[1],r[2]) for r in rows]))])
        counts.extend(map(int,n));n_strata+=len(rows);n_info+=sum(r[3] for r in rows)
    return {"delta":vals[0][0]-vals[1][0],"target_contrast":vals[0][0],
            "source_contrast":vals[1][0],"overlap_weighted_delta":vals[0][1]-vals[1][1],
            "cell_counts":counts,"matched_strata":n_strata,"informative_strata":n_info,
            "distinct_species":unique,"distinct_genera":ng,"shared_species_strata_total":shared}

def decomp_summary(rows):
    if not rows:return {"pairs":0,"status":"not_estimable"}
    out={"pairs":len(rows),"metric_order":list(METRICS),
         "median_carrier_counts":{k:float(np.median([r[k] for r in rows])) for k in (
             "matched_target_taxa","matched_source_taxa","target_carriers","source_carriers","matched_reference_strata")}}
    rng=np.random.default_rng(SEED)
    for key in ("block","donor"):
        groups=collections.defaultdict(list)
        for r in rows:groups[r[key]].append(r["values"])
        a=np.array([np.mean(groups[k],axis=0) for k in sorted(groups)])
        point=a.mean(0)
        assert np.isclose(point[3]+point[4],point[5],atol=1e-12)
        n=len(a);boot=np.array([a[rng.integers(n,size=n)].mean(0) for _ in range(1999)]) if n>=2 else None
        out[key]={"groups":n,"mean":point.tolist(),"interval":np.quantile(boot,[.025,.975],axis=0).T.tolist() if boot is not None else None}
    out["status"]="exploratory_estimable" if len(rows)>=60 and min(out[k]["groups"] for k in ("block","donor"))>=30 else "insufficient_support"
    return out

def genus_summary(rows):
    if not rows:return {"pairs":0,"status":"not_estimable"}
    # Reuse group-equal bootstrap, but give the two guild-specific changes
    # their correct labels rather than the original raw cross-guild labels.
    a=summarize(rows)
    a["metric_order"]=["island_minus_mainland_insect_minus_wind",
        "insect_island_minus_mainland","wind_island_minus_mainland","overlap_weighted_delta"]
    a["distinct_species"]=len(set().union(*(r["distinct_species"] for r in rows)))
    a["distinct_genera"]=len(set().union(*(r["distinct_genera"] for r in rows)))
    a["median_informative_strata"]=float(np.median([r["informative_strata"] for r in rows]))
    a["median_shared_species_strata_total"]=float(np.median([r["shared_species_strata_total"] for r in rows]))
    a["cell_order"]=["insect_island","insect_mainland","wind_island","wind_mainland"]
    return a

def self_test():
    s=(("1",),("2",));t=[];m=[]
    for g,y,z in (("insect",1,0),("wind",0,1)):
        for i in range(6):
            t.append((f"t{g}{i}","A" if g=="insect" else "B",s,g,y))
            m.append((f"m{g}{i}","A" if g=="insect" else "B",s,g,z))
    assert np.allclose(ref_decomp(t,m)["values"],[2,2,2,0,0,0])
    q= guild_genus(t,m)
    assert q["delta"]==2 and len(q["distinct_genera"])==2
    assert guild_genus(t,t) is None
    assert contrast(t,m,"common_reference_pattern_and_genus") is None
    print("EXTERNAL_ASSURANCE_EXACT_COHORT_SELF_TEST_PASS")

def main():
    p=argparse.ArgumentParser();p.add_argument("--prior",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    old=prior(a.prior)
    reg,rpin=fetch("regions");lst,lpin=fetch("lists");tax,tpin=fetch("taxonomy")
    env,epin=fetch("geoentities_env_misc",envvar="longitude,latitude");over,opin=fetch("overlap")
    assert [rpin,lpin,tpin,epin,opin]==old["metadata_pins"],"publisher geography/metadata drift"
    e,zid=eligible_lists(reg,lst,tax);plan,census=geo_plan(e,env,over)
    assert census==old["geography_census"],"cohort drift"
    ts=[]
    for k in ("guild_broad","self_fertilization"):
        tid,digest,n=PIN[k];v,d=trait_rows(tid,"strict")
        assert len(d["sources"])==1 and d["sources"][0]["sha256"]==digest and len(v)==n
        ts.append(v)
    # Prior artifact stores exact checklist source URL, hash and counts.
    prior_files={v["source"]["url"]:v for v in old["source_files"]}
    assert len(prior_files)==155
    data={};seen=set()
    for kind,ids in (("Island",sorted({i for i,m,d,g in plan})),("Mainland",sorted({m for i,m,d,g in plan}))):
        for i in ids:
            v,pin=roster(i,e[kind][i],zid)
            url=pin["source"]["url"]
            assert url in prior_files and pin==prior_files[url],"checklist source drift"
            data[i]=v;seen.add(url);time.sleep(.15)
    assert seen==set(prior_files)
    panels={}
    for scope in SCOPES:
        rows=[];original={k:[] for k in ENDPOINTS}
        genus_in=[];genus_dis=[];genus_common=[]
        for i,m,d,g in plan:
            t=carriers(data[i][scope],ts);s=carriers(data[m][scope],ts)
            meta={"block":str(g),"donor":m}
            for k in ENDPOINTS:
                v=contrast(t,s,k)
                if v is not None:original[k].append(v|meta)
            v=ref_decomp(t,s)
            if v is not None:rows.append(v|meta)
            vi=guild_genus(t,s)
            # Exclude every listed focal native name, not only matched carriers.
            sd=[r for r in s if r[0] not in data[i][scope]]
            vd=guild_genus(t,sd)
            if vi is not None:genus_in.append(vi|meta)
            if vd is not None:
                genus_dis.append(vd|meta)
                if vi is not None:genus_common.append(vi|meta)
        for k in ENDPOINTS:
            oldp=old["panels"][scope+"__"+k];new=summarize(original[k])
            assert new["pairs"]==oldp["pairs"]
            if new["pairs"]:
                assert np.isclose(new["mean_island_equal_delta"],oldp["mean_island_equal_delta"],atol=1e-12)
        p0=decomp_summary(rows)
        assert p0["pairs"]==old["panels"][scope+"__common_reference_pattern"]["pairs"]
        p0["raw_full_sample"]=summarize(original["raw"])
        panels[scope+"__reference_decomposition"]=p0
        for label,rr in (("inclusive",genus_in),("disjoint",genus_dis),("inclusive_on_disjoint_eligible",genus_common)):
            panels[scope+"__guild_genus_"+label]=genus_summary(rr)
    assert len(panels)==8
    out={"prior_artifact_sha256":PRIOR_SHA,"exact_source_files_reproduced":len(seen),
         "geography_census":census,"panels":panels,"warning":"Post-pilot conditional diagnostic, same source bytes. Exact reference-carrier/reweighting decomposition; same-guild genus matching is a different estimand from the frozen cross-guild within-genus test. Disjoint comparison excludes focal native names in corresponding scope; not an independent replication, direct pollination, autonomous selfing or causal evolution. No new geography and no relaxed support gate."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    for k,v in panels.items():print("EXTERNAL_ASSURANCE_DIAGNOSTIC",k,json.dumps(v,allow_nan=False),flush=True)
    print("EXTERNAL_ASSURANCE_DIAGNOSTIC_COMPLETE",json.dumps({"panels":len(panels),"exact_source_files_reproduced":len(seen)}),flush=True)
if __name__=="__main__":main()
