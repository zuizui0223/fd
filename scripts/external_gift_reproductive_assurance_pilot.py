#!/usr/bin/env python3
"""Publisher-only island–mainland reproductive-assurance pilot.

Post-coverage exploratory design; selects geography before reading traits.
Nearest eligible mainland centroid is a comparison proxy, not a reconstructed
colonist donor. GIFT aggregated self fertilization is not autonomous selfing;
reported insect/wind modes need not be observed effective pollination.
No project-derived inputs. All paired analyses preserve identical references.
"""
import argparse,collections,hashlib,json,math,pathlib,time
import numpy as np
from external_gift_floral_trait_overlap import fetch,ident,num,trait_rows
from external_gift_coverage_source_audit import PIN,flag

SEED=20261125
LIMIT=96
MAX_KM=2000
MIN_CELL=5
SCOPES=("native_no_positive_question_flags","strict_native")
ENDPOINTS=("raw","common_reference_pattern","common_reference_pattern_and_genus")

def priority(i):
    return hashlib.sha256(("fd_gift_native_assurance_20261125:"+i).encode()).hexdigest()

def eligible_lists(regions,lists,taxonomy):
    taxa={ident(x["taxon_ID"]):x for x in taxonomy}
    target=[x for x in taxonomy if str(x["taxon_name"]).casefold()=="angiospermae"]
    assert len(target)==1
    z=target[0];zid=ident(z["taxon_ID"])
    complete={i for i,x in taxa.items() if num(x.get("lft")) is not None and num(x.get("rgt")) is not None and float(x["lft"])<=float(z["lft"]) and float(x["rgt"])>=float(z["rgt"])}
    geo={ident(x["entity_ID"]):x for x in regions}
    out={"Island":{},"Mainland":{}}
    for x in lists:
        i=ident(x["entity_ID"]);g=geo[i];kind=g.get("entity_class")
        if kind not in out:continue
        if str(x.get("restricted"))!="0" or str(x.get("suit_geo"))!="1":continue
        if str(x.get("native_indicated"))!="1" or ident(x["taxon_ID"]) not in complete:continue
        if str(x.get("subset","")).strip().casefold() not in {"all","native","native and naturalized","native and historically introduced"}:continue
        out[kind].setdefault(i,[]).append(x)
    return {k:{i:sorted(v,key=lambda x:(str(x["subset"]).casefold()!="native",int(x["list_ID"])))[0] for i,v in a.items()} for k,a in out.items()},zid

def km(a,b):
    lat1,lon1=map(math.radians,a);lat2,lon2=map(math.radians,b)
    h=math.sin((lat2-lat1)/2)**2+math.cos(lat1)*math.cos(lat2)*math.sin((lon2-lon1)/2)**2
    return 6371*2*math.asin(math.sqrt(min(1,max(0,h))))

def overlap_graph(rows,candidates):
    """Publisher overlap fractions. Filter candidate IDs before validating values.

    The official GIFT_no_overlap R function filters both endpoints to candidate
    IDs before applying the 0.1 fraction threshold. Never rescale out-of-range
    values or treat a relevant missing value as zero.
    """
    edges=collections.defaultdict(set);quarantine=set();counts=collections.Counter()
    reasons=collections.Counter();minmax=[None,None]
    for r in rows:
        a=ident(r["entity1"]);b=ident(r["entity2"])
        counts["raw_rows"]+=1
        if a not in candidates or b not in candidates:
            counts["irrelevant_rows"]+=1
            continue
        counts["candidate_rows"]+=1
        x=num(r.get("overlap12"));y=num(r.get("overlap21"))
        if x is None or y is None:
            reasons["missing_or_nonfinite"]+=1
            quarantine.update((a,b));continue
        if not (0<=x<=1 and 0<=y<=1):
            reasons["out_of_unit_interval"]+=1
            quarantine.update((a,b));continue
        counts["valid_candidate_rows"]+=1
        lo,hi=min(x,y),max(x,y)
        minmax=[lo if minmax[0] is None else min(lo,minmax[0]),
                hi if minmax[1] is None else max(hi,minmax[1])]
        if a!=b and hi>.1:
            edges[a].add(b);edges[b].add(a)
    audit={k:counts[k] for k in ("raw_rows","irrelevant_rows","candidate_rows","valid_candidate_rows")}
    audit.update({"invalid_candidate_rows":sum(reasons.values()),
                  "invalid_reasons":dict(reasons),"quarantined_candidate_entities":len(quarantine),
                  "eligible_candidate_entities":len(candidates),"valid_fraction_range":minmax,
                  "threshold":.1,"rule":"candidate pairs first; quarantine both endpoints of any invalid relevant pair"})
    assert audit["raw_rows"]==audit["irrelevant_rows"]+audit["candidate_rows"]
    assert audit["candidate_rows"]==audit["valid_candidate_rows"]+audit["invalid_candidate_rows"]
    return edges,quarantine,audit

def geo_plan(eligible,environment,overlap):
    coords={};bad=set()
    for r in environment:
        i=ident(r["entity_ID"]);lat=num(r.get("latitude"));lon=num(r.get("longitude"))
        if lat is None or lon is None or not -90<=lat<=90 or not -180<=lon<=180:continue
        if i in coords and coords[i]!=(lat,lon):bad.add(i)
        coords[i]=(lat,lon)
    for i in bad:coords.pop(i,None)
    candidates=(set(eligible["Island"])|set(eligible["Mainland"]))&coords.keys()
    edges,quarantine,overlap_audit=overlap_graph(overlap,candidates)
    print("EXTERNAL_ASSURANCE_OVERLAP_AUDIT",json.dumps(overlap_audit,allow_nan=False),flush=True)
    if not overlap_audit["valid_candidate_rows"]:
        raise ValueError("no validated overlap rows among candidate entities")
    # Do not infer a zero overlap from a malformed relevant record.
    coords={i:v for i,v in coords.items() if i not in quarantine}
    main=sorted(set(eligible["Mainland"])&coords.keys(),key=priority)
    mset=set();mkeep=[]
    for m in main:
        if edges[m]&mset:continue
        mset.add(m);mkeep.append(m)
    queues=collections.defaultdict(list)
    no_donor=0
    for i in sorted(set(eligible["Island"])&coords.keys(),key=priority):
        lat,lon=coords[i];g=(math.floor(lat/10),math.floor(lon/10))
        cand=sorted(((km(coords[i],coords[m]),m) for m in mkeep if m not in edges[i]),key=lambda p:(p[0],priority(p[1])))
        cand=[p for p in cand if p[0]<=MAX_KM]
        if not cand:no_donor+=1;continue
        queues[g].append((i,cand[0][1],cand[0][0],g))
    # Geographic round-robin; globally unique nonoverlapping target entities.
    chosen=[];keep=set();skip=0
    order=sorted(queues,key=lambda g:priority(str(g)))
    while order and len(chosen)<LIMIT:
        remaining=[]
        for g in order:
            q=queues[g]
            while q:
                i,m,d,block=q.pop(0)
                if edges[i]&keep:skip+=1;continue
                # Cross-cohort overlap among any selected targets and donors
                # would invalidate a simple separate geography comparison.
                if any(i in edges[z[1]] or m in edges[z[0]] for z in chosen):
                    skip+=1;continue
                chosen.append((i,m,d,block));keep.add(i);break
            if q:remaining.append(g)
            if len(chosen)>=LIMIT:break
        order=remaining
    assert not any(edges[i]&keep for i in keep)
    donors={m for i,m,d,g in chosen}
    assert not any(edges[i]&donors for i in keep)
    assert not any(edges[m]&donors for m in donors)
    out={"overlap_schema_audit":overlap_audit,"eligible_island_entities":len(eligible["Island"]),"eligible_mainland_entities":len(eligible["Mainland"]),
         "valid_island_centroids":len(set(eligible["Island"])&coords.keys()),"valid_mainland_centroids":len(main),
         "nonoverlapping_mainland_candidates":len(mkeep),"islands_without_comparator_within_2000km":no_donor,
         "selected_island_entities":len(chosen),"selected_mainland_comparators":len(donors),
         "selected_geographic_blocks":len({g for i,m,d,g in chosen}),
         "skipped_overlap_in_selection":skip,"centroid_distance_median_km":float(np.median([d for i,m,d,g in chosen])) if chosen else None}
    return chosen,out

def roster(i,record,zid):
    rows,pin=fetch("checklists",listid=ident(record["list_ID"]),taxonid=zid,namesmatched=0,filter="native")
    by=collections.defaultdict(list)
    for r in rows:
        assert ident(r["entity_ID"])==i and ident(r["list_ID"])==ident(record["list_ID"])
        if r.get("work_ID") in (None,""):continue
        by[ident(r["work_ID"])].append(r)
    out={s:{} for s in SCOPES};excluded=collections.Counter()
    for w,rs in by.items():
        # A contradictory duplicate is not resolved by a single favourable
        # row. Reject the entire work_ID when any retrieved row has native
        # or questionable status inconsistent with the cohort.
        n={flag(r.get("native")) for r in rs};q={flag(r.get("questionable")) for r in rs};s={flag(r.get("quest_native")) for r in rs}
        if n!={"1"} or not q<={"0","missing"} or not s<={"0","missing"}:
            excluded["contradiction_or_nonliteral_native"]+=1;continue
        gen={str(r.get("genus_ID")) for r in rs if r.get("genus_ID") not in (None,"")}
        if len(gen)!=1:excluded["genus_ambiguous_or_missing"]+=1;continue
        v=next(iter(gen));out[SCOPES[0]][w]=v
        if q=={"0"} and s=={"0"}:out[SCOPES[1]][w]=v
    assert out[SCOPES[1]].keys()<=out[SCOPES[0]].keys()
    return out,{"source":pin,"native_counts":{k:len(v) for k,v in out.items()},"excluded":dict(excluded),"duplicate_work_id_rows":sum(len(v)-1 for v in by.values())}

def carriers(r,traits):
    a,b=traits
    out=[]
    for i in sorted(r.keys()&a.keys()&b.keys()):
        g=a[i]["v"];v=b[i]["v"]
        if g not in ("insect","wind") or v not in ("present","absent"):continue
        if not a[i]["refs"] or not b[i]["refs"]:continue
        sig=(tuple(sorted(a[i]["refs"])),tuple(sorted(b[i]["refs"])))
        out.append((i,r[i],sig,g,int(v=="present")))
    return out

def contrast(target,source,kind):
    def partition(rows):
        out=collections.defaultdict(lambda:{"insect":[],"wind":[]})
        for i,gen,sig,g,y in rows:
            k=0 if kind=="raw" else sig if kind=="common_reference_pattern" else (sig,gen)
            out[k][g].append(y)
        return out
    t=partition(target);s=partition(source)
    common=[]
    min_stratum=MIN_CELL if kind=="raw" else 2
    for k in t.keys()&s.keys():
        n=[len(t[k][g]) for g in ("insect","wind")]+[len(s[k][g]) for g in ("insect","wind")]
        if min(n)<min_stratum:continue
        v=[np.mean(t[k]["insect"])-np.mean(t[k]["wind"]),np.mean(s[k]["insect"])-np.mean(s[k]["wind"])]
        common.append((n,v))
    if not common:return None
    ns=np.sum([n for n,v in common],axis=0)
    if min(ns)<MIN_CELL:return None
    # Equal strata conditional estimand. Report overlap-weighted separately.
    v=np.array([v for n,v in common])
    p=v.mean(0);q=np.average(v,axis=0,weights=[min(n) for n,_ in common])
    return {"delta":float(p[0]-p[1]),"target_contrast":float(p[0]),"source_contrast":float(p[1]),
            "overlap_weighted_delta":float(q[0]-q[1]),"cell_counts":ns.tolist(),"matched_strata":len(common)}

def summarize(rr):
    if not rr:return {"status":"not_estimable","pairs":0}
    out={"pairs":len(rr),"median_cell_counts":np.median([x["cell_counts"] for x in rr],axis=0).tolist(),
         "median_matched_strata":float(np.median([x["matched_strata"] for x in rr])),
         "mean_island_equal_delta":float(np.mean([x["delta"] for x in rr]))}
    rng=np.random.default_rng(SEED)
    for tag in ("block","donor"):
        gr=collections.defaultdict(list)
        for x in rr:gr[x[tag]].append([x["delta"],x["target_contrast"],x["source_contrast"],x["overlap_weighted_delta"]])
        a=np.array([np.mean(v,axis=0) for v in gr.values()])
        n=len(a)
        boot=np.array([a[rng.integers(n,size=n)].mean(0) for j in range(1999)]) if n>=2 else None
        out[tag]={"groups":n,"mean":a.mean(0).tolist(),"interval":np.quantile(boot,[.025,.975],axis=0).T.tolist() if boot is not None else None}
    gate=len(rr)>=60 and all(out[x]["groups"]>=30 for x in ("block","donor"))
    return out|{"status":"exploratory_estimable" if gate else "insufficient_support",
         "positive_direction_both_intervals":bool(gate and all(out[x]["interval"][0][0]>0 for x in ("block","donor"))),
         "metric_order":["island_minus_mainland","island_insect_minus_wind","mainland_insect_minus_wind","overlap_weighted_delta"]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:
        z=[(str(i),"G",(("1",),("2",)),g,y) for g,y in (("insect",1),("wind",0)) for i in range(10)]
        q=[(str(i),"G",(("1",),("2",)),g,0) for g in ("insect","wind") for i in range(10)]
        for k in ENDPOINTS:assert contrast(z,q,k)["delta"]==1
        assert contrast(z[:4],q,"raw") is None
        assert abs(km((0,0),(0,1))-111.195)<.01
        edges,bad,audit=overlap_graph([
            {"entity1":"9","entity2":"1","overlap12":None,"overlap21":"-1"},
            {"entity1":"1","entity2":"2","overlap12":".2","overlap21":".01"},
            {"entity1":"2","entity2":"3","overlap12":"1.01","overlap21":".4"},
            {"entity1":"1","entity2":"4","overlap12":".1","overlap21":".01"},
            {"entity1":"4","entity2":"5","overlap12":None,"overlap21":"0"},
        ],{"1","2","3","4","5"})
        assert bad=={"2","3","4","5"} and edges["1"]=={"2"}
        assert audit["irrelevant_rows"]==1 and audit["valid_candidate_rows"]==2
        assert audit["invalid_reasons"]=={"out_of_unit_interval":1,"missing_or_nonfinite":1}
        assert audit["valid_fraction_range"]==[.01,.2]
        print("EXTERNAL_GIFT_ASSURANCE_TEST_PASS");return
    reg,rpin=fetch("regions");lst,lpin=fetch("lists");tax,tpin=fetch("taxonomy")
    assert rpin["sha256"]=="d12d9129eb8aeb7aafb295ddf99b3ffdd810dcd628ab951ca31e560ab772b64b"
    assert lpin["sha256"]=="7cbf01e767f5abc061fb3622db031c95e4b24646b98d875eed9f1d9e742ba244"
    eligible,zid=eligible_lists(reg,lst,tax)
    assert len(eligible["Island"])==1253 and len(eligible["Mainland"])==751
    env,epin=fetch("geoentities_env_misc",envvar="longitude,latitude");over,opin=fetch("overlap")
    plan,census=geo_plan(eligible,env,over)
    print("EXTERNAL_ASSURANCE_GEOGRAPHY",json.dumps(census),flush=True)
    assert plan
    tables=[]
    for name in ("guild_broad","self_fertilization"):
        tid,pin,n=PIN[name]
        v,d=trait_rows(tid,"strict")
        assert len(d["sources"])==1 and d["sources"][0]["sha256"]==pin and len(v)==n
        tables.append(v)
    sources={};data={};errors=[]
    for kind,ids in (("Island",sorted({i for i,m,d,g in plan})),("Mainland",sorted({m for i,m,d,g in plan}))):
        for i in ids:
            try:data[i],sources[i]=roster(i,eligible[kind][i],zid)
            except Exception as e:errors.append({"entity_ID":i,"list_ID":eligible[kind][i]["list_ID"],"error_type":type(e).__name__,"detail":str(e)[:180]})
            time.sleep(.25)
    panels={}
    for scope in SCOPES:
        rr={k:[] for k in ENDPOINTS}
        counts=collections.Counter()
        for i,m,d,g in plan:
            if i not in data or m not in data:continue
            counts["downloaded_pairs"]+=1
            t=carriers(data[i][scope],tables);s=carriers(data[m][scope],tables)
            for k in ENDPOINTS:
                z=contrast(t,s,k)
                if z is not None:rr[k].append(z|{"block":str(g),"donor":m})
        for k in ENDPOINTS:
            p0=summarize(rr[k]);p0["downloaded_pairs"]=counts["downloaded_pairs"]
            key=scope+"__"+k
            panels[key]=p0
            print("EXTERNAL_ASSURANCE_PANEL",key,json.dumps(p0,allow_nan=False),flush=True)
    out={"metadata_pins":[rpin,lpin,tpin,epin,opin],"trait_pins":{k:PIN[k] for k in ("guild_broad","self_fertilization")},
         "geography_census":census,"source_files":list(sources.values()),"download_errors":errors,
         "panels":panels,"primary_post_coverage_exploratory":SCOPES[0]+"__common_reference_pattern",
         "warning":"Nonoverlap uses validated publisher overlap fractions >0.1 among candidate entities; any invalid relevant pair quarantines both endpoints. No unknown overlap is converted to zero. Not proof of geographic independence. Closest mainland centroid within 2000 km is not proven source or ecological match. Native missing-question flags are not verified absence of doubts; strict sensitivity requires all flags zero. Repeated taxa and source publications remain dependent. Reported self fertilization is not autonomous selfing and GIFT modes are not verified observed pollination. Geographic and donor bootstrap intervals conditional on observed species/reference inputs; no causal or evolutionary conclusion. No project-created botanical data."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
if __name__=="__main__":main()
