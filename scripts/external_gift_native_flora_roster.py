#!/usr/bin/env python3
"""External-only GIFT 3.1 native roster feasibility for published floral traits.

Geography plan adapted from fd PR72, commit 2b2d8ea66ada2d808d519a9658416fd3fe8fddfe.
Plan selection uses metadata only, before trait files are loaded.
"""
import argparse,collections,hashlib,json,math,pathlib,time,urllib.parse
import numpy as np
from tropical_pollinator_source_preflight import get
from external_floral_replication_intake import name
GIFT="https://gift.uni-goettingen.de/api/extended/index3.1.php"
LIMIT=96
MAX_KM=2000
SCOPES=("native_no_positive_question_flags","strict_native")
def fetch(query,**params):
    url=GIFT+"?"+urllib.parse.urlencode({"query":query,**params})
    raw=get(url,maxbytes=80_000_000)
    x=json.loads(raw)
    if isinstance(x,dict):x=x.get("data",x.get("result"))
    if not isinstance(x,list) or not all(isinstance(v,dict) for v in x):raise ValueError("unexpected GIFT schema")
    return x,{"url":url,"sha256":hashlib.sha256(raw).hexdigest(),"rows":len(x)}
def num(v):
    try:
        f=float(v);return f if math.isfinite(f) else None
    except (ValueError,TypeError):return None
def ident(v):
    if v is None or not str(v).isdigit():raise ValueError("invalid publisher ID")
    return str(int(v))
def flag(v):
    return "missing" if v is None or str(v).strip()=="" else str(v).strip()
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

def native_names(i,record,zid):
    rows,pin=fetch("checklists",listid=ident(record["list_ID"]),taxonid=zid,namesmatched=0,filter="native")
    by=collections.defaultdict(list)
    cols=sorted({c for r in rows for c in r})
    if rows:
        assert {"entity_ID","list_ID","work_ID","genus_ID","species","native","questionable","quest_native"}<=set(cols),cols
    for r in rows:
        assert ident(r["entity_ID"])==i and ident(r["list_ID"])==ident(record["list_ID"])
        if r.get("work_ID") not in (None,""):by[ident(r["work_ID"])].append(r)
    out={s:{} for s in SCOPES};excluded=collections.Counter()
    for w,rs in by.items():
        n={flag(r.get("native")) for r in rs};q={flag(r.get("questionable")) for r in rs};s={flag(r.get("quest_native")) for r in rs}
        if n!={"1"} or not q<={"0","missing"} or not s<={"0","missing"}:
            excluded["nonliteral_native_or_question_flag"]+=1;continue
        names={name(r.get("species")) for r in rs};g={str(r.get("genus_ID")) for r in rs}
        if len(names)!=1 or None in names or len(g)!=1 or next(iter(g)) in {"None","","0"}:
            excluded["unresolved_binomial_or_genus"]+=1;continue
        v=(next(iter(names)),next(iter(g)))
        out[SCOPES[0]][w]=v
        if q=={"0"} and s=={"0"}:out[SCOPES[1]][w]=v
    assert set(out[SCOPES[1]])<=set(out[SCOPES[0]])
    return out,{"source":pin,"columns":cols,"excluded":dict(excluded),"counts":{s:len(v) for s,v in out.items()}}
def consistent_names(data):
    byid=collections.defaultdict(set);byname=collections.defaultdict(set)
    for v in data.values():
        for w,(n,g) in v[SCOPES[0]].items():
            byid[w].add((n,g));byname[n].add(w)
    badid={k for k,v in byid.items() if len(v)!=1}
    badname={k for k,v in byname.items() if len(v)!=1}
    clean={w:next(iter(v))[0] for w,v in byid.items() if w not in badid and next(iter(v))[0] not in badname}
    out={i:{s:{clean[w] for w in r if w in clean} for s,r in v.items()} for i,v in data.items()}
    return out,{"distinct_work_ids":len(byid),"conflicting_work_ids":len(badid),"ambiguous_binomials":len(badname),"clean_distinct_binomials":len(clean)}
def cohort():
    pins=[];rows=[]
    for query,params in [("regions",{}),("lists",{}),("taxonomy",{}),("geoentities_env_misc",{"envvar":"longitude,latitude"}),("overlap",{})]:
        v,p=fetch(query,**params)
        pins.append(p|{"canonical_rows_sha256":hashlib.sha256(("\n".join(sorted(json.dumps(x,sort_keys=True,ensure_ascii=False,separators=(",",":")) for x in v))+"\n").encode()).hexdigest()})
        rows.append(v)
    reg,lst,tax,env,ov=rows
    eligible,zid=eligible_lists(reg,lst,tax)
    plan,census=geo_plan(eligible,env,ov)
    assert plan and len(plan)<=LIMIT
    # Coordinates for the actual selected entities, not geography chosen
    # from trait availability. Region is absolute latitude <=23.5 tropical.
    coord={ident(r["entity_ID"]):(num(r.get("latitude")),num(r.get("longitude"))) for r in env}
    manifest=[{"target":i,"donor":m,"target_list":ident(eligible["Island"][i]["list_ID"]),"donor_list":ident(eligible["Mainland"][m]["list_ID"]),"region":"tropical" if abs(coord[i][0])<=23.5 else "north" if coord[i][0]>0 else "south","block":str(g)} for i,m,d,g in plan]
    # Manifest digest and metadata census recorded before reading any trait values.
    plan_sha=hashlib.sha256(json.dumps(manifest,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    print("FD_EXTERNAL_GIFT_GEOGRAPHY",json.dumps({"census":census,"metadata_pins":pins,"plan_sha256":plan_sha,"planned_regions":dict(collections.Counter(x["region"] for x in manifest))},allow_nan=False),flush=True)
    data={};sources=[];errors=[]
    for kind,ids in (("Island",sorted({x["target"] for x in manifest})),("Mainland",sorted({x["donor"] for x in manifest}))):
        for j,i in enumerate(ids):
            try:
                data[i],z=native_names(i,eligible[kind][i],zid)
                sources.append(z)
            except Exception as e:errors.append({"entity_ID":i,"list_ID":eligible[kind][i]["list_ID"],"error_type":type(e).__name__,"detail":str(e)[:180]})
            time.sleep(.25)
    native,audit=consistent_names(data)
    return manifest,native,{"metadata_pins":pins,"geography":census,"plan_sha256":plan_sha,"planned_regions":dict(collections.Counter(x["region"] for x in manifest)),"successful_lists":len(data),"checklist_sources":sources,"download_errors":errors,"taxonomy_audit":audit}
def self_test():
    d={"a":{SCOPES[0]:{"1":("plant one","2"),"2":("plant two","2")},SCOPES[1]:{"1":("plant one","2")}},
       "b":{SCOPES[0]:{"3":("plant two","2")},SCOPES[1]:{}}}
    t,a=consistent_names(d);assert t["a"][SCOPES[0]]=={"plant one"} and a["ambiguous_binomials"]==1
    assert abs(km((0,0),(0,1))-111.195)<.01
    print("FD_EXTERNAL_GIFT_ROSTER_SELF_TEST_PASS",flush=True)
