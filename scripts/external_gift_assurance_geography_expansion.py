#!/usr/bin/env python3
"""Metadata-fixed expansion of the official GIFT reproductive assurance pilot.

The original 96-pair geography and 155 checklist hashes must reproduce.
Expand the unchanged metadata-only ordering to at most 288 targets before
loading traits. New checklists are retrieved once, serially. New failures
are retained without replacement. Added islands are not independent studies.
"""
import argparse,collections,hashlib,json,pathlib,time,zipfile
import numpy as np
from scipy.sparse import csr_matrix
import external_gift_reproductive_assurance_pilot as pilot
from external_gift_assurance_exact_cohort import prior,PRIOR_SHA
from external_gift_coverage_source_audit import PIN
from external_gift_floral_trait_overlap import fetch,trait_rows

LIMIT=288
DRAWS=999
SCOPES=pilot.SCOPES
COHORTS=("original_96","added_targets","pooled")
MODES=("inclusive","disjoint","inclusive_on_disjoint_eligible")
SEED=pilot.SEED
NEW_PRIOR_SHA="85c15d81f65f29a1a8ad2ac3ad3ea2b57949d6684e2b9b08526556f56ca73800"
GUILDS=("insect","wind")

def matched_cells(t,s):
    def collect(a):
        x=collections.defaultdict(list)
        for row in a:x[(row[2],row[3])].append(row)
        return x
    a=collect(t);b=collect(s)
    keep=[]
    for sig in sorted({r[2] for r in t}&{r[2] for r in s}):
        cells=[a[(sig,g)] for g in GUILDS]+[b[(sig,g)] for g in GUILDS]
        if min(map(len,cells))>=2:keep.append(cells)
    if not keep:return None
    n=np.sum([[len(c) for c in cells] for cells in keep],axis=0)
    if min(n)<5:return None
    def m(c):return float(np.mean([r[-1] for r in c]))
    v=np.array([[m(c) for c in cells] for cells in keep])
    contrasts=np.column_stack([v[:,0]-v[:,1],v[:,2]-v[:,3]])
    p=contrasts.mean(0)
    d=float(p[0]-p[1])
    old=pilot.contrast(t,s,"common_reference_pattern")
    assert old is not None and np.isclose(d,old["delta"],atol=1e-12)
    return {"cells":keep,"delta":d,"target_contrast":float(p[0]),"source_contrast":float(p[1]),
            "overlap_weighted_delta":float(np.average(contrasts[:,0]-contrasts[:,1],weights=[min(map(len,c)) for c in keep])),
            "cell_counts":n.tolist(),"matched_strata":len(keep)}

def all_ids(rows):
    return {r[0] for row in rows for cells in row["cells"] for cell in cells for r in cell}

def coherent_genus_sensitivity(rows):
    if not rows:return {"status":"not_estimable"}
    mapping=collections.defaultdict(set)
    for row in rows:
        for cells in row["cells"]:
            for cell in cells:
                for r in cell:mapping[r[0]].add(r[1])
    conflict=sum(len(s)!=1 for s in mapping.values())
    if conflict:return {"status":"genus_identity_conflict","conflicting_species":conflict}
    genera=sorted(set().union(*mapping.values()))
    lookup={g:i for i,g in enumerate(genera)}
    rr=[];cc=[];den=[];num=[];owners=[];n_cells=0
    for i,row in enumerate(rows):
        for cells in row["cells"]:
            owners.append(i)
            for cell in cells:
                counts=collections.Counter();pos=collections.Counter()
                for r in cell:
                    j=lookup[r[1]];counts[j]+=1;pos[j]+=r[-1]
                for j,n in counts.items():
                    rr.append(n_cells);cc.append(j);den.append(n);num.append(pos[j])
                n_cells+=1
    A=csr_matrix((den,(rr,cc)),shape=(n_cells,len(genera)))
    B=csr_matrix((num,(rr,cc)),shape=A.shape)
    n=np.bincount(owners)
    W=csr_matrix((1/n[owners],(owners,np.arange(len(owners)))),shape=(len(rows),len(owners)))
    def delta(w):
        p=(B@w)/(A@w)
        q=p.reshape(-1,4)
        return W@(q[:,0]-q[:,1]-q[:,2]+q[:,3])
    base=delta(np.ones(len(genera)))
    assert np.allclose(base,[r["delta"] for r in rows],atol=1e-12)
    groups={}
    for name in ("block","donor"):
        labels=sorted({r[name] for r in rows})
        g=np.array([labels.index(r[name]) for r in rows])
        n0=np.bincount(g)
        groups[name]=(csr_matrix((1/n0[g],(g,np.arange(len(g)))),shape=(len(labels),len(g))),len(labels))
    rng=np.random.default_rng(SEED)
    draw={k:{"taxon":[],"joint":[]} for k in groups}
    for j in range(DRAWS):
        w=np.maximum(rng.exponential(size=len(genera)),1e-12)
        d=delta(w)
        for k,(H,G) in groups.items():
            q=H@d
            draw[k]["taxon"].append(float(q.mean()))
            draw[k]["joint"].append(float(q[rng.integers(G,size=G)].mean()))
    return {"status":"completed","distinct_species":len(mapping),"distinct_genera":len(genera),"draws":DRAWS,
         "groups":{k:{"groups":G,"mean":float((H@base).mean()),
             "genus_weight_interval":np.quantile(draw[k]["taxon"],[.025,.975]).tolist(),
             "joint_genus_and_group_weight_interval":np.quantile(draw[k]["joint"],[.025,.975]).tolist()}
             for k,(H,G) in groups.items()},
         "warning":"Shared positive genus weights in all carrier cells; exploratory dependence sensitivity, not calibrated CI."}

def summary(rows,do_genus=True):
    out=pilot.summarize(rows)
    if rows:
        out["distinct_matched_taxa"]=len(all_ids(rows))
        out["distinct_matched_reference_patterns"]=len({c[0][0][2] for row in rows for c in row["cells"]})
        if do_genus:
            out["genus_dependency_sensitivity"]=coherent_genus_sensitivity(rows)
    return out

def eval_cohort(plan,data,tables,scope,cohort,exclude_ref=None):
    result={m:[] for m in MODES};downloaded=0
    for i,m,d,g in plan:
        if i not in data or m not in data:continue
        downloaded+=1
        t=pilot.carriers(data[i][scope],tables);s=pilot.carriers(data[m][scope],tables)
        if exclude_ref is not None:
            t=[r for r in t if exclude_ref not in r[2][0] and exclude_ref not in r[2][1]]
            s=[r for r in s if exclude_ref not in r[2][0] and exclude_ref not in r[2][1]]
        q=[r for r in s if r[0] not in data[i][scope]]
        u=matched_cells(t,s);v=matched_cells(t,q)
        tag={"block":str(g),"donor":m}
        if u is not None:result["inclusive"].append(u|tag)
        if v is not None:
            assert u is not None
            result["disjoint"].append(v|tag);result["inclusive_on_disjoint_eligible"].append(u|tag)
    assert len(result["disjoint"])==len(result["inclusive_on_disjoint_eligible"])
    return {k:summary(v,do_genus=exclude_ref is None)|{"planned_pairs":len(plan),"downloaded_pairs":downloaded} for k,v in result.items()}

def canonical_rows(rows):
    # Sort complete row representations, preserving multiplicity and every
    # field; distinguishes ordering-only changes across current fetches.
    # Historical canonical hashes were not recorded, so equality here
    # cannot establish historical semantic equality.
    lines=sorted(json.dumps(r,sort_keys=True,ensure_ascii=False,separators=(",",":"),allow_nan=False) for r in rows)
    return hashlib.sha256(("\n".join(lines)+"\n").encode("utf-8")).hexdigest()

def load_comparison(path):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==NEW_PRIOR_SHA
    with zipfile.ZipFile(path) as z:
        names=[n for n in z.namelist() if n.endswith(".json")]
        assert len(names)==1
        return json.loads(z.read(names[0]))

def compare_pilots(old,new):
    assert new["geography_census"]==old["geography_census"],"pilot geography census differs"
    def manifests(a):return {x["source"]["url"]:x for x in a["source_files"]}
    m0=manifests(old);m1=manifests(new)
    assert len(m0)==len(m1)==155 and m0==m1,"prior checklist manifests differ"
    assert new["trait_pins"]==old["trait_pins"],"prior trait pins differ"
    assert set(new["panels"])==set(old["panels"])
    for k,p in old["panels"].items():
        q=new["panels"][k]
        assert p["pairs"]==q["pairs"] and p["status"]==q["status"] and p["downloaded_pairs"]==q["downloaded_pairs"]
        if p["pairs"]:
            assert np.isclose(p["mean_island_equal_delta"],q["mean_island_equal_delta"],atol=1e-12)
            for label in ("block","donor"):
                assert p[label]["groups"]==q[label]["groups"]
                assert np.allclose(p[label]["mean"],q[label]["mean"],atol=1e-12)
    return {"checklist_manifests_identical":len(m0),"six_panel_points_and_support_reproduced":True,
            "metadata_comparison":[{"url":a["url"],"old_sha256":a["sha256"],"new_sha256":b["sha256"],"old_rows":a["rows"],"new_rows":b["rows"],"raw_equal":a==b} for a,b in zip(old["metadata_pins"],new["metadata_pins"])]}

def metadata_snapshot(old,comparison):
    requests=[("regions",{}),("lists",{}),("taxonomy",{}),("geoentities_env_misc",{"envvar":"longitude,latitude"}),("overlap",{})]
    data=[];pins=[];audit=[]
    for j,(query,params) in enumerate(requests):
        rows,pin=fetch(query,**params)
        prev=old["metadata_pins"][j]
        fresh=canonical_rows(rows)
        exact=pin==prev
        info={"endpoint":query,"prior":prev,"current":pin,"current_canonical_rows_sha256":fresh,"raw_exact_old":exact}
        if comparison is not None:
            cmp=comparison["metadata_pins"][j]
            info["comparison_pilot"]=cmp
            info["raw_exact_comparison"]=pin==cmp
        if not exact:
            if comparison is None:raise AssertionError("publisher metadata drift: "+query)
            # Re-fetch only drifted endpoints and prove within-run canonical
            # stability; do not declare that the historic bytes are equivalent.
            rows2,pin2=fetch(query,**params)
            info["repeat_current"]=pin2
            info["repeat_canonical_rows_sha256"]=canonical_rows(rows2)
            info["current_repeat_canonical_equal"]=fresh==info["repeat_canonical_rows_sha256"]
            assert info["current_repeat_canonical_equal"],"unstable current metadata: "+query
        data.append(rows);pins.append(pin);audit.append(info)
        print("EXTERNAL_GIFT_METADATA_PIN_AUDIT",query,json.dumps(info,allow_nan=False),flush=True)
    return data,pins,audit

def self_test():
    sig=(("1",),("2",));t=[];s=[]
    for guild,ya,yb in (("insect",1,0),("wind",0,1)):
        for i in range(6):
            t.append((f"t{guild}{i}",f"A{i}",sig,guild,ya))
            s.append((f"s{guild}{i}",f"B{i}",sig,guild,yb))
    z=matched_cells(t,s);assert z["delta"]==2
    rows=[z|{"block":"1","donor":"1"},z|{"block":"2","donor":"2"}]
    a=coherent_genus_sensitivity(rows)
    assert np.allclose(a["groups"]["block"]["joint_genus_and_group_weight_interval"],[2,2])
    assert matched_cells(t[:3],s) is None
    assert canonical_rows([{"a":1},{"a":2}])==canonical_rows([{"a":2},{"a":1}])
    assert canonical_rows([{"a":1},{"a":1}])!=canonical_rows([{"a":1}])
    print("EXTERNAL_GIFT_EXPANSION_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--prior",type=pathlib.Path);p.add_argument("--comparison-pilot",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    old=prior(a.prior)
    cmp=load_comparison(a.comparison_pilot) if a.comparison_pilot else None
    comparison_audit=compare_pilots(old,cmp) if cmp is not None else None
    if comparison_audit is not None:
        print("EXTERNAL_GIFT_PRIOR_PILOT_COMPARISON",json.dumps(comparison_audit,allow_nan=False),flush=True)
    (reg,lst,tax,env,ov),metadata_pins,metadata_audit=metadata_snapshot(old,cmp)
    snapshot_class="exact_original_metadata_bytes" if metadata_pins==old["metadata_pins"] else "versioned_current_metadata_original_checklist_manifest_and_census_required"
    print("EXTERNAL_GIFT_EXPANSION_SNAPSHOT_CLASS",snapshot_class,flush=True)
    eligible,zid=pilot.eligible_lists(reg,lst,tax)
    original,census0=pilot.geo_plan(eligible,env,ov)
    assert census0==old["geography_census"] and len(original)==96
    pilot.LIMIT=LIMIT
    plan,census=pilot.geo_plan(eligible,env,ov)
    assert plan[:96]==original and len(plan)>96 and len(plan)<=LIMIT
    # Prior artifact has list URLs, not the full historical target/donor
    # pairing map. On drift, require exact historic file set and census but
    # never claim identical historic pairing solely from those checks.
    old_urls={v["source"]["url"] for v in old["source_files"]}
    planned_original_lists={str(eligible["Island"][i]["list_ID"]) for i,m,d,g in original}|{str(eligible["Mainland"][m]["list_ID"]) for i,m,d,g in original}
    from urllib.parse import urlparse,parse_qs
    old_list_ids={parse_qs(urlparse(u).query)["listid"][0] for u in old_urls}
    assert planned_original_lists==old_list_ids and len(old_list_ids)==155,"original planned checklist set differs"
    print("EXTERNAL_GIFT_EXPANSION_GEOGRAPHY",json.dumps(census),flush=True)
    tables=[]
    for name in ("guild_broad","self_fertilization"):
        tid,pin,n=PIN[name];t,d=trait_rows(tid,"strict")
        assert len(d["sources"])==1 and d["sources"][0]["sha256"]==pin and len(t)==n
        tables.append(t)
    prior_files={v["source"]["url"]:v for v in old["source_files"]};assert len(prior_files)==155
    old_ids={i for i,m,d,g in original}|{m for i,m,d,g in original}
    data={};sources={};errors=[];seen=set()
    for kind,ids in (("Island",sorted({i for i,m,d,g in plan})),("Mainland",sorted({m for i,m,d,g in plan}))):
        for i in ids:
            try:
                v,pin=pilot.roster(i,eligible[kind][i],zid)
                url=pin["source"]["url"]
                if i in old_ids:
                    assert pin==prior_files[url],"original checklist drift"
                    seen.add(url)
                data[i]=v;sources[i]=pin
            except Exception as e:
                if i in old_ids:raise
                errors.append({"entity_ID":i,"list_ID":eligible[kind][i]["list_ID"],"error_type":type(e).__name__,"detail":str(e)[:160]})
            time.sleep(.35)
    assert seen==set(prior_files)
    panels={}
    plans=dict(zip(COHORTS,(plan[:96],plan[96:],plan)))
    for scope in SCOPES:
        for label,z in plans.items():
            for mode,v in eval_cohort(z,data,tables,scope,label).items():
                k="__".join((scope,label,mode));panels[k]=v
                print("EXTERNAL_GIFT_EXPANSION_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
    for scope in SCOPES:
        p0=panels[scope+"__original_96__inclusive"];o=old["panels"][scope+"__common_reference_pattern"]
        assert p0["pairs"]==o["pairs"]
        if p0["pairs"]:assert np.isclose(p0["mean_island_equal_delta"],o["mean_island_equal_delta"],atol=1e-12)
    # Exhaustively leave out every publisher reference ID that appears in
    # the paired insect/wind trait records, not just an inconvenient source.
    A,B=tables;relevant={i for i in A.keys()&B.keys() if A[i]["v"] in GUILDS}
    ids=sorted(set().union(*(A[i]["refs"]|B[i]["refs"] for i in relevant)))
    leave={}
    for ref in ids:
        for label in ("added_targets","pooled"):
            v=eval_cohort(plans[label],data,tables,SCOPES[0],label,exclude_ref=ref)["inclusive"]
            k=label+"__exclude_ref_"+ref;leave[k]=v
            print("EXTERNAL_GIFT_EXPANSION_LEAVE_SOURCE",k,json.dumps(v,allow_nan=False),flush=True)
    out={"geography":census,"prior_artifact_sha256":PRIOR_SHA,
         "snapshot_class":snapshot_class,"comparison_pilot_sha256":NEW_PRIOR_SHA if cmp is not None else None,
         "comparison_pilot_audit":comparison_audit,"metadata_pins":metadata_pins,"metadata_audit":metadata_audit,
         "current_geography_plan":[{"target_entity":i,"comparison_entity":m,"distance_km":d,"block":list(g)} for i,m,d,g in plan],
         "original_checklist_hashes_verified":len(seen),"distinct_checklists_downloaded":len(data),
         "new_download_errors":errors,"source_manifest":list(sources.values()),
         "panels":panels,"leave_each_reference_out":leave,
         "warning":"Expansion of current metadata-fixed geographic pilot, not independent observations. Exact old metadata bytes required by default. Explicit comparison-pilot mode on metadata drift verifies two prior pilot manifests, historical checklist set and 155 checklist hashes, census and panel points; current canonical stability is not historical equivalence. Current 96-prefix preserved; historic exact target-donor pairing map not archived, so versioned run cannot prove its identity. At most 288 targets, no response-based substitution. Reference-conditioned aggregate self fertilization contrast; no verified autonomous selfing, actual effective pollination, ancestral source, or causal island effect. Added target islands can share taxa, comparators and publications. Disjoint removes focal native names within scope. Genus weights are sensitivity distributions."}
    assert len(panels)==18
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("EXTERNAL_GIFT_EXPANSION_COMPLETE",json.dumps({"panels":len(panels),"leaveout_panels":len(leave),"original_checklist_hashes_verified":len(seen),"distinct_checklists_downloaded":len(data),"new_download_errors":len(errors)}),flush=True)
if __name__=="__main__":main()
