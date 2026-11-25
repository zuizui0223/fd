#!/usr/bin/env python3
"""GIFT 3.1 direct publisher traits and bounded native-checklist pilot.

No project-created botanical data. Trait categories are publisher-coded
syndromes (possibly inferred), NOT observed effective pollinator identities.
Pilot lists are selected from metadata before trait values are read.
"""
import argparse,hashlib,itertools,json,math,pathlib,re,time,urllib.parse
from collections import Counter
from tropical_pollinator_source_preflight import get
from external_floral_syndrome_preflight import GIFT,tab

IDS={"guild_broad":"3.6.2","guild_detailed":"3.6.3","colour":"3.21.1",
     "flower_max_length_cm":"3.25.2","flower_max_width_cm":"3.26.2","self_fertilization":"3.1.1"}
SCOPES={"strict":(0,0),"inclusive":(1,1)}
GROUPS={"island":{"Island"},"mainland":{"Mainland"}}
SEED="fd_external_gift31_preflight_20261125"
N_LISTS=24

def fetch(query,**params):
    url=GIFT+"?"+urllib.parse.urlencode({"query":query,**params})
    raw=get(url,maxbytes=80_000_000)
    rows=tab(raw)
    return rows,{"url":url,"sha256":hashlib.sha256(raw).hexdigest(),"rows":len(rows)}

def num(v):
    try:
        f=float(v);return f if math.isfinite(f) else None
    except (TypeError,ValueError):return None

def ident(v):
    if v is None or not str(v).isdigit():raise ValueError("invalid publisher ID")
    return str(int(v))

def tokens(s):
    # Preserve the actual raw state census, interpret only exact recognized
    # single states. Compound strings stay unknown until delimiter audit.
    return str(s).strip().casefold()

def refs(v):
    return {x for x in re.split(r"[,;|\s]+",str(v or "")) if x.isdigit()}

def trait_rows(tid,scope):
    bias_ref,bias_deriv=SCOPES[scope]
    rows=[];sources=[]
    # Published traits API pagination, conservative cap. All intended
    # endpoints have fewer than 10,000 entries in the pinned metadata.
    for offset in range(0,50000,10000):
        x,pin=fetch("traits",traitid=tid,biasref=bias_ref,biasderiv=bias_deriv,startat=offset,limit=10000)
        sources.append(pin);rows.extend(x)
        if len(x)<10000:break
    else:raise ValueError("pagination cap reached")
    if not rows:return {},{"sources":sources,"status":"empty"}
    req={"work_ID","trait_value","agreement","references"}
    cols=set().union(*(r.keys() for r in rows));assert req<=cols
    grouped={}
    for r in rows:
        grouped.setdefault(ident(r["work_ID"]),[]).append(r)
    use={};duplicate=conflict=low=missing=0;values=Counter()
    numeric=tid in ("3.25.2","3.26.2")
    for i,rr in grouped.items():
        duplicate+=len(rr)-1
        vv={str(r["trait_value"]).strip() for r in rr}
        if len(vv)!=1:conflict+=1;continue
        v=next(iter(vv))
        if numeric:
            val=num(v)
            if val is None or val<=0:missing+=1;continue
        else:
            val=tokens(v)
            if val in ("","none","na","nan","null"):missing+=1;continue
            aa=[num(r["agreement"]) for r in rr]
            if any(a is not None and a<.66 for a in aa):low+=1;continue
            values[val]+=1
        use[i]={"v":val,"refs":set().union(*(refs(r["references"]) for r in rr))}
    return use,{"status":"retrieved","sources":sources,"columns":sorted(cols),
        "raw_rows":len(rows),"unique_work_ids":len(grouped),"eligible_work_ids":len(use),
        "duplicate_rows":duplicate,"conflicting_work_ids_excluded":conflict,
        "low_agreement_excluded":low,"invalid_value_excluded":missing,
        "categorical_states":dict(values) if not numeric else None,
        "numeric_min_max":[min(v["v"] for v in use.values()),max(v["v"] for v in use.values())] if numeric and use else None,
        "records_with_reference_ids":sum(bool(v["refs"]) for v in use.values())}

def pair_audit(tables):
    out={}
    for a,b in itertools.combinations(tables,2):
        A=tables[a];B=tables[b];common=A.keys()&B.keys()
        both_refs=[k for k in common if A[k]["refs"] and B[k]["refs"]]
        out[a+"__"+b]={"shared_work_ids":len(common),"both_reference_ids_available":len(both_refs),
            "shared_reference":sum(bool(A[k]["refs"]&B[k]["refs"]) for k in both_refs),
            "disjoint_references":sum(not (A[k]["refs"]&B[k]["refs"]) for k in both_refs)}
    return out

def pilot_selection(regions,lists,taxonomy):
    taxa={ident(x["taxon_ID"]):x for x in taxonomy}
    target=[x for x in taxonomy if str(x["taxon_name"]).casefold()=="angiospermae"]
    assert len(target)==1
    z=target[0];zid=ident(z["taxon_ID"])
    lft=float(z["lft"]);rgt=float(z["rgt"])
    # Taxon list must encompass all angiosperms, not just one family.
    complete={i for i,x in taxa.items() if num(x.get("lft")) is not None and num(x.get("rgt")) is not None and float(x["lft"])<=lft and float(x["rgt"])>=rgt}
    geo={ident(x["entity_ID"]):x for x in regions}
    out={};eligible={}
    for group,classes in GROUPS.items():
        roster={}
        for x in lists:
            i=ident(x["entity_ID"]);g=geo[i]
            if g.get("entity_class") not in classes:continue
            if str(x.get("restricted"))!="0" or str(x.get("suit_geo"))!="1":continue
            if str(x.get("native_indicated"))!="1" or ident(x["taxon_ID"]) not in complete:continue
            if str(x.get("subset","")).strip().casefold() not in {"all","native","native and naturalized","native and historically introduced"}:continue
            roster.setdefault(i,[]).append(x)
        eligible[group]={"entities":len(roster),"lists":sum(map(len,roster.values()))}
        # One checklist per entity, prefer literal native subset, stable IDs.
        picked=sorted(roster,key=lambda i:hashlib.sha256((SEED+group+i).encode()).hexdigest())[:N_LISTS]
        out[group]=[(i,sorted(roster[i],key=lambda x:(str(x.get("subset")).casefold()!="native",int(x["list_ID"])))[0]) for i in picked]
    return out,eligible,zid

def native_list(i,x,zid):
    rows,pin=fetch("checklists",listid=ident(x["list_ID"]),taxonid=zid,namesmatched=0,filter="native")
    names={};excluded=Counter()
    for r in rows:
        if ident(r["entity_ID"])!=i or ident(r["list_ID"])!=ident(x["list_ID"]):raise ValueError("list identity mismatch")
        # Literal documented nonquestionable native; never infer missing status.
        if str(r.get("native"))!="1":excluded["nonliteral_native"]+=1;continue
        if str(r.get("questionable"))!="0" or str(r.get("quest_native"))!="0":
            excluded["questionable_or_unknown"]+=1;continue
        if r.get("work_ID") in (None,""):excluded["no_work_id"]+=1;continue
        w=ident(r["work_ID"])
        g=r.get("genus_ID")
        names.setdefault(w,set()).add(str(g))
    clean={w:next(iter(g)) for w,g in names.items() if len(g)==1}
    return clean,{"source":pin,"native_work_ids":len(clean),"genus_conflict_excluded":len(names)-len(clean),"excluded_rows":dict(excluded)}

def list_coverage(sel,zid,tables):
    # Raw source retrieval is bounded and nonconcurrent. Coverage only;
    # no biological contrasts on a small, not spatially disjoint pilot.
    out={}
    for group,lists in sel.items():
        rows=[];errors=[]
        for i,x in lists:
            try:
                names,d=native_list(i,x,zid)
                n=set(names)
                d["trait_carrier_counts"]={k:len(n&v.keys()) for k,v in tables.items()}
                d["paired_carrier_counts"]={a+"__"+b:len(n&tables[a].keys()&tables[b].keys()) for a,b in itertools.combinations(tables,2)}
                rows.append(d)
            except Exception as e:
                errors.append({"list_ID":x["list_ID"],"error_type":type(e).__name__,"detail":str(e)[:180]})
            time.sleep(.2)
        cov={}
        for a,b in itertools.combinations(tables,2):
            k=a+"__"+b;v=[r["paired_carrier_counts"][k] for r in rows]
            cov[k]={"at_least_1":sum(x>=1 for x in v),"at_least_5":sum(x>=5 for x in v),"at_least_10":sum(x>=10 for x in v),"max":max(v,default=0)}
        out[group]={"selected_lists":len(lists),"success":len(rows),"errors":errors,"coverage":cov,"per_list_source_and_counts":rows}
    return out

def self_test():
    assert tokens(" Red ")=="red" and refs("12, 13;14")=={"12","13","14"}
    assert num("nan") is None and ident("123")=="123"
    a=pair_audit({"x":{"1":{"refs":{"4"}},"2":{"refs":{"5"}}},"y":{"1":{"refs":{"6"}},"2":{"refs":{"5"}}}})
    assert a["x__y"]=={"shared_work_ids":2,"both_reference_ids_available":2,"shared_reference":1,"disjoint_references":1}
    print("EXTERNAL_GIFT_TRAITS_TEST_PASS")

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    reg,rpin=fetch("regions");lst,lpin=fetch("lists");tax,tpin=fetch("taxonomy");meta,mpin=fetch("traits_meta")
    mt={r["Lvl3"]:r for r in meta}
    selected,census,zid=pilot_selection(reg,lst,tax)
    print("EXTERNAL_GIFT_LIST_ELIGIBILITY",json.dumps(census),flush=True)
    out={"metadata_sources":[rpin,lpin,tpin,mpin],"angiosperm_taxon_id":zid,"metadata_eligible_lists":census,"traits":{},"overlap":{}}
    tables={}
    for scope in SCOPES:
        ts={};out["traits"][scope]={}
        for k,tid in IDS.items():
            if tid not in mt:
                out["traits"][scope][k]={"status":"trait_id_absent"};continue
            # Check the actual official name for the self-fertilization ID.
            if k=="self_fertilization" and not re.search(r"self|fertiliz",str(mt[tid]["Trait2"]),re.I):
                out["traits"][scope][k]={"status":"unexpected_trait_definition","metadata":mt[tid]};continue
            try:
                ts[k],d=trait_rows(tid,scope)
                out["traits"][scope][k]=d|{"metadata":mt[tid]}
            except Exception as e:
                out["traits"][scope][k]={"status":"unavailable_or_invalid","error_type":type(e).__name__,"detail":str(e)[:180],"metadata":mt[tid]}
            print("EXTERNAL_GIFT_TRAIT",scope,k,json.dumps(out["traits"][scope][k],allow_nan=False),flush=True)
        out["overlap"][scope]=pair_audit(ts)
        tables[scope]=ts
        print("EXTERNAL_GIFT_OVERLAP",scope,json.dumps(out["overlap"][scope],allow_nan=False),flush=True)
    out["native_checklist_pilot"]=list_coverage(selected,zid,tables["strict"])
    out["warning"]="Pinned public publisher v3.1 intake; not latest version claim. Metadata selects one suitable unrestricted native-aiming list per pure island or pure mainland entity; no overlapping-entity exclusion yet. Metadata suitability is not proof of full floristic completeness. Work_ID overlap is not same individual or independently measured morphology/guild. Publisher syndrome labels can be inferred. Pilot is coverage only; no ecological effect or actual pollinator inference. All errors preserved."
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("EXTERNAL_GIFT_PILOT",json.dumps({k:{x:y for x,y in v.items() if x!="per_list_source_and_counts"} for k,v in out["native_checklist_pilot"].items()},allow_nan=False),flush=True)
if __name__=="__main__":main()
