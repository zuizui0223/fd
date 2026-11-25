#!/usr/bin/env python3
"""Audit external GIFT coverage loss and trait-source concentration.

Do not interpret unknown checklist flags as confirmed absences.
Same metadata-defined sample as the first official GIFT pilot.
No project-created data or biological effect estimation.
"""
import argparse,collections,json,pathlib,time
from external_gift_floral_trait_overlap import fetch,ident,pilot_selection,trait_rows

PIN={
 "guild_broad":("3.6.2","66cfa0398c3222131e642f083e567651efc1719e04b15a18c1462911988e6b75",5852),
 "guild_detailed":("3.6.3","29e088a18c738c761a0d41f67966c4b455bbe3e9d7bf792b7eecac440b34ee52",3119),
 "colour":("3.21.1","07db04eaec9c67ffb6084591e07e75be1dd942a59bae58c3620943847bc07d66",729),
 "flower_length":("3.25.2","d7d2e1918ca70c578c215fbccf0f8ea9f466a22f256558b9e717f56445f6d114",1866),
 "self_fertilization":("3.1.1","59d44b73c6267fb6bb08a677b58395c5480392415a50b5eccbfd9e111407fd17",3530)}
PAIRS=(("guild_broad","self_fertilization"),("guild_detailed","self_fertilization"),("guild_detailed","colour"),("colour","flower_length"))

def reference_index(rows):
    grouped=collections.defaultdict(list)
    for r in rows:grouped[ident(r["ref_ID"])].append(r)
    fields=("ref_long","ref_short","DOI","doi")
    meta={};conflicts=0;multiline=0
    for k,v in grouped.items():
        if len(v)>1:multiline+=1
        unique={}
        for field in fields:
            vv=sorted({str(r[field]).strip() for r in v if r.get(field) not in (None,"") and str(r[field]).strip()})
            if vv:unique[field]=vv
        ambiguous=any(len(x)>1 for x in unique.values())
        conflicts+=int(ambiguous)
        meta[k]={"ref_ID":k,"metadata_rows":len(v),"citation_field_candidates":unique,
                 "citation_conflict":ambiguous}
    audit={"raw_reference_rows":len(rows),"distinct_reference_ids":len(grouped),
           "ids_with_multiple_rows":multiline,"ids_with_conflicting_citation_fields":conflicts,
           "maximum_rows_per_id":max(map(len,grouped.values()),default=0),
           "warning":"Reference IDs are publisher records, not independent studies. Repeated IDs can have several metadata rows; citation alternatives are retained, not arbitrarily selected."}
    return meta,audit

def concentration(ids,meta):
    counts=collections.Counter()
    for rr in ids:
        if rr:
            for r in rr:counts[r]+=1/len(rr)
    n=sum(counts.values())
    order=counts.most_common()
    return {"reference_count":len(counts),"species_equivalent_weight":n,
            "effective_reference_count":1/sum((c/n)**2 for r,c in order) if n else None,
            "largest_reference_share":order[0][1]/n if n else None,
            "reference_ids_missing_from_metadata":sum(r not in meta for r in counts),
            "top_references":[{"ref_ID":r,"fractional_species_weight":c,
                "share":c/n,"citation":meta.get(r)} for r,c in order[:10]]}

def source_pairs(tables,meta):
    out={}
    for a,b in PAIRS:
        A=tables[a];B=tables[b];ids=A.keys()&B.keys()
        states=collections.Counter((str(A[k]["v"]),str(B[k]["v"])) for k in ids) if b!="flower_length" else collections.Counter((str(A[k]["v"]),"numeric_length") for k in ids)
        shared=[A[k]["refs"]&B[k]["refs"] for k in ids]
        disjoint=[k for k in ids if A[k]["refs"] and B[k]["refs"] and not A[k]["refs"]&B[k]["refs"]]
        out[a+"__"+b]={"paired_species":len(ids),"state_counts":[{"axis_a":k[0],"axis_b":k[1],"species":v} for k,v in sorted(states.items())],
            "axis_a_source_concentration":concentration([A[k]["refs"] for k in ids],meta),
            "axis_b_source_concentration":concentration([B[k]["refs"] for k in ids],meta),
            "shared_source_concentration":concentration(shared,meta),
            "both_source_sets_disjoint_species":len(disjoint),
            "disjoint_axis_a_sources":concentration([A[k]["refs"] for k in disjoint],meta),
            "disjoint_axis_b_sources":concentration([B[k]["refs"] for k in disjoint],meta)}
    return out

def flag(v):
    if v is None or str(v).strip()=="":return "missing"
    return str(v).strip()

def row_sets(rows):
    cats=collections.Counter();rosters={k:set() for k in ("all_named","literal_native","native_no_positive_question_flags","strict_native")}
    for r in rows:
        w=r.get("work_ID")
        if w in (None,""):continue
        w=ident(w)
        n=flag(r.get("native"));q=flag(r.get("questionable"));s=flag(r.get("quest_native"))
        cats[(n,q,s)]+=1
        rosters["all_named"].add(w)
        if n!="1":continue
        rosters["literal_native"].add(w)
        if q=="1" or s=="1":continue
        if q not in {"0","missing"} or s not in {"0","missing"}:continue
        rosters["native_no_positive_question_flags"].add(w)
        if q=="0" and s=="0":rosters["strict_native"].add(w)
    return rosters,[{"native":k[0],"questionable":k[1],"quest_native":k[2],"rows":v} for k,v in sorted(cats.items())]

def checklists(selected,taxon,tables):
    out={}
    for group,ls in selected.items():
        successes=[];errors=[]
        for i,(e,l) in enumerate(ls):
            try:
                args=dict(listid=ident(l["list_ID"]),taxonid=taxon,namesmatched=0)
                rows,pin=fetch("checklists",**args,filter="native")
                assert all(ident(r["entity_ID"])==e and ident(r["list_ID"])==args["listid"] for r in rows)
                sets,flags=row_sets(rows)
                rec={"source":pin,"entity_class":l["entity_class"],"subset":l["subset"],
                    "flag_census":flags,"native_filter_rows":len(rows),
                    "counts":{s:{"named_species":len(v),"pair_carriers":{a+"__"+b:len(v&tables[a].keys()&tables[b].keys()) for a,b in PAIRS}} for s,v in sets.items()}}
                # Fixed first six lists per group also retrieve no filter.
                # Audit filter membership; do not replace the native cohort
                # by all named records even if coverage looks better.
                if i<6:
                    allrows,apin=fetch("checklists",**args)
                    assert all(ident(r["entity_ID"])==e and ident(r["list_ID"])==args["listid"] for r in allrows)
                    alls,aflags=row_sets(allrows)
                    rec["unfiltered_check"]={"source":apin,"rows":len(allrows),"flag_census":aflags,
                        "native_filter_names_missing_from_unfiltered":len(sets["all_named"]-alls["all_named"]),
                        "unfiltered_literal_native_not_in_filter":len(alls["literal_native"]-sets["all_named"])}
                successes.append(rec)
            except Exception as e:errors.append({"list_ID":l["list_ID"],"type":type(e).__name__,"detail":str(e)[:180]})
            time.sleep(.2)
        agg={}
        for s in ("all_named","literal_native","native_no_positive_question_flags","strict_native"):
            a0=[r["counts"][s] for r in successes]
            agg[s]={"nonempty_lists":sum(v["named_species"]>0 for v in a0),"total_names_counted_across_lists":sum(v["named_species"] for v in a0),
                 "coverage":{a+"__"+b:{str(n):sum(v["pair_carriers"][a+"__"+b]>=n for v in a0) for n in (1,5,10,30)} for a,b in PAIRS}}
        out[group]={"selected":len(ls),"success":len(successes),"errors":errors,"summary":agg,"per_list":successes}
        print("EXTERNAL_GIFT_FLAG_AUDIT",group,json.dumps(out[group],allow_nan=False),flush=True)
    return out

def self_test():
    a=[{"work_ID":"1","native":"1","questionable":"0","quest_native":None},
       {"work_ID":"2","native":"1","questionable":"0","quest_native":"0"},
       {"work_ID":"3","native":"1","questionable":"1","quest_native":"0"},
       {"work_ID":"4","native":None,"questionable":"0","quest_native":"0"}]
    v,f=row_sets(a)
    assert v["strict_native"]=={"2"} and v["native_no_positive_question_flags"]=={"1","2"} and v["literal_native"]=={"1","2","3"}
    assert concentration([{"1","2"},{"1"}],{})["largest_reference_share"]==.75
    m,a=reference_index([{"ref_ID":"1","ref_long":"A"},{"ref_ID":"1","ref_long":"B"},{"ref_ID":"2","ref_long":"C"},{"ref_ID":"2","ref_long":"C"}])
    assert a["distinct_reference_ids"]==2 and a["ids_with_multiple_rows"]==2 and a["ids_with_conflicting_citation_fields"]==1
    assert m["1"]["citation_field_candidates"]["ref_long"]==["A","B"] and not m["2"]["citation_conflict"]
    print("EXTERNAL_GIFT_COVERAGE_AUDIT_TEST_PASS")

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    regions,rpin=fetch("regions");lists,lpin=fetch("lists");tax,tpin=fetch("taxonomy");rr,refpin=fetch("references")
    assert rpin["sha256"]=="d12d9129eb8aeb7aafb295ddf99b3ffdd810dcd628ab951ca31e560ab772b64b"
    assert lpin["sha256"]=="7cbf01e767f5abc061fb3622db031c95e4b24646b98d875eed9f1d9e742ba244"
    meta,reference_audit=reference_index(rr)
    print("EXTERNAL_GIFT_REFERENCE_ID_AUDIT",json.dumps(reference_audit,ensure_ascii=False,allow_nan=False),flush=True)
    selected,census,zid=pilot_selection(regions,lists,tax)
    assert census=={"island":{"entities":1253,"lists":1496},"mainland":{"entities":751,"lists":810}}
    tables={}
    for name,(tid,pin,n) in PIN.items():
        vals,d=trait_rows(tid,"strict")
        assert len(d["sources"])==1 and d["sources"][0]["sha256"]==pin and len(vals)==n
        tables[name]=vals
    pairs=source_pairs(tables,meta)
    print("EXTERNAL_GIFT_REFERENCE_CONCENTRATION",json.dumps(pairs,ensure_ascii=False,allow_nan=False),flush=True)
    out={"metadata_pins":[rpin,lpin,tpin,refpin],"trait_pins":PIN,
        "reference_id_audit":reference_audit,"source_pair_concentration":pairs,"native_flag_audit":checklists(selected,zid,tables),
        "warning":"Official external-source audit only. Fractional reference weights describe concentration, not independent studies or causal evidence. Unknown questionable flags never interpreted as confirmed clean native status. Unfiltered and loose cohorts are upper support diagnostics, not inferred complete checklists. No ecological effect estimation."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False))
if __name__=="__main__":main()
