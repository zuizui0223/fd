#!/usr/bin/env python3
"""Outcome-blind EMLN external spatial-pollination validation support audit."""
import argparse, collections, hashlib, io, itertools, json, pathlib, re, tempfile, urllib.parse, urllib.request
import pandas as pd
from external_euppollnet_bee_loss import BASE,SOURCES,download
from external_euppollnet_channel_retention import channel
from external_floral_replication_intake import name

ROOT=pathlib.Path(__file__).resolve().parents[1]
M=ROOT/"config/external_emln_spatial_source_manifest_20261128.json"
NA={"","na","nan","none","unknown","not recorded","unresolved"}
DOI=re.compile(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+",re.I)

def source_get(m,d,f):
    pin=m["directories"][d][f]
    url="https://raw.githubusercontent.com/"+m["source_repo"]+"/"+m["source_commit"]+"/website/data/data-raw/"+urllib.parse.quote(d)+"/"+f
    with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"FD source verification"}),timeout=45) as response:
        b=response.read(pin["size"]+1)
    assert len(b)==pin["size"],"publisher file size drift"
    assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==pin["blob"],"publisher blob drift"
    return pd.read_csv(io.BytesIO(b),dtype=str,keep_default_na=False),hashlib.sha256(b).hexdigest()

def doi_set(values):
    return {m.group(0).casefold().rstrip(".,;") for x in values for m in DOI.finditer(str(x))}

def pivot(d,key):
    assert {key,"attribute","value"}.issubset(d)
    out={}
    for k,a,v in d[[key,"attribute","value"]].itertuples(index=False,name=None):
        out.setdefault(k,{}).setdefault(a,set()).add(v)
    return out

def scalar(d,k):
    v=d.get(k,set())
    return next(iter(v)) if len(v)==1 else ""

def exact_node(d):
    if scalar(d,"taxa_verified").casefold()!="true" or scalar(d,"verification_level").casefold()!="species":return None
    s=scalar(d,"species")
    if s.casefold() in NA:s=scalar(d,"node_name")
    s=" ".join(s.replace("_"," ").split())
    if not s:return None
    # Publisher lower-case verified binomials permitted; do not strip extra tokens.
    return name(s[0].upper()+s[1:])

def schema_audit(n,l,e):
    N=pivot(n,"node_id");L=pivot(l,"layer");E=pivot(e,"interaction_id")
    lookup=collections.defaultdict(set)
    roles=collections.Counter();rank=collections.Counter();taxonomy=collections.Counter()
    strict={};trait_count=collections.Counter()
    for k,d in N.items():
        nm=scalar(d,"node_name");role=scalar(d,"type").casefold()
        if nm:lookup[nm].add(k)
        lookup[k].add(k)
        roles[role or "[unknown]"]+=1
        if role=="pollinator":
            rank[scalar(d,"verification_level") or "[missing]"]+=1
            order=scalar(d,"order").strip().casefold();family=scalar(d,"family").strip().casefold()
            if order not in NA:taxonomy["order_known"]+=1
            if order not in NA and family not in NA:taxonomy["order_family_known"]+=1
            for tr in ("proboscis_length","body_length","body_width","body_thickness","proboscis_width"):
                try:
                    if 0<float(scalar(d,tr))<float("inf"):trait_count[tr]+=1
                except ValueError:pass
        x=exact_node(d)
        if not x or role not in ("plant","pollinator"):continue
        if role=="pollinator" and (order in NA or family in NA):continue
        strict[k]=(role,x)
    layer_map=collections.defaultdict(set)
    for k,d in L.items():
        layer_map[k].add(k)
        v=scalar(d,"name")
        if v:layer_map[v].add(k)
    def resolve(v,m):
        s=m.get(v,set())
        return next(iter(s)) if len(s)==1 else None
    def layer_year(z):
        y=scalar(z,"year").strip()
        return y if re.fullmatch(r"\d{4}",y) else None
    counts=collections.Counter();by=collections.defaultdict(set)
    for d in E.values():
        a,b=resolve(scalar(d,"node_from"),lookup),resolve(scalar(d,"node_to"),lookup)
        x,y=resolve(scalar(d,"layer_from"),layer_map),resolve(scalar(d,"layer_to"),layer_map)
        if a is None or b is None or x is None or y is None:
            counts["unresolved_endpoint_or_layer"]+=1;continue
        if x!=y:
            counts["cross_layer"]+=1;continue
        if scalar(d,"type").casefold() not in ("pollination",""):
            counts["other_interaction_type"]+=1;continue
        try:
            if not 0<float(scalar(d,"weight"))<float("inf"):raise ValueError
        except ValueError:
            counts["invalid_weight"]+=1;continue
        ra,rb=scalar(N[a],"type").casefold(),scalar(N[b],"type").casefold()
        if {ra,rb}!={"plant","pollinator"}:
            counts["invalid_roles"]+=1;continue
        counts["role_resolved_positive_local_links"]+=1
        if a not in strict or b not in strict:
            counts["missing_verified_species_or_visitor_taxonomy"]+=1;continue
        p,v=(a,b) if ra=="plant" else (b,a)
        by[x].add((strict[p][1],strict[v][1]))
        counts["strict_local_links"]+=1
    # Distinct layers are only candidate sites; actual independent spatial
    # replication and comparable protocol are manual promotion gates.
    eligible={}
    for k,edges in by.items():
        plants=collections.defaultdict(set)
        for p,v in edges:plants[p].add(v)
        if len(plants)>=5 and len({v for p,v in edges})>=5:
            eligible[k]=plants
    pairs=0;shared=0;turnover=0
    for a,b in itertools.combinations(sorted(eligible),2):
        ya,yb=layer_year(L[a]),layer_year(L[b])
        if ya and yb and ya!=yb:continue
        pairs+=1
        A,B=eligible[a],eligible[b]
        for p in A.keys()&B.keys():
            shared+=1
            if A[p]-B[p] and B[p]-A[p]:turnover+=1
    return {"nodes":len(N),"layers":len(L),"interaction_records":len(E),"role_counts":dict(roles),
            "visitor_verification_levels":dict(rank),"visitor_taxonomy":dict(taxonomy),"positive_measured_trait_carriers":dict(trait_count),
            "strict_nodes":dict(collections.Counter(x[0] for x in strict.values())),
            "layer_type_counts":dict(collections.Counter(v for d in L.values() for v in d.get("type",set()))),
            "explicit_year_layers":sum(layer_year(d) is not None for d in L.values()),
            "coordinate_field_layers":sum(bool(scalar(d,"latitude")) and bool(scalar(d,"longitude")) for d in L.values()),
            "edge_attrition":dict(counts),"candidate_eligible_layers":len(eligible),
            "candidate_layer_pairs":pairs,"candidate_shared_plant_layer_pairs":shared,
            "candidate_two_sided_turnover_pairs":turnover,
            "promotion_requires_review":True}

def self_test():
    assert doi_set(["https://doi.org/10.1234/ABC."])=={"10.1234/abc"}
    assert exact_node({"taxa_verified":{"TRUE"},"verification_level":{"species"},"species":{"plant_species"}})=="plant species"
    assert exact_node({"taxa_verified":{"FALSE"},"verification_level":{"species"},"species":{"Plant_species"}}) is None
    assert scalar({"type":{"space","time"}},"type")==""
    assert pivot(pd.DataFrame({"k":["1","1"],"attribute":["type","type"],"value":["space","time"]}),"k")["1"]["type"]=={"space","time"}
    print("FD_EXTERNAL_EMLN_AUDIT_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true")
    a=p.parse_args();self_test()
    if a.self_test:return
    m=json.loads(M.read_text());out={};desc_counts=collections.Counter()
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);f="Metadata.csv";sha,blob=SOURCES[f]
        mpin=download(BASE+f,w/f,sha,5_000_000,blob)
        eu=pd.read_csv(w/f,dtype=str,keep_default_na=False)
        assert "DOI" in eu.columns or "doi" in eu.columns or any("doi" in c.casefold() for c in eu)
        col=next(c for c in eu if "doi" in c.casefold())
        eud=doi_set(eu[col])
    for d,files in sorted(m["directories"].items()):
        entry={"status":"unread","files":{}}
        try:
            meta,h=source_get(m,d,"description.csv");entry["files"]["description.csv"]=h
            assert {"attribute","value"}.issubset(meta)
            labels=set(meta.loc[meta.attribute.eq("ecological_network_type"),"value"].str.casefold())
            types=set(meta.loc[meta.attribute.eq("multilayer_network_type"),"value"].str.casefold())
            entry.update(ecological_network_type=sorted(labels),layer_types=sorted(types))
            if labels!={"pollination"} or "spatial" not in types:
                entry["status"]="not_spatial_pollination";out[d]=entry;continue
            refs,h=source_get(m,d,"references.csv");entry["files"]["references.csv"]=h
            ds=doi_set(refs.to_numpy().ravel())
            entry.update(source_dois=sorted(ds),exact_euppollnet_doi_overlap=sorted(ds&eud),no_reference_doi=not bool(ds))
            tables=[]
            for f in ("nodes.csv","layers.csv","interactions.csv"):
                t,h=source_get(m,d,f);entry["files"][f]=h;tables.append(t)
            entry["audit"]=schema_audit(*tables)
            entry["status"]="coverage_audited"
        except Exception as ex:
            entry.update(status="source_or_schema_failure",error=type(ex).__name__+": "+str(ex)[:180])
        out[d]=entry;print("FD_EXTERNAL_EMLN_AUDIT_PANEL",d,json.dumps(entry,allow_nan=False),flush=True)
    for v in out.values():desc_counts[v["status"]]+=1
    result={"source_commit":m["source_commit"],"euppollnet_metadata_pin":mpin,"euppollnet_doi_count":len(eud),
        "status_counts":dict(desc_counts),"panels":out,
        "scope":"Outcome-blind source/coverage audit only. Published layer identities are candidates, not proof of independent sites or comparable protocol; exact DOI nonoverlap not proof of observational independence. Taxonomic role direction, names, source geography and observation windows require original-paper verification. No ecological effect calculated, no own-data input, no raw interaction export."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False))
    print("FD_EXTERNAL_EMLN_AUDIT_COMPLETE",json.dumps({"status_counts":dict(desc_counts),"directories":len(out)}),flush=True)
if __name__=="__main__":main()
