"""Outcome-blind MDD v2.5 nomenclatural coverage audit.

A current accepted taxon can connect an older map and diet name, but a split
or lump is not proof of equal traits. Only reciprocal singleton connections
with an unambiguous accepted species ID are allowed in a separate sensitivity.
No fuzzy matching, no congener trait imputation, no transfer across >1 map or
diet name in an accepted ID. Keep exact matches as the original primary.
"""
from collections import defaultdict,Counter
import hashlib,io,json,re
import pandas as pd
from tropical_pollinator_source_preflight import get

BASE="https://zenodo.org/records/21654811/files/"
FILES={"Species_Syn_v2.5.csv":("018ab2df0e2f18b7aa032965511e3e11",80_000_000),
       "MDD_v2.5_6904species.csv":("533e662fd5b8f66a5f56c191e7efac44",12_000_000)}

def norm(x):return " ".join(str(x).replace("_"," ").split()).casefold()
def binomial(x):
    s=norm(x)
    return s if re.fullmatch(r"[a-z][a-z-]+ [a-z][a-z-]+",s) else None

def key(x):return re.sub("[^a-z0-9]","",x.casefold())
def read(name):
    md5,cap=FILES[name];raw=get(BASE+name+"?download=1",maxbytes=cap)
    assert hashlib.md5(raw).hexdigest()==md5
    return pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False,low_memory=False),{"md5":md5,"sha256":hashlib.sha256(raw).hexdigest()}

def col(d,*choices):
    m={key(x):x for x in d.columns}
    for x in choices:
        if key(x) in m:return m[key(x)]
    raise KeyError("No expected column "+str(choices)+" in "+str(list(d.columns)))

def index_names(checklist,syn):
    cid=col(checklist,"id","MDD_species_id","MDD_id")
    cn=col(checklist,"sciName","MDD_species","scientific_name")
    order=col(checklist,"order","MDD_order")
    bats=checklist[checklist[order].str.casefold().eq("chiroptera")]
    assert len(bats)>1000 and bats[cid].is_unique
    ids=set(bats[cid])
    acc={i:binomial(n) for i,n in zip(bats[cid],bats[cn])}
    assert all(acc.values()) and len(set(acc.values()))==len(acc)
    si=col(syn,"MDD_species_id")
    orig=col(syn,"MDD_original_combination")
    rank=col(syn,"MDD_original_rank")
    validity=col(syn,"MDD_validity")
    status=col(syn,"MDD_nomenclature_status")
    idx=defaultdict(set)
    for i,n in acc.items():idx[n].add(i)
    s=syn[syn[si].isin(ids)].copy()
    stats=Counter()
    for i,n,r,v,t in zip(s[si],s[orig],s[rank],s[validity],s[status]):
        stats["bat_synonym_rows"]+=1
        # Names originally described at species rank only. A subspecies,
        # hybrid, uncertain or unavailable name is not an automatic donor.
        if r.strip().casefold()!="species":continue
        if re.search(r"dubi|inquir|hybrid|composit|uncertain|misappl",v,re.I):continue
        if re.search(r"unavail|nomen nud|homonym|misappl",t,re.I):continue
        b=binomial(n)
        if b:idx[b].add(i);stats["eligible_original_binomial_rows"]+=1
    return idx,dict(stats)|{"accepted_bat_species":len(acc),"indexed_unique_names":len(idx),"ambiguous_names":sum(len(v)>1 for v in idx.values()),
        "validity_values":s[validity].value_counts().to_dict(),"nomenclature_values":s[status].value_counts().to_dict()}

def singleton_crosswalk(map_names,diet_names,idx):
    maps=set(map_names);diet=set(diet_names)
    def reverse(names):
        out=defaultdict(set)
        for n in names:
            for i in idx.get(n,set()):out[i].add(n)
        return out
    m=reverse(maps);d=reverse(diet)
    matches={};blocked=Counter()
    for n in sorted(maps-diet):
        ids=idx.get(n,set())
        if len(ids)!=1:
            blocked["no_or_ambiguous_accepted_id"]+=1;continue
        i=next(iter(ids))
        if m[i]!={n}:
            blocked["multiple_map_names"]+=1;continue
        if len(d[i])!=1:
            blocked["no_or_multiple_diet_names"]+=1;continue
        donor=next(iter(d[i]))
        if donor in maps or idx.get(donor,set())!={i}:
            blocked["donor_already_used_or_ambiguous"]+=1;continue
        matches[n]=donor
    assert len(set(matches.values()))==len(matches)
    assert not set(matches)&diet and not set(matches.values())&maps
    return matches,dict(blocked)

def build(map_names,diet,return_taxonomy=False):
    c,cm=read("MDD_v2.5_6904species.csv")
    s,sm=read("Species_Syn_v2.5.csv")
    idx,stats=index_names(c,s)
    matches,blocked=singleton_crosswalk(map_names,diet.index,idx)
    extra=diet.loc[list(matches.values())].copy()
    extra.index=list(matches)
    augmented=pd.concat([diet,extra])
    assert augmented.index.is_unique
    report={"sources":{"checklist":cm,"synonyms":sm},"taxonomy":stats,
            "exact_matches":len(set(map_names)&set(diet.index)),"reciprocal_singleton_added":len(matches),
            "remaining_unmatched":len(set(map_names)-set(augmented.index)),"blocked":blocked,
            "warning":"Nomenclatural sensitivity only; singleton does not certify historical taxonomic circumscription or measured donor trait. No fuzzy or genus-based trait imputation."}
    if return_taxonomy:
        from bat_family_local_completeness import taxonomy
        family_map,family_audit=taxonomy(c,s,set(map_names),augmented)
        return augmented,report,family_map,family_audit
    return augmented,report

def self_test():
    idx={"a a":{"1"},"b b":{"1"},"c c":{"2"},"d d":{"2"},"e e":{"2"},"f f":{"3","4"},"g g":{"3"},"h h":{"5"},"i i":{"5"}}
    result,_=singleton_crosswalk(["a a","c c","e e","f f","h h","i i"],["b b","d d","g g","h h"],idx)
    assert result=={"a a":"b b"}
    result,_=singleton_crosswalk(["a a"],["b b","a a"],idx)
    assert result=={}
    print("MDD_SINGLETON_CROSSWALK_TEST_PASS",flush=True)
if __name__=="__main__":self_test()
