#!/usr/bin/env python3
"""Immutable public LepTraits consensus schema and exact butterfly name linkage audit."""
import io,json,hashlib,pathlib,re
import pandas as pd
from tropical_pollinator_source_preflight import get,blobsha
URL="https://raw.githubusercontent.com/RiesLabGU/LepTraits/f73317d28b38f80923ad8e02cc06d35f5c16ed1f/consensus/consensus.csv"
BLOB="b2bc7bd87eca199b62f08ac539660aebc9fbaf79"
def norm(s):
    s=" ".join(str(s).replace("_"," ").strip().casefold().split())
    return s if re.fullmatch(r"[a-z-]+ [a-z-]+",s) else None
def audit_traits(out,inc):
    raw=get(URL,maxbytes=5_000_000);assert blobsha(raw)==BLOB
    d=pd.read_csv(io.BytesIO(raw),dtype=str,encoding="utf-8-sig",encoding_errors="strict")
    fields=list(d.columns)
    name_candidates={}
    observed={n for n in inc.final_species_name.map(norm) if n}
    for col in fields:
        if re.search("species|taxon|scientific|binomial|name",col,re.I):
            names={n for n in d[col].map(norm) if n}
            name_candidates[col]={"binomial_strings":len(names),"exact_overlap":len(names&observed)}
    # Do not choose a name column based on ecological outcomes; report all candidate overlap. No fuzzy/genus imputation.
    wing_fields=[c for c in fields if re.search("wing|forewing|fwlength",c,re.I)]
    host_fields=[c for c in fields if re.search("host",c,re.I)]
    res={"sha256":hashlib.sha256(raw).hexdigest(),"rows":len(d),"columns":fields,"name_candidates":name_candidates,"wing_field_nonmissing":{c:int(d[c].notna().sum()) for c in wing_fields},"host_fields":host_fields,"positive_incidence_rows":len(inc),"distinct_normalized_names":len(observed),"warning":"Schema and exact-name opportunity only. Confirm taxonomic column and measurement units before trait inference."}
    (out/"leptraits_linkage_audit.json").write_text(json.dumps(res,indent=2,ensure_ascii=False))
    print("LEPTRAITS_LINKAGE",json.dumps(res,ensure_ascii=False),flush=True)
