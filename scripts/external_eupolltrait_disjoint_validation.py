#!/usr/bin/env python3
"""Study-ID-disjoint bee-inclusive protocol validation of corrected ITD result."""
import argparse, collections, json, pathlib, re, tempfile
import numpy as np
import pandas as pd
from external_eupolltrait_size_turnover import load_traits,report
from external_eupolltrait_size_qc_family import corrected_traits,network,decomposition_report
from external_eupolltrait_intake import prepare_network
from external_euppollnet_bee_loss import BASE,SOURCES,download
from external_emln_spatial_audit import doi_set

# Exact publisher labels inspected before accessing these groups' ITD outcomes.
ALLOWED={
    "bees","bees, syrphids, conopids, tachinid flies and butterflies",
    "bees and syrphids","hymenoptera, syrphids, bombylids and coleoptera",
    "hymenoptera","bees, syrphids and butterflies",
    "hymenoptera, diptera and lepidoptera","hymenoptera and lpidoptera"
}
ORDERS=("hymenoptera","diptera","lepidoptera","coleoptera")

def scope(meta):
    m=meta.copy()
    m["year"]=pd.to_numeric(m.Year,errors="coerce")
    m["label"]=m.Taxa_recorded.fillna("").str.strip().str.casefold()
    m["broad"]=m.label.map(lambda x:all(k in x for k in ORDERS))
    m["allowed"]=m.label.isin(ALLOWED)
    # Exclude whole study if any metadata row qualified for discovery.
    discovery=set(m.loc[m.broad,"Study_id"])
    do=doi_set(m.loc[m.Study_id.isin(discovery),"DOI"])
    overlap=set()
    for study,g in m.groupby("Study_id"):
        if study not in discovery and doi_set(g.DOI)&do:overlap.add(study)
    # Reject inconsistent or missing metadata in each study-year.
    keys=m.groupby(["Study_id","year"],as_index=False).allowed.all()
    keys=keys.loc[keys.allowed&~keys.Study_id.isin(discovery|overlap),["Study_id","year"]]
    assert not set(keys.Study_id)&(discovery|overlap)
    a=m.merge(keys,on=["Study_id","year"],how="inner")
    # A citation DOI overlap is checked exactly; same projects or source
    # observations may still appear under other DOIs or missing references.
    suffix={s.split("_",1)[-1].casefold() for s in discovery}
    same_suffix={s for s in keys.Study_id.unique() if s.split("_",1)[-1].casefold() in suffix}
    # Same-study-name suffix is a conservative coarse grouping, not authorship.
    audit={"discovery_metadata_studies":len(discovery),"discovery_study_ids":sorted(discovery),
       "exact_doi_overlap_excluded_study_ids":sorted(overlap),
       "validation_candidate_studies":int(keys.Study_id.nunique()),"validation_study_years":len(keys),
       "validation_study_ids":sorted(keys.Study_id.unique()),
       "validation_labels":{str(k):int(v) for k,v in a.label.value_counts().items()},
       "shared_discovery_study_suffix_ids":sorted(same_suffix),
       "excluded_scope_labels":sorted(set(m.loc[~m.broad&~m.allowed,"label"])),
       "limitations":"Study ID/DOI disjointness, not certified independent field observations; scope labels document intended taxa, not complete detection."}
    return keys,audit

def panel(rows):
    return {"size":report(rows),
            "family_accounting_descriptive":decomposition_report(rows),
            "both_exclusive_ge2_descriptive":report(rows.loc[rows.min_exclusive.ge(2)] if not rows.empty else rows)}

def self_test():
    m=pd.DataFrame({"Study_id":["a","a","b","c","d","e"],"Year":["2020"]*6,
        "Taxa_recorded":["Hymenoptera, Diptera, Lepidoptera and Coleoptera","Bees",
        "Bees","Bumblebees","Bees and syrphids","Hymenoptera"],
        "DOI":["10.1234/a","10.1234/a","10.1234/a","NA","10.1234/d","Unpublished"]})
    k,a=scope(m)
    assert set(k.Study_id)=={"d","e"} and a["exact_doi_overlap_excluded_study_ids"]==["b"]
    print("FD_EUPOLLTRAIT_DISJOINT_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    tax,m,meta,pins=load_traits();tab,selection=corrected_traits(m,"bee")
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);npins={}
        for f,(sha,blob) in SOURCES.items():npins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
        raw=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False);assert len(raw)==623476
        md=pd.read_csv(w/"Metadata.csv",dtype=str)
        keys,audit=scope(md);print("FD_EUPOLLTRAIT_DISJOINT_SCOPE",json.dumps(audit,allow_nan=False),flush=True)
        d,na=prepare_network(raw,md)
        annual,ad=network(d,tab,"bee",tax);ar=report(annual)
        assert len(annual)==1261 and abs(ar["study_equal"]["size_gain"]-.11064886953661242)<1e-10
        v,va=prepare_network(raw,md,eligible_study_year=keys)
        assert not set(d.Study_id)&set(v.Study_id)
        rows,diag=network(v,tab,"bee",tax)
        main=panel(rows)
        # Exclusion sensitivity only, not an independent confirmation.
        excl=set(audit["shared_discovery_study_suffix_ids"])
        s=rows.loc[~rows.study.isin(excl)] if not rows.empty else rows
        out={"source_pins":pins,"network_source_pins":npins,"selection":selection,
          "scope_audit":audit,"network_name_audit":va,
          "annual_discovery_replay":{"comparisons":len(annual),"gain":ar["study_equal"]["size_gain"]},
          "validation_network":diag,"primary":main,
          "descriptive_exclude_shared_study_suffix":report(s),
          "scope":"Post-discovery study-ID-disjoint protocol validation in SAME publisher database. Exact whitelist of broad bee-inclusive metadata, excludes any discovery study and exact overlapping DOI; no independent source guarantee or fully blinded claim. Same corrected ITD, fixed degrees/shared counts, original crossed uncertainty and frozen gate. Family and >=2-exclusive sensitivities descriptive. No actual pollination or island claim."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
        print("FD_EUPOLLTRAIT_DISJOINT_RESULT",json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()
