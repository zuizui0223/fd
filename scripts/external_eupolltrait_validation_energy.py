#!/usr/bin/env python3
"""Replay the discovery half-energy estimand in the separate validation study cohort."""
import argparse,json,pathlib,tempfile
import numpy as np
import pandas as pd
from external_eupolltrait_size_turnover import load_traits
from external_eupolltrait_size_qc_family import corrected_traits,network
from external_eupolltrait_visitor_genus import genus_design
from external_eupolltrait_energy import compare,panel,self_test as energy_self_test
from external_eupolltrait_disjoint_validation import scope
from external_eupolltrait_provenance import components,cluster_infer
from external_eupolltrait_intake import prepare_network
from external_euppollnet_bee_loss import BASE,SOURCES,download

def component_metric(rows,mp,metric):
    # Same observation/null/gain fields and original weights; no new bootstrap.
    a=rows[["study","year","plant","site_a","site_b",
            metric+"_obs",metric+"_null",metric+"_gain"]].rename(columns={
        metric+"_obs":"size_obs",metric+"_null":"size_null",metric+"_gain":"size_gain"})
    return cluster_infer(a,mp)

def self_test():
    energy_self_test()
    x=pd.DataFrame({"study":["a","b"],"year":[2020,2020],
      "plant":["Genus species","Genus other"],"site_a":[(0,0)]*2,"site_b":[(0,1)]*2,
      "energy_half_obs":[.1,.2],"energy_half_null":[.3,.6],"energy_half_gain":[.2,.4]})
    z=component_metric(x,{"a":"a","b":"b"},"energy_half")
    assert abs(z["estimate"]-.3)<1e-12 and not z["gate"]
    print("FD_EUPOLLTRAIT_VALIDATION_ENERGY_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    tax,m,meta,pins=load_traits();tab,audit=corrected_traits(m,"bee")
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);npins={}
        for f,(sha,blob) in SOURCES.items():npins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
        raw=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False);assert len(raw)==623476
        md=pd.read_csv(w/"Metadata.csv",dtype=str)
        keys,sa=scope(md);mp,ma=components(md)
        disc={mp[s] for s in sa["discovery_study_ids"]}
        d,na=prepare_network(raw,md,eligible_study_year=keys)
        assert not set(d.Study_id)&set(sa["discovery_study_ids"])
        fm=tax.set_index("name").family.to_dict()
        q=d.loc[d.group.eq("bee")&d.name.isin(tab)&d.name.map(fm).eq(d.family)]
        names=sorted(q.name.unique());T=np.array([tab[k] for k in names])
        assert len(names)==501 and len(q)==438950
        D=np.abs(T[:,None]-T[None,:])
        groups,codes,E0,movable=genus_design(names,fm,T)
        E=E0[codes[:,None],codes[None,:]];np.fill_diagonal(E,0.)
        def pair_fn(a,b,sa,sb,tt,F,K):
            assert np.array_equal(T,tt)
            return compare(a,b,sa,sb,D,E)
        rows,diag=network(d,tab,"bee",tax,pair_fn=pair_fn)
        assert len(rows)==19621
        primary=panel(rows)
        assert primary["comparisons"]==19618
        assert abs(primary["study_equal"]["cross_gain"]-.10953097474723927)<1e-10
        sens=rows.loc[rows.min_exclusive.ge(2)]
        assert len(sens)==6220
        strict=rows.loc[~rows.study.map(mp).isin(disc)]
        component=component_metric(rows,mp,"energy_half")
        assert abs(component["estimate"]-primary["study_equal"]["energy_half_gain"])<1e-10
        out={"source_pins":pins,"network_source_pins":npins,"selection":audit,
             "scope_audit":sa,"network_name_audit":na,"metadata_components":ma,"network":diag,
             "primary":primary,"component_dependence":component,
             "descriptive_both_exclusive_ge2":panel(sens),
             "descriptive_both_exclusive_ge2_component":component_metric(sens,mp,"energy_half"),
             "descriptive_discovery_component_disjoint":component_metric(strict,mp,"energy_half"),
             "scope":"Post-discovery application of frozen univariate half-energy endpoint to separate study-ID validation cohort. Same fixed local degrees, exact shared count, 399 draws, original study/year/plant weights. Cross distance equals within-set dispersion plus half-energy. Genus-label expectation descriptive only. Discovery energy gate failed; a positive validation result does not erase that failure. No independent publisher replication, causal equivalence, multi-axis functional diversity or actual pollination claim."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
        print("FD_EUPOLLTRAIT_VALIDATION_ENERGY_RESULT",json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()
