"""Independent original EltonTraits bird diet x pre-N1 GloBI guild-catalog audit.

Outputs COUNTS ONLY, never public taxon/source/plant/island identity joins.
- Wilman et al. 2014 EltonTraits BirdFuncDat original public figshare file.
- 2026-06-04 GloBI v0.9 species-resolved flower-visitor NEE catalog,
  assembled before focal island N1 and independently of focal floral traits.
This is source validity, NOT island-level guild occupancy, effectiveness,
phenotypic evolution, or a bird-to-moth pollinator switch.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ELTON_SOURCE_SHA = "97216eb1797da077169ebb1ebea275db293b09fc62f8bb8911f9beb98c50d321"
ELTON_WEIGHTS_SHA = "ebdeb7d03669e690bcaf088165f37cb47ba7ee2cd433462a31e120d992dfa63e"
GLOBI_CATALOG_ZIP_SHA = "c946e36cd04cb319fc63a5bde90534c1cf3c7f41480fc672638f85f5b2a89b02"
GLOBI_V09_SOURCE_SHA = "095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17"


def _read_member(path: Path, member: str) -> bytes:
    with ZipFile(path) as z:
        if z.namelist().count(member) != 1:
            raise ValueError(f"Missing/duplicate expected source member: {member}")
        return z.read(member)


def _read_sources(elton_zip: Path, globi_zip: Path):
    eb = _read_member(elton_zip, "eltontraits_bird_nectar_weights.csv")
    receipt = json.loads(_read_member(elton_zip, "source_receipt.json"))
    assert receipt["contract"] == "fd_original_eltontraits_bird_nectar_v1"
    if receipt["original_sha256"] != ELTON_SOURCE_SHA:
        raise ValueError("Frozen original EltonTraits public source SHA changed")
    if hashlib.sha256(eb).hexdigest() != ELTON_WEIGHTS_SHA:
        raise ValueError("Frozen EltonTraits processed dietary fractions changed")
    if receipt["n_name_unique_valid_diet_species"] != 9993:
        raise ValueError("EltonTraits 9993 accepted dietary species not recovered")
    if receipt["has_eBird_taxonomy_or_occurrences"] is not False:
        raise ValueError("Unexpected eBird source lineage")
    if hashlib.sha256(globi_zip.read_bytes()).hexdigest() != GLOBI_CATALOG_ZIP_SHA:
        raise ValueError("Frozen original source GloBI NEE artifact SHA changed")
    gb = _read_member(globi_zip, "channel_taxon_catalog_confirmatory.csv")
    gc = json.loads(_read_member(globi_zip, "channel_taxon_catalog_receipt.json"))
    if (gc.get("contract") != "chapter1_nee_channel_taxon_catalog_v1"
            or gc.get("source_version") != "0.9"
            or gc.get("source_sha256") != GLOBI_V09_SOURCE_SHA
            or gc.get("n_confirmatory_catalog_taxa") != 6928
            or gc.get("uses_focal_plant_traits") is not False):
        raise ValueError("Original pre-N1 GloBI original source gate drifted")
    return pd.read_csv(io.BytesIO(eb)), pd.read_csv(io.BytesIO(gb)), receipt, gc


def audit(weight: pd.DataFrame, catalog: pd.DataFrame, *,
          expect_source_counts: bool = True) -> dict:
    wreq = {"original_scientific_name", "diet_nectar_fraction"}
    creq = {"channel_id", "pollinator_species", "pollinator_family",
            "catalog_tier", "n_pollination_rows", "n_flower_visit_rows"}
    if not wreq <= set(weight) or not creq <= set(catalog):
        raise ValueError("Required original bird/candidate source fields unavailable")
    if (weight.original_scientific_name.isna().any()
            or weight.original_scientific_name.duplicated().any()):
        raise ValueError("EltonTraits bird species identity cannot be aligned uniquely")
    if (catalog.pollinator_species.isna().any()
            or catalog.duplicated(["channel_id","pollinator_species"]).any()
            or not catalog.catalog_tier.eq("confirmatory").all()):
        raise ValueError("Original GloBI independent confirmatory species catalog drifted")
    q = pd.to_numeric(weight.diet_nectar_fraction,errors="coerce")
    if not np.isfinite(q).all() or not q.between(0,1).all():
        raise ValueError("Nonfinite or out-of-bounds nectar fraction")
    if expect_source_counts and len(weight) != 9993:
        raise ValueError("Elton original accepted species denominator changed")
    w = weight[["original_scientific_name","diet_nectar_fraction"]].rename(
        columns={"original_scientific_name":"pollinator_species"})
    birds = catalog.loc[catalog.channel_id.eq("flower_visiting_birds")].copy()
    moths = catalog.loc[catalog.channel_id.eq("lepidoptera")].copy()
    if expect_source_counts and (len(birds) != 370 or len(moths) != 2360):
        raise ValueError("Original GloBI independent channel denominators changed")
    # Literal accepted original scientific names, no coarse genus mapping,
    # synonym imputation, or outcome-informed mapping. Values missing -> NA.
    birds["pollinator_species"] = birds.pollinator_species.astype(str).str.strip()
    w["pollinator_species"] = w.pollinator_species.astype(str).str.strip()
    joined = birds.merge(w,on="pollinator_species",how="left",
                         validate="one_to_one",indicator=True)
    if joined.pollinator_species.duplicated().any():
        raise ValueError("Unexpected ambiguous bird-name match")
    n_unmatched = int(joined._merge.eq("left_only").sum())
    matched = joined.loc[joined._merge.eq("both")].copy()
    if matched.empty:
        raise ValueError("No exact accepted species matches between original sources")
    positive = matched.diet_nectar_fraction.gt(0)
    zero = matched.diet_nectar_fraction.eq(0)
    if not (positive | zero).all():
        raise ValueError("Unresolved matched bird diet fraction")
    claims = pd.to_numeric(matched.n_pollination_rows,errors="coerce")
    visits = pd.to_numeric(matched.n_flower_visit_rows,errors="coerce")
    if claims.isna().any() or visits.isna().any() or claims.lt(0).any() or visits.lt(0).any():
        raise ValueError("Original GloBI claim/visit counts invalid")
    # Species-level 2x2 is a post-hoc *documentation* association; it is
    # phylogenetically and publication non-independent, not causal evidence.
    a = int((positive & claims.gt(0)).sum())
    b = int((positive & claims.eq(0)).sum())
    c = int((zero & claims.gt(0)).sum())
    d = int((zero & claims.eq(0)).sum())
    odds, p = fisher_exact([[a,b],[c,d]])
    moth_family = moths.pollinator_family.fillna("").astype(str).str.strip()
    sph = moths.loc[moth_family.eq("Sphingidae")].copy()
    claims_sph = pd.to_numeric(sph.n_pollination_rows,errors="coerce")
    visits_sph = pd.to_numeric(sph.n_flower_visit_rows,errors="coerce")
    if claims_sph.isna().any() or visits_sph.isna().any():
        raise ValueError("Sphingidae original claim/visit counts invalid")
    if expect_source_counts and len(sph) != 92:
        raise ValueError("Original GloBI Sphingidae records changed")
    total = len(birds)
    return dict(
        contract="fd_original_elton_globi_guild_catalog_bridge_v1",
        role="original_source_independent_species_diet_vs_guild_documentation_audit_not_island_test",
        original_eltontraits_diet_species=int(len(weight)),
        original_globi_confirmatory_bird_species=int(total),
        original_globi_confirmatory_lepidoptera_species=int(len(moths)),
        n_exact_name_matched_bird_species=int(len(matched)),
        n_bird_catalog_unmatched_diet_species=n_unmatched,
        n_matched_bird_diet_nectar_positive=int(positive.sum()),
        n_matched_bird_diet_nectar_zero=int(zero.sum()),
        n_pollination_claim_positive_diet=a,
        n_no_pollination_claim_positive_diet=b,
        n_pollination_claim_zero_diet=c,
        n_no_pollination_claim_zero_diet=d,
        matched_bird_positive_nectar_share=float(positive.mean()),
        catalog_bird_nectar_positive_lower_bound=float(positive.sum()/total),
        catalog_bird_nectar_positive_upper_bound_if_unmatched_all_positive=float(
            (positive.sum()+n_unmatched)/total),
        exploratory_species_fisher_odds=float(odds),
        exploratory_species_fisher_p=float(p),
        n_globi_confirmatory_sphingidae=int(len(sph)),
        n_sphingidae_with_pollination_claim=int(claims_sph.gt(0).sum()),
        n_sphingidae_with_flower_visitation_claim=int(visits_sph.gt(0).sum()),
        original_elton_raw_sha256=ELTON_SOURCE_SHA,
        original_globi_raw_sha256=GLOBI_V09_SOURCE_SHA,
        raw_species_island_plant_source_ids_exported=False,
        independent_spatial_guild_maps_available=False,
        no_eBird_data_used=True,
        biological_nondetection_or_shift_inferred=False,
    )


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--elton-zip",type=Path,required=True)
    p.add_argument("--globi-zip",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    w,g,_,_=_read_sources(args.elton_zip,args.globi_zip)
    result=audit(w,g)
    body=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    forbidden=('"pollinator_species"', '"accepted_species"', '"island_id"',
               '"first_source"', '"latitude"', '"longitude"', '"entity_ID"')
    if any(f in body for f in forbidden):
        raise ValueError("Source taxonomy or geography escaped count-only report")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n",encoding="utf-8")
    print(json.dumps({"n_matched":result["n_exact_name_matched_bird_species"],
                      "n_positive":result["n_matched_bird_diet_nectar_positive"],
                      "n_sphingid_claims":result["n_sphingidae_with_pollination_claim"],
                      "spatial_guild_test":False}))


if __name__=="__main__":
    main()
