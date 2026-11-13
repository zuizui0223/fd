"""Real original EMLN60 and EMLN71 site-matched flower-visitor turnover.

Only original source CSVs exported *losslessly* from pinned EMLN .rda are
allowed. Outputs anonymized location/period summary (not named plant taxa,
source interaction links, or exact coordinates). No historical loss claim.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from fd.emln_island_visitor_turnover import (
    original_emln_source,canary_same_plant_partner_turnover,
    mauritius_native_partner_rewiring,
)

EXPECTED_ORIGINAL = {
   "60":"e4f046f2cce1f84b4513e79611c731fe7cda37c94be1a0dab1726b901a30ddc2",
   "71":"89f5acf9a746975c42583a3cd329d5aa604a9b5e5f49534bdeb5168e4fa13e2b",
}
EXPECTED_BLOB = {
  "60":"a1bde277e2ac9cc269ea85f6ba9d55dd4a70b19f",
  "71":"4f1cf162871341b4b37116d83fc3e2bc728df58b",
}


def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--original-source-dir",type=Path,required=True)
    a.add_argument("--source-csv-dir",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    z=a.parse_args()
    data={}
    for network,nlayer in (("60",14),("71",24)):
        original=next(z.original_source_dir.glob(f"emln{network}_*.rda"),None)
        if original is None or sha256(original)!=EXPECTED_ORIGINAL[network]:
            raise ValueError("Original EMLN source SHA mismatch, stop")
        files={
            name:z.source_csv_dir/f"emln{network}_{name}.csv"
            for name in ("nodes","interactions","layers")
        }
        for file in files.values():
            if not file.is_file() or file.stat().st_size==0:
                raise ValueError("Missing original pinned EMLN original source table")
        n,e,l=original_emln_source(
          pd.read_csv(files["nodes"],dtype=str),
          pd.read_csv(files["interactions"],dtype=str),
          pd.read_csv(files["layers"],dtype=str),
          expected=nlayer)
        data[network]=(n,e,l)
    canary=canary_same_plant_partner_turnover(*data["60"])
    mauritius=mauritius_native_partner_rewiring(*data["71"])
    assert canary["n_geographical_islands"]==5
    assert canary["n_island_study_locations"]==6
    assert canary["n_independent_mainland_study_locations"]==1
    assert canary["n_total_replicate_source_sites"]==14
    assert canary["n_original_visitor_links"]==969
    assert canary["n_island_locations_where_strict_site_matched_within_exceeds_cross"]<=6
    assert mauritius["n_two_week_periods_with_two_survey_sites"]==12
    assert mauritius["n_independent_managed_sites"]==1
    assert mauritius["n_independent_unmanaged_invaded_sites"]==1
    assert mauritius["n_original_visitor_links"]==1953
    output={
      "status":"SOURCE_PINNED_OBSERVED_PARTNER_TURNOVER_AND_SITE_PERTURBATION",
      "reference_source_repo":"Ecological-Complexity-Lab/emln",
      "reference_source_commit":"e66e8702910728358e71f7814d01afd3c7b00c30",
      "source_original_sha256":EXPECTED_ORIGINAL,
      "source_original_git_blob":EXPECTED_BLOB,
      "canary_actual_plant_visitor_partner_turnover":canary,
      "mauritius_native_plant_rewiring":mauritius,
      "does_not_match_fd_source_original_botanical_trait_carrier_ids":True,
      "source_visitor_taxa_are_observational_not_pollen_effectiveness":True,
      "no_historical_pollinator_population_decline_or_functional_rescue_proven":True,
      "source_not_representative_global_island_sampling":True,
    }
    raw=json.dumps(output,sort_keys=True)
    for forbidden in ('"node_name"','"visitor"','"plant"','"node_id"',
       '"layer_from"','"node_from"','"Latitude"','"Longitude"'):
        if forbidden in raw:
            raise ValueError("Identifying original source taxon or exact site leaked")
    z.output.parent.mkdir(parents=True,exist_ok=True)
    z.output.write_text(json.dumps(output,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print(json.dumps({
       "status":output["status"],
       "n_canary_source_locations":canary["n_island_study_locations"],
       "n_canary_physical_islands":canary["n_geographical_islands"],
       "n_strict_canary_matching_within_greater_than_cross":
           canary["n_island_locations_where_strict_site_matched_within_exceeds_cross"],
       "n_mauritius_periods":mauritius["n_two_week_periods_with_two_survey_sites"],
       "mauritius_mean_native_plant_difference":mauritius[
          "period_equal_mean_original_metrics"][
          "n_native_visited_plant_difference_CMA_minus_Control"]
    },sort_keys=True))
if __name__=="__main__":
    main()
