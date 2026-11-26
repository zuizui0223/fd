# External floral replication intake — 2026-11-26

The original FD results are in main (merge fb07764971ea54074e3ab84b7ecb50204165c0c4). The common-carrier descriptive color/form contrast survives both candidate source definitions; the isolation-gradient primary is unsupported (p=.726). This intake cannot rescue that negative result.

Publisher-only inputs:
- NaturalHistoryMuseum/plant-traits, pinned commit 4843c22bb2fc91018d9828dcddebf067c83b1097, CC BY 4.0, project data DOI 10.5519/p3dm31kc. Publisher repository directly distributes the public combined CSV and separate Flora of China, Flora of North America and Ecoflora tables. Verify exact Git blob hashes; previous NHM portal request in PR72 returned HTTP 403. This is a separately published public distribution, not a bypass of restricted data.
- AusTraits 5.0.0, DOI 10.5281/zenodo.10156222 (fixed version, not latest), original published ZIP, expected MD5 a4301cd000153533c0391abd9f60b889 and SHA256 recorded after retrieval. Inspect actual traits member, trait definitions, contributing source IDs and species-level consistency. Perianth fusion is not floral form; perianth symmetry is not automatically homologous to every floral symmetry definition.
- Wang et al. (2022/2023), DOI 10.5061/dryad.ghx3ffbrh, global symmetry dataset is a candidate next source. Audit direct versus assigned genus/family states and taxonomy before use. No download yet.
- Official GIFT publisher native checklists are a potential independent roster (PR72 source intake); previously limited published flower colour coverage. Need public metadata-only cohort selection, nonoverlap, native-status validation and name/ID reconciliation. No checklist rows or project data used in this run.

This source-only readiness step does not estimate an ecological effect. Exact-binomial intersections are conservative, not accepted-name harmonization. String consistency across source values is not validation of measurement or primary floral colour. Source texts may overlap; independent file provenance does not establish independent underlying observations. A later external test must lock comparable endpoint definitions, source references, region and eligible cohorts before fitting, and report missingness. Do not equate corolla length with depth/accessibility or inferred syndrome with realized pollinators. No raw data redistribution in this branch.

Implementation: scripts/external_floral_replication_intake.py. Outcome pending.
