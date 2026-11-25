# FD — floral functional diversity and pollinator channels

**Scientific ownership (effective upon merging PR #76):** FD owns floral functional diversity, syndrome composition and pollinator-channel research regardless of source. island retains canonical data custody, protected execution and Chapter 1 reproductive-assurance research. See [shared-source research and ownership](shared_db_research/README.md). Shared-source snapshots are archival, not active external workflows; no raw data copied. Historical external-only safeguards remain active for external workflows.

## Research here
- Independent bird occurrence/checklist, diet, and flight morphology; bee/butterfly source geometry and trait audits; documented source limitations.
- Independent EMLN observed flower-visitor turnover, source-only bird nectar metadata and independent MDD bat source/trait audit. GloBI external catalog can be used if independent of focal plant records.
- Shared-database floral functional diversity, functional reassembly, pollinator-channel opportunity, and syndrome composition. See [FD evidence spine](shared_db_research/FD_EVIDENCE_SPINE.md). No biological claim of pollinator replacement follows from the external source audits.

## Boundary checks
```sh
python scripts/check_external_source_boundary.py
python -m pip install -e '.[test]'
python -m pytest -q
```
Static path allowlist and known-source marker checks are a safeguard, not complete proof of all transitive data provenance. A new or changed research source requires explicit dependency review and manifest update. No original data are newly published by this reorganization.

Historical fd reports may refer to the botanical analyses as context; ownership follows the scientific estimand in the shared-source registry. Preserve original results and uncertainty and cite the owning repository.
