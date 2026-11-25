# FD — independent external-database research

**Ownership decision, 25 November 2026:** any analysis using our own `island` GBIF-derived flora, integrated species trait ledger or derived aggregates now belongs to [island PR #347](https://github.com/zuizui0223/island/pull/347). This includes every mixed analysis containing those inputs, irrespective of whether an old report named its parent "GIFT4505". Official GIFT-only analyses using independently sourced trait data remain eligible here, but none of the historical mixed 4505 analyses is relabelled independent.

This branch is the external-only active worktree proposal. See [source index](docs/EXTERNAL_SOURCE_RESEARCH_INDEX.md) and [file ownership manifest](config/external_db_ownership.json). Legacy external source workflows from unmerged branches are preserved under `docs/legacy_external_workflows`, not autoactivated. The previous main and all old PRs remain retrievable from immutable Git history; their own/mixed research is to be archived in island, not cited as a second independent result here. Do not merge legacy mixed PRs into this worktree.

## Research here
- Independent bird occurrence/checklist, diet, and flight morphology; bee/butterfly source geometry and trait audits; documented source limitations.
- Independent EMLN observed flower-visitor turnover, source-only bird nectar metadata and independent MDD bat source/trait audit. GloBI external catalog can be used if independent of focal plant records.
- Future external-only functional-diversity research with explicitly named official source releases and provenance. No biological claim of pollinator replacement follows from the external source audits.

## Boundary checks
```sh
python scripts/check_external_source_boundary.py
python -m pip install -e '.[test]'
python -m pytest -q
```
Static path allowlist and known-source marker checks are a safeguard, not complete proof of all transitive data provenance. A new or changed research source requires explicit dependency review and manifest update. No original data are newly published by this reorganization.

Historical fd reports may refer to the botanical analyses as context; those figures/claims belong to island. Preserve original results and uncertainty and cite the owning repository.
