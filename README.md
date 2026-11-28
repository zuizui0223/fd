# FD — external-source functional diversity and pollinator channels

**Source-first boundary (28 November 2026):** this repository owns only verified external biological inputs, independently of analysis technique. Project GBIF island incidence, integrated floral/reproductive trait ledger, derived aggregates or cohorts selected from those data belong to [island](https://github.com/zuizui0223/island), even when mixed with external source pools or animal observations. This supersedes the theme-first assignment in PR #76. Prior shared-source results and code are preserved for provenance only, not active FD research or independent evidence. See [island's source-first contract](https://github.com/zuizui0223/island/blob/docs/source-first-ownership-20261128/docs/source_first_research_ownership_20261128.md).

## Active research

- Publisher GIFT native lists with external NHM, Wang and AusTraits floral measurements; [coverage and effect audit](docs/EXTERNAL_FLORAL_REPLICATION_INTAKE_20261126.md). Raw length direction is promising but holdout geographic and conditional genus gates failed.
- Independently sourced bird, bee, butterfly and bat occurrence/range/checklist, diet and functional trait studies.
- Independent EMLN visitation networks, external nectar and reward datasets, and GloBI only without own-DB cohort dependence.
- No project botanical inputs, including public aggregates, in new active FD workflows. Shared geography and generic code are neutral if not fitted to own biological observations. Source-file independence is not proof of independent original observations.

## Archive and boundary checks

[Shared-source archive](shared_db_research/README.md) is read-only historical provenance; no deletion or new raw transfer. Own-dependent analyses should continue in island. [External source index](docs/EXTERNAL_SOURCE_RESEARCH_INDEX.md).

```sh
python scripts/check_external_source_boundary.py
python -m pip install -e '.[test]'
python -m pytest -q
```

Static file allowlist and own-source marker screens do not prove all transitive source ancestry. Audit source chains and cohort selection before activating new workflows. Preserve all nulls and avoid duplicate publication or claiming measured pollinator replacement from morphology.
