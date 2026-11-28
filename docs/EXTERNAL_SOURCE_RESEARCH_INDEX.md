# External-source research index

## Retained source families

| Family | Source basis | Active/archive assets |
|---|---|---|
| Bird/checklist/diet/flight | independent DARs, Baldwin/Myers, EltonTraits, AVONET | `scripts/global_bird_source_audit.py`, `scripts/tropical_bird_hwi_occupancy.py` |
| Native bees and butterflies | Marshall source, Yau occurrence archive, LepTraits | `scripts/tropical_bee_bird_crosswalk.py`, `scripts/tropical_butterfly_trait_filtering.py`, associated source reports |
| Independent empirical networks | EMLN60 and EMLN71, original pinned source | `src/fd/emln_island_visitor_turnover.py`, `scripts/run_emln_island_visitor_turnover.py`, tests; source workflow retained under docs |
| Bird nectar source | EltonTraits/AVONICHE and GloBI external taxon catalog without own plant inputs | `src/fd/eltontraits_bird_nectar.py`, `src/fd/avoniche_nectar_preflight.py`, `src/fd/elton_globi_guild_bridge.py`, original aggregate |
| Bat metadata and independent ranges | MDD and original EltonTraits | `scripts/audit_mdd_bat_diet_metadata.py`, `scripts/bat_mdd_synonym_audit.py`, `scripts/export_mdd_bat_ranges.R` |

External animal records that are used as predictors of our own plant-flora aggregates are **mixed**. Under the 28 November source-first policy, own-dependent functional-diversity and channel-response estimands are island-owned; FD retains independently sourced external-only research. External source acquisition or stand-alone animal analyses are retained here. Historical source artifact IDs and immutable original results remain traceable; no raw restricted records are copied. Original mixed source documents have not been rewritten as external-only evidence.

The external-only path list is intentionally conservative. Historical variants not source-audited are held in the source-pinned archive; do not treat them as proven external-independent. Current own/mixed archive checkpoint: [source-first island policy (PR #361)](https://github.com/zuizui0223/island/pull/361); the earlier contract is historical. No deletion of Git history.


## Source-first result ownership and official GIFT distinction (2026-11-28)

The controlling unit for **both source classification and scientific ownership**, is the actual biological input ancestry of a result, including selection of the analysis cohort, training data, labels, and derived aggregates. Neither island count nor the historical label "GIFT4505" identifies a data source.

- A project-built GBIF island-species roster or project-integrated floral/reproductive ledger is own data; any result dependent on either must be labelled shared-source and use protected execution in island. Scientific ownership belongs to island regardless of the estimand. Official GIFT mainland donor pools, animal data, and external experiments do not change that owner.
- A publisher GIFT island checklist and independently retrieved publisher trait values, with no own biological data used even for cohort selection, are external-only and are eligible for FD's functional-diversity research. Shared generic code and outcome-blind geography/climate alone are neutral.
- External-only is a statement about the provenance of the inputs, not proof of independent underlying species observations or source publications. Unknown dependencies must be audited before promotion.
- The previously archived EuPPollNet guild-removal simulation selected eleven islands from the own 928-island morphology-complete cohort, hence its selected result is shared-source, and therefore belongs to island despite the external network endpoint. A fresh external-only network cohort would be eligible here.
- Standalone GloPL experimental isolation evidence is external-source context in the existing frozen island H3, not a new independent fd result. The GloPL × own-trait H4 and FD bridge are mixed.
- [Official GIFT source pilot (fd #72)](https://github.com/zuizui0223/fd/pull/72) is coverage/readiness, not a fitted ecological effect. Source audit run 45560668861 retrieved official publisher values: strict colour × maximum flower length shared 404 work IDs globally, but in the bounded 24-island-list pilot only two lists had any joint coverage (maximum four carriers). Do not treat publisher self-fertilization as automatically equivalent to self-compatibility or syndrome categories as measured pollinators.

The [historical source-classification registry and validator in island PR #353](https://github.com/zuizui0223/island/pull/353) retain source references, support and claim limits. Its inventory is not an independent external replication and does not alter the frozen GEB submission.

## New external-only network research (2026-11-28)
[Same-plant conditional bee loss](EXTERNAL_EUPPOLLNET_BEE_LOSS_20261128.md): direct EuPPollNet publisher records and public geography only; corrected island–mainland primary nonestimable (zero same study/year/plant pairs).

- [Channel retention](EXTERNAL_EUPPOLLNET_CHANNEL_RETENTION_20261128.md): +.039528, initial frozen gate passed; two-exclusive sensitivity unresolved.
- [Floral moderation](EXTERNAL_EUPPOLLNET_ARCHITECTURE_20261128.md): both frozen morphology gates failed.
- [Fixed-margin conditional null](EXTERNAL_EUPPOLLNET_CURVEBALL_20261128.md): +.044582, original degree and overlap controlled, supported; two-exclusive sensitivity unresolved.
- [Geographic scaling and crossed dependence](EXTERNAL_EUPPOLLNET_SPATIAL_SCALING_20261128.md): crossed study/genus uncertainty and positive spatial slope unsupported.

- [Independent EMLN source readiness](EXTERNAL_EMLN_SPATIAL_AUDIT_20261128.md): 42 source directories, 16 audits, only one strict candidate study; no independent macro replication.

- [Measured EMLN morphology](EXTERNAL_EMLN_TRAIT_TURNOVER_20261128.md): primary single-study morphology case opposite direction; no global inference.
- [EuPollTrait bridge](EXTERNAL_EUPOLLTRAIT_INTAKE_20261128.md): new external trait intake for same-network common-carrier comparison.

- [EuPollTrait group-specific size turnover](EXTERNAL_EUPOLLTRAIT_SIZE_TURNOVER_20261128.md): intake verified; fixed-margin size test with crossed dependence, pending results.
