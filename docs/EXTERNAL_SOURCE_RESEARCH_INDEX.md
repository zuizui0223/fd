# External-only research source index

## Retained source families

| Family | Source basis | Active/archive assets |
|---|---|---|
| Bird/checklist/diet/flight | independent DARs, Baldwin/Myers, EltonTraits, AVONET | `scripts/global_bird_source_audit.py`, `scripts/tropical_bird_hwi_occupancy.py` |
| Native bees and butterflies | Marshall source, Yau occurrence archive, LepTraits | `scripts/tropical_bee_bird_crosswalk.py`, `scripts/tropical_butterfly_trait_filtering.py`, associated source reports |
| Independent empirical networks | EMLN60 and EMLN71, original pinned source | `src/fd/emln_island_visitor_turnover.py`, `scripts/run_emln_island_visitor_turnover.py`, tests; source workflow retained under docs |
| Bird nectar source | EltonTraits/AVONICHE and GloBI external taxon catalog without own plant inputs | `src/fd/eltontraits_bird_nectar.py`, `src/fd/avoniche_nectar_preflight.py`, `src/fd/elton_globi_guild_bridge.py`, original aggregate |
| Bat metadata and independent ranges | MDD and original EltonTraits | `scripts/audit_mdd_bat_diet_metadata.py`, `scripts/bat_mdd_synonym_audit.py`, `scripts/export_mdd_bat_ranges.R` |

External animal records that were used as predictors of our own plant-flora aggregates are **mixed** and move with that analysis to island. External source acquisition or stand-alone animal analyses are retained here. Historical source artifact IDs and immutable original results remain traceable; no raw restricted records are copied. Original mixed source documents have not been rewritten as external-only evidence.

The external-only path list is intentionally conservative. Historical variants not source-audited are held in the source-pinned archive; do not treat them as proven external-independent. Current own/mixed archive checkpoint: [island ownership contract](https://github.com/zuizui0223/island/blob/64a58cd4817ccf1c54c179aa7fd1c5585f57dd43/docs/fd_data_ownership_20261125.md). No deletion of Git history.
