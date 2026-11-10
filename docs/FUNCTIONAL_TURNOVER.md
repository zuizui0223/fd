# Functional-trait turnover under island taxonomic turnover

**Status: exploratory, descriptive-only input/coverage gate (2026-11-10). No new empirical functional-diversity, null-model, or evolutionary finding is claimed.**

## Scientific question and independence from island Chapter 1

The existing fd category-profile analyses measure **recorded** colour/form configurations, not actual pollination performance. The upstream island Chapter 1 already tests reproductive-assurance, floral architecture, and pollen-limitation associations with isolation. This independent fd extension asks:

**When different islands replace constituent plant taxa, is reproductive-assurance or floral-access structure retained more than lineage membership, geographic source availability, and observation coverage predict?**

That question is distinct from reprinting the Chapter 1 isolation slopes. Functional similarity among different species is a hypothesis, not a demonstrated biological equivalence.

## Input dimensions

Use seven source-verified, binary plant traits. Analyse each separately; do not invent a universal mixture of colour, attraction and selfing:

| Domain | Input columns | Interpretation |
| --- | --- | --- |
| Reproductive assurance | self_compatibility, selfing_mating_system, autonomous_selfing | Capacity to reduce dependence on mates/pollen vectors; **not** directly measured seed set |
| Floral accessibility | generalized_form, actinomorphic_symmetry, shallow_open_tube | Visitor access proxies; **not** realized visitor identity |
| Signal proxy | plain_colour | Visual trait proxy; **not** measured pollinator response |

Independent GloPL pollen-addition experiments may be used for **external, species-matched functional interpretation** but are not contemporaneous fd island-pair reproductive outcomes.

The upstream corrected island database includes 106,295 accepted species and 222,688/318,885 resolved species-axis cells. The raw reproductive-assurance axis resolves only 48,497/106,295 species. These are upstream denominators, not proof of usable complete support in fd near-island matched comparisons. Avoid global complete-case requirements that leave tiny, unrepresentative support.

## Rights-protected staging contract

A local, approved private CSV must have **one row per island_id x taxon_name** in an identical, explicitly documented All *or* Direct evidence scope. Required fields are:

- island_id, taxon_name (accepted species/name identifier), genus;
- seven listed binary traits (0/1/empty). Empty means **unknown**, not absent or zero.

The paired-island CSV has island_a, island_b, comparison, with donor categories defined *without* looking at floral outcomes. Never combine All/Direct as independent replications. A shared taxon must not have inconsistent trait values or genus between islands within one staging scope.

A runnable diagnostic that preserves these rules is in src/fd/functional_turnover.py. It measures observed trait prevalence gaps for **all taxa**, **each island's species-exclusive taxa**, and **each island's genus-exclusive taxa**. The last two explicitly remove shared taxonomic identities. Every species-level datum stays in local memory; by default each trait requires >=20 observed name carriers per island. Unresolved pairs remain missing, never scored zero.

Example:

    python -m fd.functional_turnover \
      --private-traits /protected/accepted_all_traits.csv \
      --private-pairs /protected/outcome_blind_island_pairs.csv \
      --output work/private_fd_functional_audit \
      --min-carriers 20 --minimum-pairs 10

Output is only a grouped trait-scope x source-comparison table and status flag in an **ignored local work/ directory**. The code never exports island names, species names, genus names, or raw pair-level data. Do not add private input to any public workflow or release; conduct rights review before releasing any further aggregates.

## Three gates to a genuinely new macroecological finding

1. **Real coverage and attrition first.** Check how many physically distinct islands, named archipelagos, spatial blocks, unique species and genera actually survive for each trait *after* shared-taxon/genus removal. Confirm identical outcome-blind pairs and All/Direct scope sensitivity; do not assume the earlier 606 colour/form islands all have seven traits.
2. **Geography/source negative controls.** Compare local islands with geographically distance-matched other-archipelago islands *at equal taxonomic and trait-observation opportunity*, controlling area/climate and regional native source pools. Absolute trait-prevalence differences are descriptive; a negative gap is not by itself evidence of stabilizing selection or functional retention.
3. **Conditional inference, not volume chasing.** Hold observed island x taxon membership and shared-taxon structure fixed. Reassign *whole observed trait vectors* within defensible genus x source/region opportunity strata, retaining missingness, and use spatial components/source archipelagos as inferential units. A genuine functional-retention result must survive species identity removal, source/lineage control, and geographic independence. Document support loss and preregister tests before examining outcomes.

A color/form PCoA area or 17D KDE Gaussian volume is **not** a substitute for reproductive function. Pollinator guilds inferred from colour/shape and GloBI interaction occurrence are not measured pollination effectiveness. Chapter 1's already reported H3/H4 results cannot count as independent fd confirmations.

## Provenance and inference ceiling

- island Chapter 1 current: submission/chapter1_current/H1_H4_EVIDENCE_MAP.md, docs/DATABASE_BUILD.md; self-compatibility recurrent, corrected experimental pollen-limitation association, post-hoc function-trait consistency.
- fd's local floral profile result: results/20261109_archipelago_flower_transfer/README.md. Its 606-island local prediction does **not** exclude the same plant taxa.
- The continuing source/genus and joint-colour/form prediction work is in fd PR #11; its small coupling residual is a **phenotype prediction**, not measured functional replacement.
- The bundled synthetic tests validate correct missingness handling, species/genus exclusion, shared taxon consistency and only aggregated output. **They do not establish any real-data significance.**

**Stop rule:** No language about ecological niche replacement, pollination equivalence, adaptive syndrome evolution or universal function conservation without an independently verified biological functional endpoint and spatially/source-robust real-data result.
