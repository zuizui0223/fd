# Same-species pollinator-shift feasibility: original GloBI v0.9 visitor evidence

**12 November 2026 · fd PR #14 · original-source, evidence-gated extension; ecological guild switching NOT established.**

## Why this is stricter than the original 928-island comparison

A geographical overlap between island plant species and a candidate GIFT native mainland list establishes a *source-opportunity NAME MATCH*, not that a documented island plant was visited by an independently identified pollinator or that a mainland population was ancestral. Similarly, a globally documented GloBI visitor guild associated with the **plant GENUS** is not direct evidence that an identified plant species or island population received pollen from that visitor.

We therefore rebuilt the **exact source-matched species** candidate set using ONLY:
- original Direct high/medium evidence with complete colour, form and symmetry;
- original source-verified native island flora, retaining only species present in one of the five geographically closest candidate GIFT native mainland source checklists;
- island record cohort with ≥30 directly trait-complete species;
- **qualified GBIF Bombus nondetection** and independent recorded flower-visiting bird animal occurrence.

This deterministic original frozen source selection reproduces **928 original eligible islands**, **21** Bombus-adequately-nondetected/bird-occurred islands, **197** native source-matched island×species rows, **54 accepted species / 49 genera**, on **13 focal islands**. Individual species, original island/source IDs and original GIFT checklist identities remain in the protected study workspace.

## Baseline direct floral / reproductive evidence among all 54 candidates

| Trait | Original high/medium Direct known species | Missing | Comment |
| --- | ---: | ---: | --- |
| Flower colour | **54** | 0 | Species-wide categorical record only, not island population photo phenotype |
| Floral form | **54** | 0 | Open radial in 9/54, not universal across candidate pool |
| Floral symmetry | **54** | 0 | Actinomorphic in 38/54 |
| Flower size class | 30 | 24 | Direct categorical size, not flower-specific tube/fit measure |
| Tube depth class | 22 | 32 | Not all source-island pairs can compare potential access morphology |
| Self-compatibility / incompatibility | **17** | 37 | Only 17/54 direct SC/SI records; **not selfing rate or reproductive assurance** |
| Autonomous selfing capacity | 12 | 42 | One original state is mixed/variable, thus 11 unambiguous single categories |
| Mating system | **8** | 46 | Too sparse to estimate macro-scale selfing-shift coevolution |
| Inflorescence display | 42 | 12 | One direct state is mixed, not force-coded |
| Cleistogamy | 7 | 47 | Incomplete |

Using the previously frozen **GloBI global plant-genus** reference breadth data only, 33/54 candidates have a genus with ≥3 independent documented visitor-interaction references, and 15/54 have a genus with some bird documentary records. **Only 3 of those 15** have original open-radial floral form. These are taxonomic/documentation overlaps, not observed visitors of the 54 individual plant species.

### Cross-evidence attrition inside the native 54-species cohort

The old *globally documented GENUS-level* GloBI visit/breadth ledger yields ≥3 reference support for **33/54** species' genera and any global bird-channel reference for **15/54** species' genera. Among those **15** genus-bird-documentation candidates, only **6** have Direct species-level SC/SI evidence and **5** have Direct floral tube-depth evidence. Of the full 54, **9** have original Direct open-radial form (just **3** of the 15 genus-bird documentation candidates). Strict Direct high/medium reproductive record coverage is therefore inadequate for a bird-versus-bee pollination-service or mating-system adaptation model. These statements do **not** indicate actual bird interactions with any of the 54 named species.

The three older morphology-template candidates are poorly discriminating in exactly this 54-species cohort: only **5** have ≥3 Direct informative traits for `bird_like` and **4** have `bird_like` concordance ≥0.75 on those ≥3 traits. Only **3** candidates have ≥3 Direct informative traits for EACH of bee-like, butterfly-like and bird-like, and all three have high score in **all three templates**. These are overlapping morphology templates, not verified guild membership. Anonymous cohort aggregate: [`native_54_morphology_template_ambiguity.csv`](native_54_morphology_template_ambiguity.csv) and [`native_54_species_direct_functional_trait_coverage.csv`](native_54_species_direct_functional_trait_coverage.csv).

## Exact plant species × visitor reference retrieval implemented

To avoid projecting genus labels onto species, we implemented original GloBI v0.9 source-stream parsing in [island Draft PR #323](https://github.com/zuizui0223/island/pull/323): `src/island_v2/fd_globi_species_shift_v09.py`, `tests/test_fd_globi_species_shift_v09.py`, and the existing v0.9 SHA-pinned workflow.

The scanner:
1. Independently reconstructs and freezes the **same 54 native source-matched accepted plant names** from the original progressive island source and archived GBIF observation-state files; fails closed on sample drift.
2. Reads the exact **GloBI 2026-06-04 v0.9 `interactions.tsv.gz`** (SHA256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`) in the SAME source pass used by the established provider audit.
3. Matches the **full plant species name**, not merely genus, on the correctly oriented source or target side of each of the four frozen flower-visitation/pollination RO predicates. Rejects records lacking a real plant kingdom, an identified visitor species, or one of the five preselected visitor guilds (Bombus, other bees, undivided Lepidoptera, birds, Diptera).
4. Separately counts raw exact-species hits, available original **provider AND bibliographic reference** evidence, duplicate references shared across provider labels, distinct cited interactions by species and guild, records explicitly labeled `pollination` vs `flower_visit`, plus **availability of record coordinates/event times**. A coordinate is NOT treated as exact island polygon membership.
5. Returns ONLY anonymous counts: number of the 54 species with any source-cited interaction, ≥3 distinct cited references, ≥2 recorded guilds, ≥2 source-archive labels, any georeference/date, and each of five guilds' count. **Zero documented matching records means GloBI evidence absent, not true pollinator absence**.

Source/role/missing-documentation/reference reuse and candidate selection synthetic tests all passed locally; the GitHub workflow performs source-pinned original-data tests before publishing a restricted-ID-free receipt. At the time of this report, **the original v0.9 stream was still running; no empirical exact-species GloBI support count is claimed yet**. A positive global species-linked record, when found, is STILL NOT proof it belongs to the focal island or the candidate GIFT mainland source population.

## Biological question and decision conditions

The paper should treat three nonidentical hypotheses separately:
- **Flora sorting / source-pool filtering:** existing floral trait states and taxon frequencies differ between the island flora and geographically plausible mainland plant pools. Already observable, but no ecological mechanism established.
- **Real ecological pollinator switching:** the *same native plant species* or a supported island sister lineage has standardized and *location-specific* pollen-effective visits from different guilds in the source mainland and on the focal island.
- **Functional compensation / generalized visitation:** source pollinator reduction does not induce a one-for-one switch; multiple visitors could share pollination service on the island, possibly with no large change in floral structure.

The exact 54-species GloBI scan is **evidence feasibility**, not the endpoint: an independently verified island and mainland coordinate assignment or direct population interaction dataset is required to claim switching. Additional original flower phenotype **by population**, effective pollen deposition/fruit set or controlled pollen-limitation tests are required to make a claim about syndrome evolution or service compensation.

Only the first hypothesis currently has source-pinned descriptive results; the other two are not rejected, simply **unidentified in existing input schemas**. We keep fd PR #14 Draft and do not count globally genus-linked bird references as 15 actual bird-pollinated island species.

## References and provenance
- Frozen `island` progressive source archive `10058653212`, SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.
- Exact-island GBIF qualified channel evidence artifact `10314984302`.
- GloBI source genus v2 artifact `10328696801`, SHA256 `475dda8b0eef4ce22ad1aa736fa82cc049042d801217c65d52641489288ca254`.
- Original GloBI v0.9 pinned DOI `10.5281/zenodo.20546682`, SHA256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`.
- Original GloBI documented interaction records include species identity, source/provider and bibliographic citation keys and optional `decimalLatitude`, `decimalLongitude`, `eventDateUnixEpoch`. Cite original GloBI contributing datasets and do not redistribute protected raw input; source-specific anonymity thresholds apply to future releases.
