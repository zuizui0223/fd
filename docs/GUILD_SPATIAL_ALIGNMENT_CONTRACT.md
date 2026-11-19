# Frozen test: tropical radial × deep-tube flowers versus mapped pollinator guilds

**Status 2026-11-19:** implemented with synthetic tests; **not executed on original plant×island or pollinator-map data**. Existing GIFT 4,505-island positive results are unchanged. This is an exploratory ecological association, not an observation of visits, pollen transfer, fitness or a historical pollinator shift.

## Primary estimand and single joint test

The existing original tropical actinomorphic × deep-tube, source-relative Direct motif excess is approximately +10.99 percentage points (666 eligible floras), driven largely by between-genus botanical composition. The new, independent test asks whether island-minus-eligible-mainland candidate SOURCE differences in this actual same-species pair motif co-vary with mapped, source-relative nectarivorous-bird and Sphingidae guild availability.

Freeze: all 4,505 original GIFT floras as parent; tropical analysis region; all-recorded botanical status; first five ranked GIFT native mainland candidate source regions ("geo_k5"); strict **same accepted species** Direct actinomorphic and deep-tube sector; require at least five measured same-species pairs both for island and source, plus 60 eligible islands, 15 disjoint 10-degree physical blocks, and 15 distinct first-source groups. No cross-products of separate species-level single-trait margins.

Response: Direct known-pair sector proportion on island minus the corresponding source pool proportion. Unmeasured species are **missing**, not sector-negative. The primary two independent covariates are log(1+mapped qualified nectarivorous-bird range species richness) and log(1+qualified Sphingidae modelled species richness), island-minus-source and standardized within the eligible tropical test cohort. Guild variables come only from independently acquired animal data, without focal flower traits. Neither zero GBIF detections nor generic Aves/Lepidoptera record totals count as qualified absence or moth/nectarivory measures.

Fit one joint 2-df null of no guild association, with both guild terms entered together, controlling log(1+distance), log(1+area), absolute latitude, temperature and precipitation PCs, and separately recorded original Direct island/source pair coverage. Original disjoint 10-degree physical blocks and first GIFT source group each receive equal cluster weight with separate CR1 F approximations and 499 group-resampled intervals. A directional positive *exploratory gate* requires both partial guild slopes nonnegative, joint p<0.05 by both group definitions, and the lower CI for at least one guild greater than zero under both groupings. The joint p-values are exploratory, not family-wise confirmatory nor tests of causality.

## Source-independent map requirements

**Birds:** identify nectar resource dependence from AVONICHE quantitative diet/foraging niche (https://doi.org/10.5061/dryad.1zcrjdg56) and check native resident/breeding verified range polygons on exact island geometry, not centroid radii. BirdLife polygons require noncommercial research data access/permission; do not publish them or species-specific locality/identity joins.

**Sphingidae:** require species-resolved moth range/SDM surfaces with spatial and temporal GBIF detection-effort checks. Plain GBIF occurrence absence is not sufficient to assign zero richness; unknown must remain NA. Distinguish adult nectarivory from merely being in Sphingidae.

**Study population:** genuine original GIFT native source floras and original GSHHG polygon join, with independent 10-degree blocks. Map scores must be built blind to focal floral response. No eBird-based data stream is required.

## Recording coverage and limitations

The reported tropical symmetry coverage difference (island ~47%, source/mainland ~16%) is substantial, but is NOT automatically the same as Direct symmetry×depth **joint** coverage in the analysis. Use the genuine observed pair count and full accepted-species denominator separately for each island and source. The runner adjusts for measured coverage but explicitly returns worst-case bounds that permit every missing pair to be sector-positive or sector-negative. If zero lies in these bounds, a positive guild coefficient cannot establish robust direction of the total sector difference under arbitrary MNAR. A separate source-blind genus-conditioned sensitivity, phylogenetic and endemism checks are needed before interpreting guild-associated botanical sorting as lineage-level adaptation.

## Executable status and privacy

The ready-to-run module accepts a **private** one-row-per-original-island CSV with exactly 4,505 rows and the columns enforced by the code. The private join must be prepared inside authorized source CI using original plant and separate independently reviewed guild maps. Unknown and non-qualified guild states must not be converted to zero. Execute:

    python -m fd.guild_sector_alignment --private-island-join /secure/original_4505_guild_pair_join.csv --bootstrap 499 --output /secure/anonymous_guild_pair.json

Only anonymous aggregate JSON, an input source manifest and code/test receipts may be published after real execution. A non-estimable support status is a valid result. No protected island, plant or animal identities should be committed to the public fd repository. Do not merge this mechanism into Main until original map and GIFT source hashes, model support, uncertainty and missingness bounds have all been checked.
