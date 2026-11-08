# Native-mainland genus opportunity and 17-dimensional floral hypervolume

**2026-11-08 JST — independently executable `fd` extension.** [Methods and rights contract](../../docs/SOURCE_OPPORTUNITY_HYPERVOLUME17.md). Geographic/flower-profile sources pinned to [island PR #321](https://github.com/zuizui0223/island/pull/321). The results are **recorded flora / potential-mainland-source associations**, not actual plant colonization, confirmed absence, or adaptive pollinator response.

## 1. Source-only prediction of which genera are recorded on geological islands

GIFT native mainland source lists contain **748 regions and 1,666,186 native source-region×taxon records**. Each island has five independently assigned nearest mainland regions. We ranked the union of source-eligible genera without fitting any island occurrence outcome, using (i) genus prevalence in all mainland regions, (ii) number of local regions listing the genus, (iii) inverse-distance weighted local presence, or (iv) inverse-distance weighted **native genus name density normalized for source checklist size**.

**No-record is NOT a real absence**. ROC AUC and AP evaluate how these source-only ranks distinguish recorded from unrecorded genera **within the source-defined candidate universe**; they are not occupancy/establishment probabilities. Any genus missing from all five source lists is excluded from that ranking and included only in source coverage diagnostics.

The global geological source intersection includes **434 independently geological-typed islands**. Of those, **406** have adequate original all-recorded genus outcome support (174 continental, 232 oceanic) and **306** have enough explicitly native records (141 continental, 165 oceanic).

| Original island outcome | Geological island type | Number | AUC: global mainland genus prevalence | AUC: distance-weighted local genus density | Mean AP global → local |
|---|---|---:|---:|---:|---:|
| All recorded | Continental | **174** | 0.6760 | **0.7624** | 0.4143 → **0.5072** |
| All recorded | Oceanic | **232** | 0.7238 | 0.6989 | 0.2850 → 0.2922 |
| Explicit native | Continental | **141** | 0.6546 | **0.7769** | 0.3172 → **0.4517** |
| Explicit native | Oceanic | **165** | 0.7262 | 0.7266 | 0.1933 → 0.2385 |

Local source-density improvement in all-recorded continental islands is ΔAUC=**+0.0864**, island-bootstrap 95% [+0.0728,+0.1003]; in oceanic islands ΔAUC=**−0.0249**, 95% [−0.0410,−0.0080]. However this unpaired full-sample contrast mixes island-region/coverage differences and **does not** establish that geological origin changes arrival success.

**Matched geography interaction check:** Among the fixed [37 origin-matched pairs](../20261108_matched_flower_repertoire/README.md), the oceanic-minus-continental **local_density AUC improvement** is **−0.0274** (pair bootstrap 95% [−0.0503, −0.0048], exploratory sign-flip p=0.0246). Under **12 total paired comparisons** (three local scores × AP/AUC × two floristic-status definitions), **none survives BH/BY** (minimum BH q=0.2004). The explicitly native sensitivity has 31 status-eligible pairs and also no adjusted origin-specific gain. A smaller nominal difference cannot identify geological source-channel loss.

The median recorded-genus coverage by five nearest potential mainland source regions is **84.2% continental vs 70.2% oceanic** in all-recorded source-defined islands, compared with **96.6% vs 88.3%** in explicitly native records. These fractions reflect source-data sampling and floristic-status coverage, not proven biological access.

## 2. Full two-stage floral repertoire prediction without known island genus identities

This test goes beyond the prior [PR #6 held-out taxon-genus prediction](../20261108_heldout_mainland_prediction/README.md), which used the *observed island genus*. Here the potential source distribution over genera is derived from mainland GIFT regions **without using known test-island genera**, and its species' flower-colour/form profile distributions are mixed to predict a taxon-carrier-weighted island floral-profile frequency vector.

The source reference is learned from GIFT native mainland taxon names with every name from any of the **same 74 matched test islands** removed globally. All/Direct trait sources and the fixed 496-profile reference are kept identical to previous `fd` stages. Models are global mainland floral profile prior; source-only candidate-genus-weighted mix using global region prevalence; source-only mix using nearby genus density; and a separately labeled **known-island-genus oracle** benchmark (not an independent prospective prediction).

| Scope / model | Mean log loss for test taxon profile ↓ | Mean **floral-repertoire category-similarity-kernel** error ↓ |
|---|---:|---:|
| All / global mainland | 5.3563 | 0.020476 |
| All / source-only genus-density mixture | 5.3432 | **0.017567** |
| All / known-genus oracle | 5.2956 | 0.012699 |
| Direct / global mainland | 5.4687 | 0.024446 |
| Direct / source-only genus-density mixture | 5.4687 | **0.021290** |
| Direct / known-genus oracle | 5.4473 | 0.018572 |

The source-only two-stage mixture improves the **full original colour/form categorical-distance repertoire approximation** over the global mainland prior by an absolute **0.002909 All (14.2%)**, matched-pair bootstrap 95% [+0.001748,+0.004097], and **0.003156 Direct (12.9%)**, 95% [+0.001865,+0.004560]. The All exact profile log-loss improvement is only +0.0131 [0.0028,0.0235]; Direct is approximately **zero** (−0.000009). Thus the source-only predictor improves *assemblage-level structural composition* more clearly than *exact test taxon profile identities*.

The **known-genus oracle** outperforms source-only predictors but does not solve arrival: it is allowed to know precisely which genera were recorded on the island. Most potential mainland source genera still lack independent floral trait donors and fall back to the global prior. The approximation is *taxon-carrier-weighted profile frequency*, not true realized pollination-function occupancy or a demonstrated forecast of new colonists.

## 3. Hypervolume extension: 17 dimensions instead of only 2

Reconstructed common PCoA from the **original mean colour/form state-set Jaccard distance** across 496 profiles. The first **17 positive axes retain 80.4899%** of common positive PCoA inertia (vs **15.353%** for 2 axes). A single scalar preserves the relative contributions of those 17 coordinates rather than independently whitening small axes. Evaluated the **90% highest-density region of a 17D isotropic Gaussian KDE mixture**, with globally identical bandwidths 0.55 and 0.85, **12,000 draws × 2 independent Monte Carlo repetitions per island**.

All-scope geology-labelled complete-floral analysis: **260 islands, 130 geological continental / 130 oceanic**. Direct: **243 islands, 123 continental / 120 oceanic**. Total: **1,006 scope×island×bandwidth analyses**. The median Monte Carlo relative SE is around **2.0%** and the median absolute replicate log-volume discrepancy is about **0.04**. These are conditional numerical integration errors, not natural sampling uncertainty.

| All, 17D projected 90% KDE log hypervolume | Continental median | Oceanic median |
|---|---:|---:|
| Bandwidth h=0.55 | 19.54 | 19.37 |
| Bandwidth h=0.85 | 26.46 | 26.32 |

The enormous change in absolute log volume with bandwidth is deliberate evidence that Gaussian smoothing **dominates high-dimensional continuous-volume scale**. The 17D KDE is not biological full-category hypervolume, and 17 axes still omit **19.51%** positive inertia.

On the original 37 matched physical island pairs, raw projected volume trends negative for oceanic islands at both bandwidths and both All/Direct scopes. The four exploratory scope/bandwidth tests have BH q≈0.031 but **none survives BY** (minimum BY q≈0.064). **The critical support check is profile richness**: once each pair's difference in the number of observed distinct colour/form profiles is included, the residual geometric-origin volume difference becomes near zero. The source data pilot yielded All h=.55 richness-adjusted intercept +0.0003 (pair-bootstrap 95% [−0.0133,+0.0142]) and h=.85 −0.0141 ([−0.0353,+0.0069]). **Neither supports a unique geometric packing effect** after ordinary observed profile-count differences.

These projected volumes must **not** be compared in raw size to the earlier 2D KDE areas, and cannot be called the physical volume of pollination phenotypes. The original distance-based Rao Q/nearest-profile metrics remain the more interpretable primary evidence.

## 4. What the evidence changes

Previously we knew oceanic/continental island floras differed in exact recorded colour×form combinations (37 matched pairs, All/Direct BY q≈.00042), but not in categorical-distance-kernel multivariate deviation under the simple paired null. Genus×broad regional-presence conditioning did not establish any global excess contrast (PR #5).

The present analysis adds a **positive predictive validation**: *native mainland genus opportunity helps rank which genera are recorded, and its lineage floral repertoires predict island assemblage trait composition better than a genus-free global source model*. This integrates dispersal opportunity and floral phenotypic structure without declaring observed island records to be proven colonizations.

**Not established**: species arrival or persistence probabilities, a source-driven geographical difference after multiple-comparison controls, within-lineage evolutionary floral convergence, pollinator selection, native-only global representativeness, or a source-independent continuous functional hypervolume reduction.

A defensible working thesis is:

> **Island floral phenotype repertoires arise partly from the geographic and taxonomic constraints on which plant lineages are represented, while neither categorical distance nor high-dimensional smoothed volume alone demonstrates adaptive floral simplification.**

## 5. Execution and rights

- [Source-only genus ranker](../../src/fd/source_genus_opportunity.py) — private-input-only, observed versus unrecorded opportunity metrics.
- [Two-stage floral source model](../../src/fd/two_stage_assembly.py) — non-target donor profile mixing and known-genus oracle, private inputs only.
- [17D PCoA KDE estimator](../../src/fd/hypervolume17.py) — 496-profile public aggregate input, strict Monte Carlo quality checks, source-audited geological crosswalk recreated at runtime.
- [Source-opportunity/HV methods and caveats](../../docs/SOURCE_OPPORTUNITY_HYPERVOLUME17.md).
- Aggregate results: [source genus ranking](source_genus_rank_scores.csv), [12 paired origin-gain controls](paired_source_gain_contrasts.csv), [two-stage floral model scores](two_stage_flora_scores.csv); full corresponding provenance, CI KDE source output and sensitivity series are retained as downloadable artifacts and source SHA digests.

**Private and source-protected:** island accepted-taxon names, GIFT source entity checklists, geography source IDs, matched island identities, and original trait/source lineages. Never publish them in repository or workflow artifacts. Synthetic CI tests validate source-only model implementation. Full GIFT original-source calculations were executed locally under authorized source access; **the public 17D KDE pipeline uses only legally public aggregate floral profiles and the independent published geology**.

## Next distinct inference gate

Validation against island-native checklist completeness / repeat-detection efforts, with independent mainland lineages and dispersal opportunity, is necessary to distinguish ecological establishment from recording. Another binary-model fit to the same unrecorded labels without repeated observations will not identify colonization success.
