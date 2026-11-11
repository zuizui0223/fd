# Plant-family-disjoint prediction of GloBI documented visitor-channel composition

**11 November 2026 | exploratory, source-pinned fd PR #14.** This is an original-data follow-up to `GLOBI_CHANNEL_COMPOSITION_TRANSFER.md`. The dependent variable is **the share of GloBI-DOCUMENTED flower-interaction reference-channel units**, not real island pollinator abundance, current island visitation, effective pollen transfer, pollen limitation or seed set.

## What is held out

The previous 486-**genus** fivefold prediction showed that three coarse morphology fractions (actinomorphy, open-radial flower form, green/brown inconspicuous colour) improved genus-heldout five-channel MSE relative to the frozen source-documentation-only baseline.

The new test holds out **EVERY genus belonging to the same botanical family** together. There are **486 matched genera**, originating from **1,597 Direct/high-or-medium, strictly single-state flower morphology species**, in **88 GBIF-backed candidate families**. Mapping uses the frozen publicly recorded `island` `data/v2/staging/gbif/collected/island_taxa.csv` on commit `3bfa818ffbc7a0dea4319933c4b0a0eb41cfabb5`, blob `5ae92313517118f0dbb29d8863aa2e379ce88747`. All 486 matched genera had a recorded family assignment, but **two genera have competing historical vs current family labels** (majority label selected only as an exploratory test, and a two-genus-exclusion sensitivity is required). This is *not* an independently curated species phylogeny or authoritative WCVP full-family crosswalk.

The source GloBI v2 archive is the exact earlier `island` artifact `10328696801` (run `34790842133`); the protected Direct trait source is the earlier SHA-verified artifact `10058653212`. None of these original identities or row-level source interactions is published here.

Identical evaluated genera, five-channel GloBI reference proportions and covariates are used by every model:
1. **Evidence only:** log GloBI independent reference count, logged represented source plant count and Direct flora measurement count per genus.
2. **Three morphological proportions:** evidence + actinomorphic, open-radial and inconspicuous flower proportions.
3. **All marginal trait proportions:** evidence + 22 single-trait colour/form/symmetry state frequencies.
4. **Detailed joint profiles:** marginals + floral colour×form×symmetry joint category frequency features, 115 features total.

Fivefold family-disjoint CV, with five outcome-blind family fold assignments and fixed ridge strengths 10 and 100. Fit standardization on training families only, then project predicted five-channel compositions onto the unit simplex. Reporting uses equal genus weight; report family-balanced and leave-largest-family sensitivities. No family or genus identity enters the predictive feature matrix or is exported in results.

## Measured heldout-family result

For the first family-disjoint seed, ridge 100, the pooled five-channel MSE is:

| Model | Previous entire-genus holdout MSE | New entire-family holdout MSE | Family-heldout improvement vs evidence only |
| --- | ---: | ---: | ---: |
| Evidence only | 0.041396 | **0.041970** | reference |
| + Three coarse morphological rates | 0.039785 | **0.040451** | **+0.001519 (3.62% lower)** |
| + All 22 individual category marginals | 0.038272 | 0.041065 | +0.000905 (fragile) |
| + 115-dimensional joint floral profiles | 0.043707 | **0.049705** | **−0.007735 (worse)** |

**Three-trait gain remains positive for all five family-fold assignments** and both ridge strengths. Ridge 100 gains **+0.001240–+0.001519**, ridge 10 gains **+0.001145–+0.001473**. Excluding both ambiguously assigned-family genera yields ridge-100 gains **+0.001207–+0.001422** across five assignments. The largest individual family contains 58/486 evaluated genera; dropping that family from scored genus rows retains gain **+0.001977** in the first main split. Only **48/88** families show a positive family-specific mean gain for the first split, so the finding is heterogeneous.

A conditional **4,999 whole-family bootstrap** for the first genus-weighted gain gives **[+0.000403,+0.002716]** at the nominal 95% interval. This does not yield independent phylogenetic replicates: family histories and source publication biases may remain correlated.

**The detailed category-marginal model is not robust**: ridge-100 gain varies from **−0.005122 to +0.001658** over family-fold assignments, and its sign changes if two taxonomy-conflicted genera are excluded, partly because the balanced fold membership changes. All high-dimensional full-joint models are worse (ridge-100 gains −0.010519 to −0.005384). One cannot claim learned detailed or exclusive bee/moth/bird floral syndromes.

### Independent numerical cross-check

An independent implementation used scikit-learn's `GroupKFold(5)`, `StandardScaler`, `Ridge` (alpha100), and the same simplex projection, distinct from the custom outcome-blind family-balancing implementation. Its evidence-only MSE was **0.041835**, coarse 3-trait MSE **0.040630**, **gain +0.001205**. The joint 115-feature model worsened by **0.004984 MSE**. This confirms numerical direction across two different blocked-fold procedures, **not new independent botanical evidence**.

## Falsification by within-family floral-label reassignment

For the first family-heldout fivefold test (alpha100), shuffle the **WHOLE three-dimensional coarse morphology vector** among different genera *inside the same botanical family*, keeping each family-level floral distribution, all GloBI documentation evidence/counts and the source-channel outcome unchanged. Among 88 families, **51 have >=2 genera**, allowing **449/486** genus morphological vectors to be shuffled.

Over **999** within-family label shuffles:
- actual pooled gain **+0.001519**;
- posthoc shuffled median **+0.000650**;
- shuffled central 95% interval **[+0.000291,+0.001078]**;
- upper-tail fraction with +1 correction = **0.001**.

Thus the gain is not reproduced by genus-level morphology rearrangements that preserve each family's overall coarse morphology. This is a useful **family-conditioned predictive placebo**, but still not a fully conditioned phylogenetic/source-provider null, ecological selection test, or actual-island pollinator measurement. The 37 families with only one sampled genus necessarily cannot be shuffled.

## Robust scientific claim and limits

**Supported at predictive/evidence level:** directly recorded, coarse floral morphology carries a **small** amount of additional information about the relative global GloBI-DOCUMENTED five visitor-channel references, even when test genera come from entirely untrained botanical families and when shuffles preserve family-level morphology. The original genus-heldout finding is therefore **not solely a same-family sharing artifact**, and it is not trivially the broad literature-reference-count gradient.

**Not established:** a classification of actual bee, butterfly, moth, fly or bird pollination of a plant species or island; moth vs butterfly separation (GloBI channel is combined Lepidoptera); visitor abundances or effective pollen transfer; independent floral syndrome evolution; pollen limitation mediation; pollinator-driven adaptation. Unrecorded GloBI channel evidence is NOT biological absence. The 928-island GBIF channel evidence has only 33 strict mixed-status islands, and 0 of 20 region×guild contrast support cells pass the earlier ecological gate. Whole-family taxon mapping is from a candidate GBIF table, and the matched trait-observed genus subset is highly selected. Many genera can share evidence providers across family boundaries. Neither bootstrap nor 999 feature permutations creates new external biological observations.

**Next identification gate:** use independent real island visitor samples with explicit study effort, or a separate visit-effectiveness/pollen deposition database; before any causal claim audit WCVP/phylogenetic family labels and GloBI original publication-source overlap. The 486-genus GloBI predictive endpoint remains a documentation endpoint, not ecological service.

## Reproduction and rights

- Public, pure `src/fd/globi_family_holdout_cv.py`: unique-genus/family validations, five balanced whole-family test folds, leak-free ridge, simplex constraint, equal-family/cluster sensitivity and within-family three-trait label placebos.
- Synthetic regression tests in `tests/test_globi_family_holdout_cv.py`.
- Public nonidentifying aggregate tables: family-vs-genus heldout scores across five family-fold assignments, independent GroupKFold check, whole-family conditional intervals, per-channel output and within-family shuffles. The **genus-to-family mapping and original species, GloBI reference/interaction and island rows remain in a private local source-audited workspace**.
- Original GitHub Actions source artifact and blob identifiers allow auditing under the original data rights. Public GitHub CI checks synthetic algorithms and confidentiality only; it **does not independently rerun the source-restricted original-data calculations**.
