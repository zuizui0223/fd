# Archipelago-level held-out prediction of island floral category repertoires

**9 November 2026 — exploratory fd extension.**

## Main problem

The prior 37 geology-matched island pairs belong to only five spatially connected source groups. Their previously BY-positive sign-swap comparison is NOT a verified set of 37 independent geological replicates. Strict geography-disjoint re-matching (10–13 pairs) also did not establish a repeated geological-origin effect.

The next independently testable question is not geological treatment, but **whether the floras of other islands in a named archipelago predict the recorded colour × form composition of a held-out island** more accurately than the remainder of the same climatic geographic regime.

## Input eligibility

Use the original rights-cleared public island × flower-profile carriers from island PR321, the same 496 original colour/form state-set reference, and the independently geocoded Weigelt archipelago labels. All and Direct must both contain at least 50 complete-floral name units. A source group must have at least **four islands within one original analysis regime**. Missing Weigelt labels — including literal NA — do not imply a source group and must not cause artificial shared membership.

An earlier outcome-blind eligibility audit of these archives found **711 shared All/Direct islands with source archipelago labels**, **112 unique named groups**, and **45 group × regime cells having at least four islands**. The grouping is observational and Weigelt labels can be broad rather than phylogenetic source regions.

## Source-only prediction and near-distance control

For each eligible island target, remove the **entire target archipelago** from the regional reference profile pool. Define the baseline as the island-equal average of the other-archipelago profile vectors, within the same original geographic regime. Every vector has sum one across original *distinct* recorded profile IDs; do not expand multistate descriptions.

- **Same archipelago**: take the three geographically nearest *other* islands in the archipelago/region and average their profile vectors. Blend that prediction 50% with the other-archipelago regional baseline. The target island is never a donor.
- **Distance-matched other archipelago**: find three different-archipelago donors from the same regime with geographic distances to the target matched, by a one-to-one Hungarian assignment on log(1+distance), to the three in-group donors. Every distance discrepancy must be at most ln(3), and the third in-group neighbor must be no farther than 1500 km. Blend identically with the regional baseline.
- If no matched external group exists, retain the island for the broad regional baseline comparison but **exclude it from the distance-matched contrast**, never imputing a favorable control. A prior labels-and-distance-only audit found about **51 islands from 18 named source archipelagos** in the strict near-distance control.

The model does **not** know test island profile values before prediction. It may train on taxa shared by different islands; therefore this is leave-one-ISLAND-out, not leave-one-SPECIES-out, leave-one-archipelago-out or historical source-lineage prediction.

## Outcomes and inference

Score predicted profiles against the heldout island using two different metrics across the **same 496 original state-set profiles**:

1. Exact-profile squared Euclidean error: (y−p) transpose I (y−p).
2. Original categorical-similarity-kernel error: (y−p) transpose (1−D) (y−p), where D is the original mean Jaccard distance of recorded colour/form state sets, NOT projected 2D/17D KDE.

The comparison statistic is **baseline error minus within-archipelago error** (more positive means the group label improves prediction), or **distance-matched other-archipelago error minus within-archipelago error** (positive means archipelago information beyond matched geographic proximity). The latter is the decisive *location-caliper* control.

Before uncertainty calculation, average eligible island-level gains by the **named archipelago** across all regions. Use the group, not the island, as the unit in a 9,999-draw directional sign-swap sensitivity and 4,999-draw archipelago bootstrap. Compare 2 scopes × 2 metrics × 2 contrasts as a single eight-slot BH/BY family; an untestable slot receives p=1 rather than silently reducing multiplicity.

Group sign exchangeability may remain imperfect due to shared regional donors, geography, shared species and administrative breadth of source labels. This is a **source-assemblage prediction diagnostic**, not a randomized historical source experiment. A positive within-group gain without a matched-distance gain is compatible with ordinary spatial autocorrelation. Significant exact-profile gains without original-kernel gains imply taxonomic category redistribution rather than verified functional niche differentiation.

## Data rights and reproducibility

The implementation reads only the pinned public aggregate table and independent source crosswalk. It writes aggregate source support, predictive errors, eight statistical comparisons and source SHA256 digests, **not** individual island identities, archived archipelago labels, donor assignments or restricted species-to-trait observations.

Run within the existing source-pinned fd origin workflow:

    python -m fd.archipelago_cv --input work/island_public --g0 work/fd_g0 --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv --output results/archipelago_cv --permutations 9999 --bootstrap 4999

This work deliberately shifts the paper away from an unsupported general oceanic-island adaptive trait-volume contraction to the measurable, geography-dependent **predictability of recorded floral assemblage composition**. Actual colonization, establishment, plant-source ancestry, pollination effectiveness and species-wise detection processes remain unknown.
