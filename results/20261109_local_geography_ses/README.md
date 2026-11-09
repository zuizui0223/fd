# Floral profile proximity after conditioning on nearby island repertoires (N1c)

**9 November 2026 — observed floral category support; no confirmed geological origin effect.**

[Public-source full replay #42584462875](https://github.com/zuizui0223/fd/actions/runs/42584462875) successfully ran the pinned original 496-profile state-set Jaccard D, independent geological crosswalk, and the new [N1c method](../../docs/LOCAL_GEOGRAPHIC_NULL.md). Aggregate-only raw files are in its `fd-origin-strata-analysis` Actions ZIP artifact. The final post-hoc within-island matched-N1b-versus-N1c report is an additional refinement; do not infer paired island-level changes from group median shifts alone.

## Scientific question

The earlier **N1b** null shows residual nearest-profile proximity for recorded oceanic island floras after accounting for the *regional prevalence* of each floral category. Is that simply the expected feature of sharing floral phenotypes with geographically neighboring recorded island assemblages, or does this spatial reference also leave a residual?

N1c chooses the **ten nearest OTHER eligible islands in the same broad biogeographic regime**, regardless of geological origin. Among those ten neighbors, weight each complete colour/form state-set profile by the number of neighboring islands that record it. The focal island is excluded. Draw an equal number of distinct profiles without replacement, calculate within-support mean nearest-OTHER-profile distance and off-diagonal mean pairwise distance on ORIGINAL D, and compute a sample-SD standardized reference departure.

This is a geographic profile-incidence null; exact shared taxa, observed genera, actual mainland-source opportunities, detectability and selection are **not conditioned away**. Nothing uses 17D KDE or fills unobserved hypothetical trait combinations.

## Real data: equal30 near-neighbor distances

All evidence scope (30 observed distinct profiles uniformly selected per island, 99 thinning replicates; both reference baselines use 999 draws):

| Metric | Continental (130 islands) | Oceanic (130 islands) |
|---|---:|---:|
| Median observed mean nearest-profile original D | 0.290527 | 0.275905 |
| Median N1b broad-regional expected nearest D | 0.286164 | 0.290890 |
| Median **N1b** nearest-D SES | **+0.165571** | **-0.554704** |
| Fraction with negative N1b SES | 42.31% | 84.62% |
| Median **N1c** ten-neighbor expected nearest D | 0.294432 | 0.282709 |
| Median **N1c** nearest-D SES | **-0.050883** | **-0.231046** |
| Fraction with negative N1c SES | **52.31%** | **70.77%** |
| Median N1b observed-minus-null original D | +0.004079 | -0.013566 |
| Median N1c observed-minus-null original D | -0.001216 | -0.005534 |

The **Direct** nested evidence scope is qualitatively consistent: equal30 nearest D median N1b SES **+0.163399 continental** (n=123) and **-0.543248 oceanic** (n=119); N1c SES **-0.048296 continental** and **-0.258306 oceanic** (same n). Oceanic negative N1c SES fraction is 68.91%.

For equal30 mean *pairwise* original D, All oceanic N1b SES **-0.482856** changes to N1c **-0.220519**; continental N1b **-0.078290** changes to **-0.007521**. Thus attenuation is not unique to the one nearest-neighbor statistic.

**Read this carefully:** On the identical equal30 All focal population, moving from a broad regional profile-frequency distribution to the ten geographically closest recorded island floras substantially attenuates *descriptive* oceanic category packing. It does not fully eliminate the negative median SES. Differences between group medians are **not** proof of geographic mediation or selection; the donor reference definition itself changes.

## Geographic-support audit: a major asymmetry

| Scope / origin | Eligible geography-typed islands | Nearby donor sets found | 30-profile support | Median distinct reference types | Median distance to 10th donor |
|---|---:|---:|---:|---:|---:|
| All continental | 130 | 130 | 130 | 185 | **386.7 km** |
| All oceanic | 130 | 130 | 130 | 172.5 | **1,149.2 km** |
| Direct continental | 123 | 123 | 123 | 183 | **408.8 km** |
| Direct oceanic | 120 | 120 | 119 | 173 | **1,330.2 km** |

The ten nearest eligible islands are **three times farther away** (comparing these median radii) for oceanic than for continental test islands in All, so "nearest ten" is NOT a common absolute geographic reach. The oceanic reference remains spatially broader; a lower oceanic proximity score could partly reflect this asymmetry. An explicitly matched absolute-distance reference or richer independent islands would be needed before stronger claims.

At full observed support (n not rarefied), some ten-neighbor pools have too few distinct profiles to provide n without replacement; this is **missing N1c support**, not evidence of no effect. Therefore the 30-profile analysis is the cleaner same-sample-size headline.

## Regional variation, no geological p-value

All equal30 nearest-D N1c SES medians are **-0.201 northern-midlatitude oceanic (n=48)** and **-0.226 tropical oceanic (n=76)**, but only **-0.073 continental northern-midlatitude (n=85)** and **+0.220 continental tropical (n=31)**. The southern extratropical oceanic cell contains only six islands, so its `-0.758` median is a fragile regional description, not a third independent confirmation.

No N1c geological-origin permutation p or new BH/BY claim is provided. The earlier N0/N1 32-slot spatially cautious tests were not significant (minimum BH q=0.446), and the original globally geography-disjoint 10-pair exact-profile repertoire tests did not establish a geological effect. This new N1c reference is an **explicitly post-hoc falsification**, not an additional independent test family.

## Decision gate

The defensible observation is: **the recorded floral profile neighborhoods of oceanic islands remain somewhat closer than expected under both a broad regionally frequency-weighted and a geographically-nearest-island weighted empirical reference, but are much less unusual under the latter**. This is compatible with geography/source-lineage structure, not sufficient evidence of pollinator-mediated adaptive trait compression.

To identify a distinct source-lineage-adjusted residual, the next gate is the [private-input-only taxon-held-out N2](../../docs/N2_SOURCE_DISTANCE.md), conditional on recorded genus counts and independently held-out plausible mainland native donor profiles. It requires source rights, donor-coverage and name-level predictive calibration before any real N2 inferential result is admissible.

## Reproduction and rights

```bash
python -m fd.local_geography_ses \
  --input work/island_public \
  --g0 work/fd_g0 \
  --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
  --output results/local_geography_ses \
  --neighbors 10 --draws 999 --thinnings 99
```

Only population-aggregated results and source digests are exported; source island identifiers, exact GPS points, donor neighbor membership, species and genus identities are not released. All/Direct overlap, and repeated thinnings of one island are not independent replications.
