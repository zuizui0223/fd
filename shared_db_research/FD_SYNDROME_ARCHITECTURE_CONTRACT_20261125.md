# FD: selective reassembly of floral architecture — 25 November 2026

## Evidence and ownership

FD main PR #76 assigns floral functional diversity, syndrome composition and pollinator channels to FD. Shared canonical botanical inputs, protected execution and reproductive assurance remain in island. This branch contains reviewed archival aggregates and a prospective analysis specification only. It does not activate protected workflows, alter Chapter 1, or merge legacy PR #69 or the conflicting reproductive-assurance PR #75.

The original paired syndrome atlas is the actual completed job log at fd run 45555736864 (pinned Git file at 919aca256109ae4e6eb64f113fd0f7cb788d694d, blob d6730e04db0020889ed3847d622eac58a7bfbeba). Audit all 80 panels, 256 template/variant records, 64 scope/pair/template groups. Preserve all results, including unsupported cells, in results/20261125_fd_syndrome_architecture/pair_atlas_all_endpoint_robustness.json.

Retrospective consistency filter: both geographic and first-source group intervals exclude zero in the same direction in all four donor/overlap variants; the standardized phi additionally requires >=200 defined-variance island pairs, >=30 groups and >=1800 valid draws per grouping. 17/64 meet the covariance screen, 10/64 also meet the phi screen (6 all-recorded, 4 literal-native). No endpoint satisfies both screens in both origin scopes. These are descriptive robustness screens, not multiple-testing-adjusted tests or independent replications. Phi uses its own defined-variance cohort. Shared species remain dependent. No template is an observed animal guild.

| Scope | Morphology proxy | Pair | Geographic covariance difference range | Phi difference range |
|---|---|---|---:|---:|
| all recorded | bird | colour × depth | +.012345 to +.014278 | +.07629 to +.09605 |
| all recorded | bird | colour × size | +.005265 to +.006865 | +.05027 to +.05481 |
| all recorded | bat | colour × size | -.008324 to -.007274 | -.07006 to -.06064 |
| all recorded | bird | form × depth | +.016867 to +.018712 | +.07190 to +.08216 |
| all recorded | specialist bee | symmetry × depth | +.008761 to +.012134 | +.03731 to +.05120 |
| all recorded | open generalist | symmetry × depth | -.009988 to -.008488 | -.03949 to -.03620 |
| literal native | specialist bee | colour × symmetry | -.016813 to -.015932 | -.05746 to -.05395 |
| literal native | specialist bee | form × symmetry | -.033051 to -.028099 | -.11535 to -.08899 |
| literal native | open generalist | form × symmetry | -.032684 to -.028423 | -.14960 to -.12919 |
| literal native | specialist bee | symmetry × size | +.009280 to +.010918 | +.07350 to +.08544 |

All ten have same-sign between-genus×grade component by the retrospective screen. Five have a same- or opposite-sign within-stratum component passing the pointwise screen; these within components must not be treated as independently supported lineage effects without variable-stratum and influence audits. The key literal-native form×symmetry pair has both specialist-bee and open-generalist negative; it is not a specialist-to-generalist substitution. The all-recorded bird proxies do not establish a native bird signal. The later color-size coherent-genus sensitivity in island #351 failed directional support; do not revive this as a primary result.

## Prospective FD mechanism-independent test: depth-conditional architecture

A more direct morphology question is whether symmetry composition changes *conditional on tube depth*, as distinct from changes in depth prevalence. Prior tropical paired results: radial×deep +10.99 pp, bilateral×deep -1.37 pp (666 islands); literal-native subset +12.63/-6.58 pp (96). These are joint frequency differences, not conditional probabilities. The archived source-relative covariance results are between-genus dominated. The 96 native islands are not a direct animal-channel study. A radial/deep increase is not evidence of nectarivore replacement.

Frozen primary estimand: tropical all-recorded island minus candidate-mainland difference in P(actinomorphic | deep) on genuinely same-species high/medium Direct, single-state symmetry×depth carriers. Retain all 3×4 categories: actinomorphic, zygomorphic, asymmetric × absent_or_open, shallow, intermediate, deep. Report P(deep), P(symmetry), P(symmetry|deep), P(symmetry|nondeep), and joint frequencies for all states. Secondary planned contrast is the difference-of-differences in radial proportion deep versus nondeep; report bilateral and asymmetric analogues to close the compositional identity. No invented missing states. Exact accepted-name joins and original pinned archive, geography, donor assignments; donor modes geo_k5 and geo50_climate10, inclusive and source disjoint from all focal island names. Native scope, all-recorded on the native-eligible island cohort and all-high quality sensitivity. Pair target and both source variants on identical support; no selection by effect or pollinator template.

Require n>=10 deep and n>=10 nondeep carriers per target and both source variants for conditional contrast; report attrition and every insufficient panel. Fixed 200 islands and 30 groups in both geographic-block and first-source group bootstrap for exploratory directional support, retain smaller native descriptive panels. Two conditional rows independently sum to one. Decompose joint frequency J_{sd}=p_s p_d+C_{sd}; C into within and between genus×exact two-axis grade via population covariance on the same roster. For conditional probabilities, compare taxonomic compositions among deep carriers, using the symmetric common-stratum identity:
  Δq = sum_common (wT+wS)/2*(qT-qS) + sum_common (wT-wS)*(qT+qS)/2 + sum_target_only wT*qT - sum_source_only wS*qS.
This is descriptive partition of a conditional distribution, not a causal decomposition. Report number of common and one-sided strata, shared mass, informative varying strata and unique genera. Preserve zero-support strata in diagnostics; no evidence of absent within-lineage effect where variation is unavailable. Reweight each accepted genus coherently across all island and mainland occurrences and both groupings; assess leave-dominant-genus sensitivity rather than treating repeated species as independent. Compare reference/grade distributions; distinguish status/quality selection from ecological sorting.

Next animal-channel prediction only after morphology result and independently documented local animal opportunity with matching geography and source support. Existing animal tests do not identify causal pollinator shifts. Do not reinterpret bird, bat or bee morphology proxies as actual pollination, ancestry, local adaptation, or evolution. No further post-hoc search for positive templates.

## Execution

The exact algebra kernel is archived in scripts/fd_depth_conditional_architecture.py. Synthetic algebra check only; full protected source integration and empirical execution require the island data custodian, as stated in main. No raw protected rows or credentials to FD. No new biological result from the proposed conditional test. The audit is actual computation on previously verified anonymous aggregate panels, not a new protected-model run.
