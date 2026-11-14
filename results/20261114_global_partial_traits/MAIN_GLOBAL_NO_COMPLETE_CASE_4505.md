# 14 November 2026 — Global island flora WITHOUT floral complete-case selection

**Original full GIFT island-flora source actually executed**: [island Actions #43586493511](https://github.com/zuizui0223/island/actions/runs/43586493511) both numerical seeds SUCCESS. Original anonymous aggregate artifacts #14110797789 (20261114), #14110252694 (20261115). This supersedes the obsolete practice of treating only 928 triple-complete floras as the global macro baseline.

## The actual global population and why 8,264 != 4,505 != 928

- Frozen GSHHG source global geographic universe contains **8,264** non-continental island units as of the audited chapter-1 geography. The independently classified Roeble/Zizka geological origin identifies **434** as truly geologically continental or oceanic. Disconnected since the Last Glacial Maximum does NOT prove geologically oceanic.
- The restored original GIFT status-flora member list, irrespective of individual floral trait availability, contains **4,505 distinct source-island botanical floras**; **one** original island flora does not join the separately audited GSHHG crosswalk. That flora remains represented in the global region-unresolved group; it must NOT be discarded nor called zero botanical diversity.
- **425** original botanical floras are represented among the 434 independently geologically typed islands. Their use for geology-specific comparisons is a SEPARATE explicitly bounded stratum.
- The historically restrictive *three different species-direct flower traits all jointly known in ≥30 plants per island* subset remains exactly **928** original island floras, as a backward-compatible sensitivity result only.
- Main data structure is **4,505 islands × 2 original source floristic-status scopes (all-recorded and strict original source-native) × 13 endpoints (10 single traits and three two-trait pairs) = 117,130** original island/status/endpoint ledger rows before anonymous regional grouping. They are NOT 117,130 independent island floras.

## Ten independent species-Direct botanical trait rosters

Each original accepted species×trait must have unambiguous one-valued high/medium Direct evidence; multiple proposed states and conflict tokens are excluded rather than imputed. The full original species roster is not reselected to require complete colour/form/symmetry.

| Original directly recorded trait | Distinct source-Direct species worldwide | Number of distinct ORIGINAL island floras with ≥5 source-Direct species for this individual trait |
| --- | ---: | ---: |
| Flower primary colour | **26,553** | **3,132** |
| Floral form | **9,173** | **2,806** |
| Floral symmetry | **32,465** | **3,145** |
| Flower size class | **7,373** | **2,781** |
| Tube depth class | **4,855** | **2,549** |
| Inflorescence display | **21,784** | **3,212** |
| Self-incompatibility label | **6,857** | **2,958** |
| Autonomous selfing capacity | **1,593** | **2,411** |
| Mating system | **3,003** | **2,376** |
| Cleistogamy | **1,323** | **2,024** |

Counts above use **all-recorded floras** (native, introduced and currently unresolved plant statuses are NOT silently treated as indigenous). Strict original GIFT source-native subsets are reported separately, but considerably smaller and geographically uneven. No WCVP regional-native record is promoted to native status on a specific island without independent source evidence.

## True pairwise evidence — NO THIRD trait required

| Actual original single species must carry | Number of original island floras with ≥5 directly recorded accepted species |
| --- | ---: |
| Flower colour + form | **2,324** |
| Flower colour + symmetry | **2,699** |
| Floral form + symmetry | **2,611** |
| The historical ≥30 joint triple flower category criterion | **928** |

Pairwise association is the exact categorical one-pair source-shuffle departure, using only the TWO relevant source-Direct traits. The Rao Q is the exact finite-population expected five independently sampled original accepted distinct species value, derived from the source population same-category probability without inventing interpolation or stochastic KDE hypervolume. A 5-species threshold is specific to the trait endpoint; 0–4 actual Direct species are retained as documented missing/insufficient opportunity, not estimated as zero variation. A genus-equal alternative analytically averages category frequencies within each recorded genus before computing between-genus source expected categorical mismatch. This equalizes genus representation but does NOT supply phylogenetic ancestry.

## Geological origin — post hoc source-relative exploratory patterns, NOT confirmed adaptation

Geological-type subset has 425 botanical floras, not 148 triple+30 source-restricted floras; 999 physical spatial-block and first-GIFT source-region cluster-bootstrap repeats were computed in TWO numerical seeds. These are the SAME original species and island observations (not independent biological replications). 624 strongly correlated post-hoc origin/outcome/scope and cluster comparisons were explored: **do not treat unadjusted intervals excluding zero as genome-wide familywise evidence or biological causal effects**.

Examples that retained the same sign in BOTH independent numerical seeds, under BOTH cluster groupings:

1. Tropical all-recorded **autonomous-selfing CAPACITY category diversity** (source-relative, genus-equal finite five-species Rao): oceanic minus continental mean **−0.06006**, geographic-block CI approximately **[−0.097,−0.026]**, first-GIFT source CI **[−0.106,−0.013]**, raw n=119 oceanic /37 continental. This estimates diversity of DOCUMENTED selfing capacity categories, NOT the actual proportion of species selfing, reproductive assurance, local autonomous seed set or selection.
2. Tropical strict original-source-native **flower-size CLASS diversity** (source-relative genus-equal finite five-species Rao): mean approximately **+0.04513**, block CI **[+0.020,+0.071]**, GIFT source CI **[+0.017,+0.068]**, n=77 oceanic /29 continental. This points the opposite direction from an assumed universal oceanic simplification, yet native labels are incomplete and this is selected among many outcomes.
3. Tropical all-recorded **flower-colour × floral-form categorical pairing departure** (GIFT native mainland source-relative): oceanic-minus-continental ~**+0.00617**, geographic block 95% **[+0.00002,+0.01228]**, first-source CI ~**[+0.00007,+0.01343]** on **120 oceanic vs36 continental** island floras. Numerically fragile at zero and multiple-comparison sensitive. This pairing departure is NOT pollinator behaviour, mutual information or evidence of independent syndrome evolution.

The broader geology-original null of morphology and pollinator niches remains uncertain: smaller independent-geography origin-matched island pair comparisons previously yielded only 3 pairs per major regime. Source-relative floristic contrasts show observational botanical *assemblage composition*, not real dispersal arrival or historical extinction.

## Mandatory next audit: data coverage and missingness

Reported original species-wide Direct trait coverage is uneven by class and geographic background. Record-rate differences may induce the appearance of floral category changes. The follow-up version expands the same source run to geological origin contrasts in **observed Direct record fraction per eligible island flora**, including zero observed Direct traits only when that island contains actual recorded plants, and NA when native-subset denominator is absent; 0 must never be imputed for missing true phenotype.

It will be necessary to show whether the exploratory associations survive taxonomic source/coverage weighting, geographic robustness, suitable multiplicity control and independent native source. The current descriptive results do not prove adaptation, selection, causal pollinator reduction or the island syndrome.

## Scope and reproducibility

- Full global analysis [fd source code](../../src/fd/global_partial_trait_flora_fd.py) and [tests](../../tests/test_global_partial_trait_flora_fd.py); original private [island runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_global_partial_trait_fd.py), [source workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-global-partial-traits.yml).
- Do not conflate native GIFT mainland candidate pools with evolutionary source lineages, or all-recorded historical floras with actual locally indigenous species. Keep old 928-triple Rao as a nested sensitivity, not as an eligibility prerequisite for marginal and pairwise global floral analysis.
- Do not export original source accepted plant names, island IDs, geography, GIFT source IDs or original botanical trait records from the protected runner. Only anonymous regional counts and conditional aggregate estimates can be published. Keep PRs DRAFT.
