# fd — Equal donor opportunity separates floral symmetry from SC/SI predictability

**2026-11-11 · real source-pinned exploratory analysis · PR #13.** This is a prediction of recorded binary trait states conditional on **known genus**, not a test of realized pollination performance, trait evolution, or independently replicated source ancestry. Only aggregate output is released.

## What changed from the previous source-disjoint test

The previous stringent global-name/source holdout found persistent genus-associated predictability for floral actinomorphy versus zygomorphy, while Direct self-compatibility (SC/SI) was unstable. But it used **different sets of training species and evaluated species for the two traits**, giving vastly different donor opportunities.

We now require EVERY training mainland species AND EVERY evaluated island species to have **BOTH** strictly single-state original traits documented. Therefore the donor names, evaluated species, evaluated genera, number of trained genus examples, heldout GIFT regions, and global island-name exclusion are the SAME for both predictions. Each trait retains its own observed state and trait-specific global baseline.

This audit uses the same original 1,419 G0 study islands and 9 within-region source-overlap components with >=40 islands. Before training, every name recorded on ANY of the 1,419 island floras (**110,587 distinct accepted names**) is excluded from ALL mainland donor entities, including unmeasured island taxa; the full five-nearest native GIFT mainland source-entity union linked to each evaluated component is excluded as well. Five-source GIFT labels describe **possible checklist opportunity**, not actual ancestors.

## Main result: controlled opportunity, different trait prediction

At tau = 1 pseudo-donor name, all-recorded island targets:

| Evidence | Floral symmetry Brier improvement vs global | SC/SI Brier improvement vs global |
| --- | ---: | ---: |
| All, snapshot ledger | **9/9** positive component comparisons | **9/9** |
| All, source-adjusted ledger | **9/9** | **9/9** |
| Direct, snapshot ledger | **9/9** | **7/9** |
| Direct, source-adjusted ledger | **9/9** | **5/9** |

All/Direct and both ledgers are overlapping sensitivities, **not 36 independent replications**. The nine components are geographically separated *within their broad regions* but still share species, genera and donor species at broader scales. There is one globally dominant donor-overlap network, so these signs are not independent evolutionary hypothesis p-values.

In the strictest Direct/source-adjusted cohort, the common mainland donor pool is only **54–55 distinct known-for-both-traits species** per test component. The median evaluated island genus-donor coverage is **5.8%** (mean) and the median number of covered test genera is **19** (minimum 4). Both traits have EXACTLY the same target and donor support.

For the 9 Direct/source-adjusted comparisons, mean across component-specific name-equal Brier improvements is **+0.00771 symmetry** versus **+0.00114 SC/SI**. The **symmetry minus SC/SI improvement is positive in all nine components** (mean difference +0.00657), and remains positive in all nine after giving each evaluated genus equal weight. However, the top three donor-covered genera account for a median **62.7% of the covered species**. Removing these outcome-blind largest genera retains positive symmetry improvement in **8/9** and SC improvement in **8/9**; the paired difference remains positive in **7/9**. This weakens any claim of broadly independent genus-level replication, especially at northern high latitudes where only **four genera** have a matched mainland donor.

### Native-record-only target sensitivity

The same all-island-name and source-region exclusions are retained. Seven of the nine source components have at least 50 known-for-both-traits **source-backed island-native** targets in the relevant strict scopes. Direct symmetry prediction remains positive **7/7** in each original ledger; SC/SI **5/7** in each. Native status comes from the original recorded source backing; unresolved is NOT relabeled native.

## Is SC/SI merely under-observed?

An additional **non-source-disjoint**, name-level five-fold donor-number simulation evaluates the full both-traits-recorded original species catalog (Direct/source-adjusted: **3,895 species**) with exactly the SAME randomly selected donor names for SC/SI and symmetry. Each test species is held out in its own test fold. A donor-size condition is repeated 150 times per fold; those 750 correlated draws are **data-availability sensitivities**, not evolutionary replications.

| Matched donor species per fold | Donor-genus test coverage (median) | SC/SI mean Brier gain | Symmetry mean Brier gain | SC / symmetry mean skill fraction of global Brier |
| --- | ---: | ---: | ---: | ---: |
| 55 | 12.8% | +0.00127 | +0.02041 | 0.60% / 8.68% |
| 100 | 19.8% | +0.00336 | +0.03347 | 1.59% / 14.36% |
| 200 | 30.3% | +0.00808 | +0.05433 | 3.84% / 23.42% |
| 500 | 47.3% | +0.01928 | +0.08999 | 9.18% / 38.91% |
| 1,000 | 60.3% | +0.03031 | +0.11916 | 14.44% / 51.57% |
| 2,000 | 71.7% | +0.04233 | +0.14656 | 20.18% / 63.47% |
| All eligible training names | 77.6% | +0.04996 | +0.16145 | 23.81% / 69.93% |

SC/SI improves with increased same-genus donor observations, showing that its weak sparse-source performance is **partly a data-availability issue**. At the SAME amounts of donor evidence, however, symmetry has stronger conditional genus-associated predictability across the entire diagnostic curve. This is consistent with the previously observed greater within-genus variability of recorded SC/SI. Because catalog evidence providers, trait ontology, and studied lineages differ, this does **not** estimate separate evolutionary rates or selection strengths. This simulation does **not** estimate the number of new independent GIFT SC/SI records that would be required.

## Exactly where new SC/SI observations would help

The strict all-1,419-island name exclusion and source-region holdout yield **104–105 Direct SC/SI donor names** per component (compared with 54–55 names already known for BOTH traits). For the Direct/source-adjusted ledger, the current fraction of evaluated SC/SI names with a same-genus donor ranges from **4.2% to 12.7%** (mean ~7.9%).

Two outcome-blind candidate groups remain in the eligible mainland flora: **185–195** taxon names per component have SC/SI recorded in All but not Direct, and **751–831** have Direct symmetry but no Direct SC/SI state. These are overlapping candidate lists across components, NOT thousands of independently available new studies.

The theoretical upper bound on same-genus evaluation coverage if every candidate received newly validated Direct SC/SI evidence is:
- Review/verify **All-only SC/SI** candidates: mean **18.9%** coverage
- Obtain true SC/SI data for **Direct-symmetry-only** candidates: mean **21.1%**
- Both candidate classes combined: mean **29.1%**

These are upper bounds on *donor availability*, not modeled prediction gains; unknown SC/SI cannot be set to either state. All evidence must not be auto-promoted to Direct, and source-specific redistribution restrictions remain in force. Independent published reproductive studies and experimentally measured SC/SI are the correct data-acquisition priority.

## Verification and inference limits

- All three original ZIP SHA256 hashes were independently verified against fd's frozen provenance:
  - progressive original: `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`
  - independent trait-source: `5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd`
  - fd G0/G1 public comparison: `2efa7d2bb40adc60cc7d22b3c0c1043eb63c7c608e492fb0a422a14b3a563801`
- An independent implementation of the previous **72** strict source-disjoint single-trait results reproduced all donor/test denominators exactly and all Brier gains with maximum absolute discrepancy **<1e-15**. This validates cohort reconstruction and calculations, NOT independent biological evidence.
- The equal-opportunity rerun yields **128** component × ledger × evidence × all/native × trait aggregate results; the cross-component summary includes **16** rows. All target/donor identity-level data remain in the locally protected workspace and are not committed or distributed.
- Public helper `src/fd/trait_pair_opportunity.py` with tests `tests/test_trait_pair_opportunity.py` checks identical species identity for both traits, impossible states, duplicate names and global donor-name leakage. GIFT region exclusion and source-rights compliance must be verified by the private cohort builder; synthetic CI cannot replay licensed raw inputs.
- The 9 components are not globally independent, Native-only is a selected status subset, and trait categories do not directly measure floral mechanical fitness or reproductive seed set. Further optimization after inspecting the results is exploratory.

**Conclusion:** recorded floral symmetry is more consistently genus-predictable than strict SC/SI **even at exactly matched donor and evaluated species opportunity**, but SC/SI gains grow with more valid donors. The next scientifically useful step is **new independent Direct reproductive compatibility evidence**, not increasing the number of overlapping island pairs or KDE hypervolume tests.
