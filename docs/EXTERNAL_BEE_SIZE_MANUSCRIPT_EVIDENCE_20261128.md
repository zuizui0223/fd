# Evidence freeze — spatial conservation of bee body size despite visitor turnover (2026-11-28)

Working title: **Bee body-size similarity persists across spatial turnover of flower visitors**.

Central observation: distinct bee species visiting the same plant at separate sites have smaller absolute log2 female/worker intertegular-distance differences than expected under a local binary fixed-margin null conditioned on observed partner overlap. This measures similarity of species-level morphology, not actual individual size, pollen transport, fitness or causal compensation.

## Main result (two nonoverlapping study-ID cohorts)
| Cohort | Contributing study IDs | Focal plant species | Exact comparisons | Observed | Null | Gain | Crossed 95% interval |
|---|---:|---:|---:|---:|---:|---:|---|
| Broad four-order discovery | 14 | 178 | 1261 | .628644 | .739293 | .110649 | [.024189,.191188] |
| Bee-inclusive validation | 29 | 405 | 19618 | .545831 | .655362 | .109531 | [.039549,.172565] |

Do not add unique species across groups. Validation metadata allowed 34 non-discovery study IDs; 29 contributed. 43 contributing study IDs are not certified independent projects. Both cohorts same EuPPollNet and EuPollTrait publisher data. The original female-worker correction was made after original outcome inspection; validation scope frozen before validation outcome but not a pristine external holdout. No population equivalence test; similar point estimates not proof of equality. Main citations to verified aggregate ledgers: [corrected discovery](../results/20261128_external_euppollnet/eupolltrait_size_qc_family_from_verified_job_log.json), [study-disjoint validation](../results/20261128_external_euppollnet/eupolltrait_disjoint_validation_from_verified_job_log.json).

## Secondary evidence and limits
- Exact common-support family accounting: discovery 1117 pairs, total .094756 = .052563 family composition + .042194 standardized residual; validation 17133 pairs, total .109945 = .078054 + .031891. Crossed intervals for standardized residual positive in both, discovery lower bound .000270; arithmetic partition, not causal mediation. Validation secondary >=2 exclusive species per side (6220 pairs) gain .096933 [.033572,.156077].
- Discovery visitor genus-coherent label audit: p=.506, residual .009047 [−.033523,.052011]. No claim of species-specific morphological sorting beyond visitor genus.
- Distribution half-energy and three-axis common-carrier morphology unsupported under crossed uncertainty. No generalized maintenance of functional diversity. Chacoff single-study opposite-direction morphology result retained as limitation, not pooled.
- Calendar-month discovery and out-of-study plant prediction unresolved. Validation temporal opportunity audit and provenance component audit pending; final claims conditional on verified results.
- Island vs mainland matched comparison nonestimable. No island-specific inference; island own-DB manuscript untouched.
- Hoverfly length file likely contains corrupted numeric/range representations; QC leaves zero eligible comparisons. No inferred repair.

## Figure and manuscript plan
Fig 1: same plant, sites, exclusive visitors, fixed degree and shared-partner null; show coverage and study-ID cohort split.
Fig 2: two cohort gains with crossed intervals; identical estimand and gates, report observed/null means and separate support.
Fig 3: family composition versus family-standardized residual on exact common support; visitor-genus null alongside as discovery-only limitation.
Fig 4: sensitivity summary with same endpoint (>=2 exclusives, provenance component, recorded month), explicitly distinguish unsupported and nonestimable. No figure or claim for unobserved functional equivalence.
Supplement: source pin receipts, original failed sex/caste mapping and correction, attrition, date/trait QC, full non-support ledger, fixed randomization and sample dependence, all exploratory analyses. Archive figure data from verified logs, not invented values; original artifact ZIP bytes not independently verified.

Potential target: Journal of Animal Ecology or Global Ecology and Biogeography depending on final geographic breadth and provenance audit; not yet a submission-ready manuscript. Prioritize robustness/independent source over adding trait endpoints. FD external-only, no own data or shared protected inputs.
