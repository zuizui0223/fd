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
- Calendar-month discovery and out-of-study plant prediction unresolved. Validation same recorded month 4005 pairs/177 plants/13 studies, gain .149014 [.045913,.285947], passes study-ID/plant-genus gate. Same spatial keys (3603): annual .142277, monthly .139100, difference −.003177 [−.051457,.046846]; not equivalence. Component-cluster gate fails (8 groups, required 10), both-exclusive sensitivity unsupported. Within-month phenology and sampling effort remain unresolved.
- Island vs mainland matched comparison nonestimable. No island-specific inference; island own-DB manuscript untouched.
- Hoverfly length file likely contains corrupted numeric/range representations; QC leaves zero eligible comparisons. No inferred repair.

## Figure and manuscript plan
Fig 1: same plant, sites, exclusive visitors, fixed degree and shared-partner null; show coverage and study-ID cohort split.
Fig 2: two cohort gains with crossed intervals; identical estimand and gates, report observed/null means and separate support.
Fig 3: family composition versus family-standardized residual on exact common support; visitor-genus null alongside as discovery-only limitation.
Fig 4: sensitivity summary with same endpoint (>=2 exclusives, provenance component, recorded month), explicitly distinguish unsupported and nonestimable. No figure or claim for unobserved functional equivalence.
Supplement: source pin receipts, original failed sex/caste mapping and correction, attrition, date/trait QC, full non-support ledger, fixed randomization and sample dependence, all exploratory analyses. Archive figure data from verified logs, not invented values; original artifact ZIP bytes not independently verified.

Potential target: Journal of Animal Ecology or Global Ecology and Biogeography depending on final geographic breadth and provenance audit; not yet a submission-ready manuscript. Prioritize robustness/independent source over adding trait endpoints. Apply only the previously defined half-energy distribution estimand to the larger validation cohort to distinguish size-distance conservation from distribution conservation; no replacement of failed discovery gate. FD external-only, no own data or shared protected inputs.

## Verified robustness update
| Sensitivity | Pairs | Studies | DOI/ORCID components | Gain | Crossed 95% interval | Frozen gate |
|---|---:|---:|---:|---:|---|---|
| Full validation, component x plant genus | 19618 | 29 | 14 | .109531 | [.041240,.174949] | pass |
| Exclude all discovery-linked components | 15329 | 22 | 12 | .094601 | [.026666,.165620] | pass |
| Equal component weighting | 19618 | 29 | 14 | .108746 | [.032450,.184530] | pass, descriptive |
| Novel focal plant genera | 4936 | 28 | 14 | .092491 | [−.058961,.251585] | fail |
| Novel genera and exclude discovery-linked components | 4342 | 21 | 12 | .060401 | [−.099869,.268701] | fail |
| Recorded month, study x plant genus | 4005 | 13 | 8 | .149014 | [.045913,.285947] | pass at study-ID level |
| Recorded month, component x plant genus | 4005 | 13 | 8 | .149014 | [.035758,.273074] | fail: fewer than 10 components |
| Recorded month, discovery-linked exclusion | 3676 | 8 | 6 | .097219 | [.010409,.182045] | fail: insufficient components |

Sources: [provenance audit](../results/20261128_external_euppollnet/eupolltrait_provenance_from_verified_job_log.json), [validation recorded-month audit](../results/20261128_external_euppollnet/eupolltrait_validation_month_from_verified_job_log.json). Missing identifiers and repeated taxa remain, and component grouping is not certified field-project independence. Main claim now supported under study and metadata component dependence, with study-level recorded-month sensitivity. No claim of novel-genus transport or generalized FD maintenance. No independent publisher replication. 
