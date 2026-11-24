# Butterfly wing and larval host traits: verified source-filtering results — 2026-11-24

[Actions 45358140791](https://github.com/zuizui0223/fd/actions/runs/45358140791) SUCCESS. Artifact 14902684326, SHA256 af11f280e929b105965bbff229a8c3e851ed1432788f33b19a4e883618df93f7. Original incidence artifact and LepTraits source digests validated. Replayed 18,283 incidence rows -> 18,280 valid normalized incidences, 2,312 distinct normalized names (one fewer than prior linkage log), 1,723 matches before duplicate audit, 1,678 retained unambiguous names. Observed valid wing-span bounds 905 species, forewing-length bounds 11, host-family count 1,088, host-account count 1,088. All source masks nested/correlated; no ecological independence assumed.

Analysis: equal-archipelago within-group slope of island species-equal mean log trait vs log2(1+ mainland distance), adjusted area and optional log record count. Minimum 5 observed trait species and >=50% trait coverage; minimum 10 or 50 butterfly record rows. 999 coherent species-level trait permutations within family/genus, preserving observed island x species incidence and trait missingness. No imputed absences. Region is broad Asian/Pacific geography, NOT strictly tropical. These are exploratory unadjusted p-values, not multiplicity-controlled evidence; too few informative archipelagos for many comparisons. Record count and documentation are not independently controlled survey effort.

## Wingspan
All/10: n58 islands, 7 informative groups, 897 species, raw slope -0.004677, record-count-adjusted -0.002454; genus p=.569/.642. Precision <=1km/10: n40, 6 informative groups, raw +.009200, adjusted +.009063, genus p=.785/.781. Non-GBIF-ID with reference/10: n35, 4 informative groups, raw -.007545, adjusted -.011331, genus p=.240/.239. No consistent wing-size filter. Forewing n too small. Wing size is not measured dispersal.

## Host plant family breadth — exploratory opposite-direction signal
| Source/threshold | Islands | Informative groups | Raw slope | Adjusted slope | Genus-shuffle p raw/adjusted |
|---|---:|---:|---:|---:|---:|
| All/10 | 75 | 8 | -.026071 | -.020830 | .004/.005 |
| All/50 | 56 | 6 | -.005521 | -.000507 | .202/.386 |
| Precision <=1km/10 | 50 | 7 | -.013777 | -.014146 | .182/.174 |
| Precision <=1km/50 | 31 | 3 | -.026200 | -.005068 | .073/.877 |
| Non-GBIF-ID with reference/10 | 42 | 4 | -.042402 | -.035159 | .002/.004 |
| Non-GBIF-ID with reference/50 | 25 | 3 | -.024272 | -.018025 | .001/.001 |

Positive predictions for broad host use on isolated islands were not supported: several panels show narrower DOCUMENTED host-family breadth, including all/10 and non-GBIF-ID/ref panels after genus-constrained shuffle and record-count adjustment. Raw/adjusted all/10 and non-GBIF-ID/ref/10 leave-one-archipelago slopes remain negative. However all/50 and precision subsets do not retain evidence. Host-account count also declines in all/10 raw and certain non-GBIF-ID/ref panels; documentation potentially confounds host breadth. Do not claim specialists colonize better or that island hosts have been lost; actual local host flora not joined. Non-GBIF-ID/ref not necessarily independent; intersection with precision has no usable comparison. Next separate documentation from host breadth on exactly matched species/islands, use host-family presence in independent floras rather than global host family count, and require strict tropical geography.

## Floral bridge
[Run 45358282275](https://github.com/zuizui0223/fd/actions/runs/45358282275) launched to reconstruct exact original GSHHG IDs and crosswalk bee polygons with >=90% mutual overlap (unique match). Public historic PR321 direct colour/form aggregates only; NOT a substitute for protected original 4,505 GIFT radial-symmetry x deep-tube analysis. No result until success.
