# Original EltonTraits × pre-N1 GloBI pollinator-guild construct audit (2026-11-19)

**Status:** successfully inspected two actual, independent taxon-level original-source archives; anonymous aggregate only. **This is NOT an island-level guild–flower association, a pollination syndrome shift, or an evolutionary adaptation test.** The original 4,505 GIFT botanical island floras are not joined here.

## Provenance and exact independent source gate

- **Bird diet:** Wilman et al. (2014) EltonTraits 1.0, Ecology, https://doi.org/10.1890/13-1917.1, Figshare collection https://doi.org/10.6084/m9.figshare.c.3306933.v1; original BirdFuncDat.txt file **5631081**, successfully fetched from official Figshare in [fd Actions #44544141410](https://github.com/zuizui0223/fd/actions/runs/44544141410). Raw original SHA-256 `97216eb1797da077169ebb1ebea275db293b09fc62f8bb8911f9beb98c50d321`. Precisely **9,995 raw lines**, **2 unnamed/missing-diet rows quarantined**, **9,993** accepted name-unique finite nectar diet fractions, **1,270** fractions >0, 8,723 equal to 0. Values are 0–100 original `Diet-Nect` converted to fractional [0,1], **not field nectar visits**. File is CC BY 4.0. No eBird observations/taxonomy are used.
- **Independent visitor evidence catalog:** protected source original artifact `zuizui0223/island` **#10291643834** (`channel_taxon_catalog_confirmatory.csv`), exact ZIP SHA-256 `c946e36cd04cb319fc63a5bde90534c1cf3c7f41480fc672638f85f5b2a89b02`. It was built using GloBI v0.9 original source DOI https://doi.org/10.5281/zenodo.20546682, source SHA-256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`. The source receipt certifies **6,928** confirmatory taxa across five potential visitation channels and `uses_focal_plant_traits=false`; this original NEE catalog predates the current floral test. GloBI `n_pollination_rows>0` is a documented pollination **claim**, not proof of effective pollen service.
- **Taxonomic match:** literal original accepted scientific name after trimming whitespace; no fuzzy joins, forced family match, same-genus substitutions, taxonomic synonym rescue, or outcomes used to choose a name. Every unmatched name remains unknown. The public result contains no original source species, island, location, occurrence, plant or provider identifiers.

## Actual original source result

| Independent bird catalog diagnostic | Original source count |
| --- | ---: |
| Confirmatory GloBI flower-visiting bird species | **370** |
| Exact species name matched to EltonTraits | **281** (75.9% of 370) |
| Matched and Elton nectar/pollen fraction > 0 | **179** |
| Matched and Elton nectar/pollen fraction = 0 | **102** |
| GloBI bird catalog not name-matched to Elton | **89** |
| Positive Elton diet **and** ≥1 GloBI pollination claim | **121 / 179** |
| Zero Elton diet **and** ≥1 GloBI pollination claim | **22 / 102** |
| All GloBI confirmatory Lepidoptera | **2,360** |
| Of those, accepted Sphingidae taxa | **92** |
| Sphingidae with ≥1 GloBI pollination claim | **58** |
| Sphingidae with ≥1 GloBI flower-visitation record | **63** |

For the **281 actually name-matched birds**, nonzero recorded nectar/pollen diet is associated with recorded GloBI pollination claim: the post-hoc 2×2 species count is `[[121, 58], [22, 80]]`, exploratory Fisher odds ratio **7.586** and two-sided p **7.06×10^-14**. This is **a source-documentation association only**: the 281 species are phylogenetically clustered, GloBI selects previously documented visitors, primary literature may overlap, and the two evidence sources and observation effort may be systematically non-independent. Do not treat this p-value as a test of island ecology, visitation rate, or selection.

**Conservative guild-definition uncertainty:** 179/370 = **48.4%** of the bird-visitor source catalog has a cross-dataset *documented* positive nectar/pollen diet. At most (179+89)/370 = **72.4%** could be diet-positive if all 89 unmatched names were positive; that upper scenario is not a prediction. The 102 Elton recorded zero fractions describe 2014 **coarse dietary attribution**, not biological inability to carry pollen. Therefore never use all 370 GloBI bird visitors as interchangeable strong nectarivores, and never turn unmatched or diet-zero names into spatial guild absence.

**Sphingidae requirement:** the 92 catalog taxa are Sphingidae-associated species-level flower-interaction candidates. There is still **no independent adult nectar-feeding trait validation** for all 92, and GloBI pollination assertions for 58 do not prove effective pollination.

## Next true ecological test (still blocked)

To use these sources in the frozen [GUILD_SPATIAL_ALIGNMENT_CONTRACT](../../docs/GUILD_SPATIAL_ALIGNMENT_CONTRACT.md), map **continuous nectar use weight × independently qualified bird spatial occupancy** for every native resident breeding bird with a reconciled BirdLife range or other independently validated range. The original EltonTraits taxonomy predates BirdLife2020; the 89 unresolved source names already demonstrate why unreviewed synonyms and missing range polygons cannot be treated as zeros. For Sphingidae, distinguish independently documented **adult** nectar feeders and effort-adjusted occurrence/SDM coverage, rather than substituting all Lepidoptera or 92 GloBI candidate species as definite nectar-pollinating moths. Distinct GIFT source opportunities, taxonomic genus-expected flower baseline, original Direct joint trait completeness, spatial blocks and first-source uncertainty remain mandatory.

**Do not merge an ecological positive claim to main** before the original botanical-island × mapped-guild join and model are run. This is a valid, completed *source measurement/construct* advance but not yet the 4505-island syndrome test.

Original aggregate (safe to publish): [ORIGINAL_AGGREGATE.json](ORIGINAL_AGGREGATE.json). Reproducible code/source guards: `src/fd/elton_globi_guild_bridge.py` and synthetic tests. The two original ZIP files and all species identifiers remain outside the public Git repository.
