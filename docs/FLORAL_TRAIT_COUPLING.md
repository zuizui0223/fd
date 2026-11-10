# Floral trait-source prediction: marginal colour, marginal form and joint combinations

**2026-11-10. Status: mathematical redefinition + source-pinned public-kernel audit + synthetic tests; private N2 joint-prediction re-estimation is NOT yet a verified empirical result.**

## Previously overlooked identification problem

The original 496-profile colour/form dissimilarity in fd is

    D(a,b) = 1/2 [1 - J(C_a,C_b)] + 1/2 [1 - J(F_a,F_b)]

where C and F are complete recorded state SETS, not artificially expanded
Cartesian-product phenotype states. Thus its positive similarity is

    Kadd(a,b) = 1 - D(a,b) = [Kcolour(a,b) + Kform(a,b)] / 2.

A quadratic *proper* kernel score on categorical probabilities p and observed
profile y is

    S_K(p,y) = p^T K p - 2(Kp)_y + K[y,y].

This score is **linear in K**. Consequently, **ALL** original Kadd source
prediction gains equal exactly one half the sum of colour-only and form-only
prediction gains.

**Important:** Kadd does NOT identify whether a native-mainland source
improves prediction of the *combination* of a flower colour and form. For
example, half red+tubular and half blue+open has identical separate colour
and form marginals to half red+open and half blue+tubular, even though
the two joint flower phenotype distributions do not overlap.
Kadd-score treats these two forecasts identically in expectation.

## New testable kernels

On the *same frozen 496 intact profile IDs* define

- Kcolour: Jaccard similarity of the **recorded colour state sets**.
- Kform: Jaccard similarity of the **recorded form state sets**.
- Kadd: their equal-weight average. Its score must reproduce the previous
  source-gain result to machine precision or the analysis must abort.
- Kjoint: the elementwise (Hadamard) product Kcolour * Kform.
  If both component Jaccard kernels are positive semidefinite, their
  product is positive semidefinite by the Schur product theorem.
  This is a **joint-pairing-sensitive proper score**: a model predicting
  the right separate marginals but wrong colour/form pairing can score
  substantially worse than the correct pairing model.

A direct independent **pivoted-Cholesky rank audit** of the frozen
496-profile public reference (Git blob `982ecb538eea67fff290c6025be7f5c5a093596d`,
rank threshold 1e-8) gave:

| Proper similarity kernel | Numerical rank among 496 profiles |
|---|---:|
| Colour-state-set Jaccard Kcolour | 31 |
| Flower-form-state-set Jaccard Kform | 46 |
| Original equal-weight Kadd | 76 |
| Colour/form Hadamard product Kjoint | **496 (full)** |

The reference contains **31 distinct recorded colour state sets** and
**46 distinct recorded form state sets**. The rank of additive K is
exactly 31+46-1=76 because marginal feature space overlap contributes
the common constant function. Kjoint has **full rank on all 496 observed
complete trait profiles**, so the associated quadratic K-score is
**strictly proper for the full probability distribution over those 496
categories** (subject to independent PSD/eigenvalue confirmation in
public CI). This changes which estimand we are capable of evaluating:
the old Kadd cannot identify exact colour/form pairing probability
structure, whereas Kjoint can.

**The mathematical full rank is not an empirical ecological effect.**
No donor-source predictive improvement or selection inference has yet
been verified with the joint kernel. The public GitHub Action additionally
calculates symmetric eigenvalues and SHA256 for source-pinned confirmation.

**Do not overclaim:** Positive Kjoint improvement does not automatically
prove an ecological 'colour–form interaction'. Joint predictions change
when colour marginals or form marginals change; a pure statistical
interaction requires a further carefully constructed **marginal-preserving**
source null (not a naive difference of non-identifiable K-scores).
The source candidate might remain genus/observation-confounded.

## Inference design on exactly the same observations

Use the previously frozen five heldout geographic folds, whole-fold
accepted-name exclusion from *every* GIFT mainland donor source, known
observed island genus counts, fixed shrink=8 and the same potential
native-mainland source assignment.

Calculate per-island (all four kernel types):
- improvement from **no-genus global prior** to global genus predictions;
- incremental improvement from **global genus** to local mainland
  source genus predictions.

Separately within the same 1,007 All / 950 Direct support-and-donor-size
matched source assignment placebo cohorts, compare the **actual** nearby
source assignment against two outcome-blind, same-regime, different-block
candidate source assignments matching independent same-genus donor coverage
to +/-2 percentage points and complete donor name pool size within a
0.5–2x ratio.

Use equal-island averages, whole-spatial-block resampling, All/Direct
as nested scope sensitivities, and controls frozen without considering
colour/form score outcomes. Publish only group score aggregates and
kernel input digests. **Do not** persist island, donor source, species,
genus or spatial block identifiers or all private per-name score rows.

## Predeclared post-hoc reporting order (not a new experiment)

1. Is the previous Kadd improvement exactly the average of colour and form
   proper-score improvements? (must pass algebraic identity)
2. Which marginal predicts better: colour or form, and does either
   individually retain a geographical source-placebo gain under strict
   opportunity matching? (descriptive, correlated comparisons)
3. Does a **joint product** kernel on the same matched reference groups
   retain positive source skill? If it does not, the old additive score
   is a marginal signal, **not demonstrated transferable joint flower
   architecture**. If it does, report it cautiously as a joint-sensitive
   association; marginal changes alone can also generate that result.
4. Examine original scores alongside genuine full-category strictly
   proper Brier/log-loss. A joint similarity gain cannot override
   zero exact-log-loss improvement, dependence of taxa across islands,
   sparse source donor coverage, or unverified ancestral sources.

All these trait-specific checks are **post-hoc explorations of one
overlapping data set**, not four independent ecology tests. Never select
a kernel based on which happens to produce the smallest p. The original
2000/1000 km geographic packing references and N1b/N1c are supporting
descriptive stress tests, not equal evidence for source or selection.

## Running the public mathematical audit

    python -m fd.trait_kernels \
      --reference work/island_public/floral_space_20261106/all/profile_coordinates.csv \
      --output results/trait_kernel_audit

This writes **public only** rank/PSD/identity diagnostics and source
SHA256. No private names or mapped islands.

The actual N2 per-score extensions live in
`src/fd/n2_crossfit_macro.py` and
`src/fd/n2_source_assignment_placebo.py`.
Their outputs are respectively
`n2_trait_component_crossfit.csv` and
`n2_source_placebo_trait_kernels.csv`, produced **only** by an
authorized local replay with nonredistributable native-mainland taxon
sources. A successful public CI does not mean this private empirical
analysis has been run.

## Manuscript decision gate

Without the source-placebo **joint-product** check, the strongest
defensible statement is: *native mainland donor geography carries
some additional information about similarity among recorded floral
traits beyond known genus composition, although strict donor-support
matching greatly attenuates that advantage*. It is not yet proof that
specific **colour–form combinations** or pollinator-relevant floral
configurations transfer consistently from mainland source pools to
islands.
