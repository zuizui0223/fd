"""Marginal versus joint floral set-state similarity kernels.

The original D is exactly the equally weighted mean of Jaccard distances on
COLOUR SETS and FORM SETS, so original K = 1-D is ADDITIVE in these traits.
The additive K-score cannot distinguish colour/form pairings when both
marginals are identical. A Schur-product Jaccard kernel can probe joint
colour/form association without multiplying recorded multistate descriptions
into fabricated profiles.

All kernels are exploratory ecological proxy scores, not pollinator fitness.
All profile identifiers and state sets here come from the PUBLIC 496-row
reference. This module never reads or emits private island/taxon data.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Mapping

import numpy as np

from .core import Profile, load_profiles
from .pipeline import file_sha256
from .proper_scoring import inspect_kernel

KERNEL_NAMES = ("colour", "form", "additive", "joint_product")


def _set_jaccard(a: frozenset[str], b: frozenset[str]) -> float:
    if not a or not b:
        raise ValueError("Reference trait-state set cannot be empty")
    return len(a & b) / len(a | b)


def phenotype_kernels(profiles: Mapping[int, Profile]) -> dict[str, np.ndarray]:
    """Return K_colour, K_form, original K_additive and K_colour * K_form.

    A multistate value is one unchanged recorded SET (not the Cartesian
    product of species-level phenotype combinations).
    """
    if len(profiles) < 2:
        raise ValueError("Need at least two distinct reference profiles")
    ids = sorted(profiles)
    n = len(ids)
    colour = np.eye(n, dtype=float)
    form = np.eye(n, dtype=float)
    for i, id_i in enumerate(ids):
        c_i, f_i = profiles[id_i]
        for j in range(i + 1, n):
            c_j, f_j = profiles[ids[j]]
            a = _set_jaccard(c_i, c_j)
            b = _set_jaccard(f_i, f_j)
            colour[i, j] = colour[j, i] = a
            form[i, j] = form[j, i] = b
    additive = (colour + form) * .5
    joint = colour * form  # Schur product of positive-semidefinite set kernels
    return {"colour": colour, "form": form,
            "additive": additive, "joint_product": joint}


def kernel_quality(profiles: Mapping[int, Profile],
                   *, verify_psd: bool = True) -> tuple[dict, dict]:
    """Validate additive D identity and strictly distinguishable joint K."""
    kernels = phenotype_kernels(profiles)
    info = {}
    for name, kernel in kernels.items():
        if not np.isfinite(kernel).all() or not np.allclose(
                np.diag(kernel), 1., atol=1e-12):
            raise ValueError(f"Invalid reference state-set kernel: {name}")
        if not np.allclose(kernel, kernel.T, atol=1e-12):
            raise ValueError(f"Asymmetric reference state-set kernel: {name}")
        if verify_psd:
            audit = inspect_kernel(1 - kernel)
            info[name] = {
                "rank_positive_1e8": audit["rank_positive_1e8"],
                "minimum_eigenvalue": audit["minimum_eigenvalue"],
                "n_reference_profiles": len(kernel),
            }
    return kernels, info


def quadratic_kernel_score(prob: np.ndarray, observed_id: int,
                           kernel: np.ndarray) -> float:
    """Proper quadratic kernel score for one categorical heldout observation."""
    p = np.asarray(prob, dtype=float)
    k = np.asarray(kernel, dtype=float)
    if p.ndim != 1 or k.shape != (len(p), len(p)):
        raise ValueError("Predicted probability and kernel dimensions differ")
    if not (0 <= observed_id < len(p)) or np.any(p < 0) or (
            not np.isfinite(p).all()) or not np.isclose(p.sum(), 1):
        raise ValueError("Invalid categorical prediction/target")
    return float(p @ k @ p - 2 * (k @ p)[observed_id] + k[observed_id, observed_id])


def audit_reference(reference: Path, output: Path) -> dict:
    """Public-only CLI run; log only kernel diagnostics, never taxon records."""
    reference, output = Path(reference), Path(output)
    profiles = load_profiles(reference)
    if len(profiles) != 496 or sorted(profiles) != list(range(496)):
        raise ValueError("Frozen original 496-profile reference changed")
    kernels, info = kernel_quality(profiles)
    additive_error = float(np.max(np.abs(
        kernels["additive"] - (kernels["colour"] + kernels["form"]) / 2
    )))
    if additive_error > 1e-12:
        raise AssertionError("Original additive colour/form kernel identity failed")
    # Kadd = 1 - 0.5*(Jaccard colour distance + Jaccard form distance).
    result = {
        "schema": "fd_public_original_floral_trait_kernel_audit_v1",
        "status": "public_reference_mathematical_identity_only",
        "reference_sha256": file_sha256(reference),
        "reference_profiles": len(profiles),
        "kernel_quality": info,
        "additive_formula": "Kadd=(K_colour+K_form)/2=1-original_D",
        "joint_formula": "Kjoint=K_colour * K_form (Hadamard/Schur)",
        "additive_max_identity_error": additive_error,
        "scope": "trait scoring algebra only; no island-level effect, no source prediction, no adaptation claim",
        "interaction_limit": "product kernel can detect joint colour/form association, but a joint K-score gain is not uniquely pure interaction after conditioning on changing marginals",
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "trait_kernel_audit.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"reference_profiles": result["reference_profiles"],
                      "kernels": result["kernel_quality"],
                      "additive_max_identity_error": additive_error}))
    return result


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    audit_reference(a.reference, a.output)


if __name__ == "__main__":
    main()
