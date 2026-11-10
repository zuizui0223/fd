"""Proper held-out probability scores on original discrete floral profiles.

Brier is strictly proper for full category distributions. Kernel quadratic
score is proper for a positive semidefinite K = 1 - original colour/form
state-set Jaccard D; it is not strictly proper when K is rank deficient.
These scores evaluate predictions of names NOT used for donor training.
Do not infer causal adaptation or pollination function.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

LABELS = ("global_prior", "global_genus", "local_genus")
LOSS_TYPES = ("logloss", "trait_distance_loss", "brier", "kernel")


def inspect_kernel(distance: np.ndarray) -> dict:
    """Check finite 0..1 reference D, PSD trait similarity, and kernel rank."""
    d = np.asarray(distance, dtype=float)
    if d.ndim != 2 or d.shape[0] < 2 or d.shape[0] != d.shape[1]:
        raise ValueError("Expected square original categorical dissimilarity")
    if not np.isfinite(d).all() or np.max(np.abs(d-d.T)) > 1e-10:
        raise ValueError("Invalid/asymmetric original dissimilarity")
    if np.max(np.abs(np.diag(d))) > 1e-10 or np.min(d) < -1e-10 or np.max(d) > 1+1e-10:
        raise ValueError("Original dissimilarity out of bounds")
    kernel = 1.0 - d
    eig = np.linalg.eigvalsh(kernel)
    if eig.min() < -1e-8:
        raise ValueError("K=1-D is NOT positive semidefinite: kernel score not proper")
    return {
        "n_categories": len(d),
        "minimum_eigenvalue": float(eig.min()),
        "rank_positive_1e8": int(np.count_nonzero(eig > 1e-8)),
        "kernel": kernel,
    }


def heldout_probability_scores(part: pd.DataFrame, model, island: str,
                               similarity: np.ndarray) -> dict[str, float]:
    """Island-averaged scores, grouping shared genus to avoid repeated kernels.

    The global and local genus distributions are conditional on KNOWING the
    island genus identities, not prospective lineage colonization predictions.
    K must have been PSD-validated against the reference D once per run.
    """
    if part.empty or part.accepted_species.duplicated().any():
        raise ValueError("Need unique held-out taxon names on one island")
    n = model.n_profiles
    k = np.asarray(similarity, dtype=float)
    if k.shape != (n, n) or not np.isfinite(k).all():
        raise ValueError("Similarity reference dimension mismatch")
    loss = np.zeros((3, 4), dtype=float)
    for genus, records in part.groupby("genus", sort=False):
        counts = records.profile_id.astype(int).value_counts()
        ids = counts.index.to_numpy(dtype=int)
        weights = counts.to_numpy(dtype=float)
        if (ids < 0).any() or (ids >= n).any():
            raise ValueError("Heldout category outside profile reference")
        prob = np.asarray(model.probabilities(island, genus)[[0, 2, 3]],
                          dtype=float)
        if prob.shape != (3, n) or not np.isfinite(prob).all() or (prob < 0).any():
            raise ValueError("Malformed donor probability vector")
        if not np.allclose(prob.sum(axis=1), 1, atol=1e-9):
            raise ValueError("Genus probabilities must sum to one")
        # Shape (3, n); each distinct genus incurs one matrix application.
        kp = prob @ k
        n_g = float(weights.sum())
        loss[:, 0] += -(np.log(np.maximum(prob[:, ids], 1e-15)) @ weights)
        loss[:, 1] += ((1.0 - kp[:, ids]) @ weights)
        loss[:, 2] += n_g * (1 + (prob*prob).sum(axis=1)) - 2 * (prob[:, ids] @ weights)
        loss[:, 3] += n_g * (1 + (prob*kp).sum(axis=1)) - 2 * (kp[:, ids] @ weights)
    loss /= len(part)
    return {f"{name}_{metric}": float(loss[i, j])
            for i, name in enumerate(LABELS)
            for j, metric in enumerate(LOSS_TYPES)}
