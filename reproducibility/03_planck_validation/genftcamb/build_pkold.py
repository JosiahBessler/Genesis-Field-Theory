"""Construct primordial curvature power spectra from GenesisFT outputs."""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np

from .io import warn_and_filter_positive

logger = logging.getLogger(__name__)


@dataclass
class SpectrumParams:
    """Parameters controlling the primordial spectrum construction."""

    kappa: float = 1.0
    As: float = 2.1e-9
    ns: float = 0.965
    k0: float = 0.05
    mu_ref: str | float = "median"


def _resolve_mu_ref(mus: np.ndarray, mu_ref: str | float) -> float:
    if isinstance(mu_ref, str):
        key = mu_ref.lower()
        if key == "median":
            return float(np.median(mus))
        if key == "mean":
            return float(np.mean(mus))
        if key == "max":
            return float(np.max(mus))
        try:
            return float(mu_ref)
        except ValueError as exc:
            raise ValueError("mu_ref must be 'median', 'mean', 'max', or a numeric value") from exc
    return float(mu_ref)


def _average_duplicates(k: np.ndarray, pk: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Average pk for identical k values (after sorting)."""
    unique_k, inverse_indices = np.unique(k, return_inverse=True)
    accum = np.zeros_like(unique_k, dtype=float)
    counts = np.zeros_like(unique_k, dtype=int)
    for idx, inv in enumerate(inverse_indices):
        accum[inv] += pk[idx]
        counts[inv] += 1
    # counts should never be zero, but guard anyway:
    with np.errstate(divide="ignore", invalid="ignore"):
        pk_avg = np.divide(accum, counts, out=np.zeros_like(accum), where=counts > 0)
    return unique_k, pk_avg


def build_primordial_spectrum(
    lambdas: np.ndarray, mus: np.ndarray, params: SpectrumParams
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    """Build a CAMB-ready primordial power spectrum (k, Pzeta).

    Steps:
      - filter nonpositive lambda/mu
      - map k = kappa*sqrt(lambda)
      - compute Pzeta from mu weights + global tilt
      - sort, average duplicate k, enforce strictly increasing
      - clamp k to [1e-5, 50] for CAMB stability
      - clip pk to [1e-30, 1e-6] to prevent numerical blowups
    """

    lambdas = np.asarray(lambdas, dtype=float)
    mus = np.asarray(mus, dtype=float)
    
    if lambdas.shape != mus.shape:
        raise ValueError(f"Input arrays must have the same shape before filtering: lambdas {lambdas.shape}, mus {mus.shape}")
    
    mask = np.isfinite(lambdas) & np.isfinite(mus) & (lambdas > 0) & (mus > 0)
    removed = int(np.size(lambdas) - np.count_nonzero(mask))
    if removed > 0:
        logger.warning("Removed %d nonpositive/nonfinite paired (lambda, mu) entries.", removed)
    
    lambdas = lambdas[mask]
    mus = mus[mask]
    
    if lambdas.size == 0:
        raise ValueError("No valid (lambda, mu) pairs remain after filtering.")
    

    k = params.kappa * np.sqrt(lambdas)

    mu_ref_value = _resolve_mu_ref(mus, params.mu_ref)
    if not np.isfinite(mu_ref_value) or mu_ref_value <= 0:
        raise ValueError("mu_ref must resolve to a positive finite value.")
    if not np.isfinite(params.k0) or params.k0 <= 0:
        raise ValueError("k0 must be a positive finite value.")
    if not np.isfinite(params.As) or params.As <= 0:
        raise ValueError("As must be a positive finite value.")

    pk = params.As * (mus / mu_ref_value) ** 2 * (k / params.k0) ** (params.ns - 1)

    # Sort by k
    order = np.argsort(k)
    k_sorted = k[order]
    pk_sorted = pk[order]

    # Drop any non-finite points before uniqueness handling
    finite = np.isfinite(k_sorted) & np.isfinite(pk_sorted) & (pk_sorted > 0)
    k_sorted = k_sorted[finite]
    pk_sorted = pk_sorted[finite]
    if k_sorted.size < 50:
        raise ValueError("Too few finite positive spectrum points after filtering; check inputs or parameters.")

    # Average duplicates
    k_unique, pk_unique = _average_duplicates(k_sorted, pk_sorted)

    # Enforce strictly increasing k
    for i in range(1, len(k_unique)):
        if k_unique[i] <= k_unique[i - 1]:
            k_unique[i] = k_unique[i - 1] * (1.0 + 1e-12)

    # Clamp to a safe CAMB band
    clip_mask = (k_unique >= 1e-5) & (k_unique <= 50.0)
    if np.count_nonzero(clip_mask) < 50:
        raise ValueError("Too few k-points remain after clamping to [1e-5, 50]; adjust kappa or input scaling.")
    k_clamped = k_unique[clip_mask]
    pk_clamped = pk_unique[clip_mask]

    # Clip pk for numerical stability
    pk_clipped = np.clip(pk_clamped, 1e-30, 1e-6)

    meta: Dict[str, float] = {
        "mu_ref_value": float(mu_ref_value),
        "min_k": float(np.min(k_clamped)),
        "max_k": float(np.max(k_clamped)),
        "num_points": int(k_clamped.size),
    }

    return k_clamped, pk_clipped, meta
