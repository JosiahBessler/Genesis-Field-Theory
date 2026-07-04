"""CAMB interface helpers (with version-compatible custom primordial spectrum support)."""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Callable, Tuple

import numpy as np
import camb

logger = logging.getLogger(__name__)


@dataclass
class Cosmology:
    """Cosmological parameters used by CAMB."""
    H0: float = 67.4
    ombh2: float = 0.0224
    omch2: float = 0.12
    tau: float = 0.054


def _validate_table(k: np.ndarray, pk: np.ndarray) -> None:
    if k.ndim != 1 or pk.ndim != 1:
        raise ValueError("k and pk must be one-dimensional arrays.")
    if k.size != pk.size:
        raise ValueError("k and pk must have the same length.")
    if not np.all(np.isfinite(k)) or not np.all(np.isfinite(pk)):
        raise ValueError("k and pk must be finite.")
    if not np.all(pk > 0):
        raise ValueError("pk values must be strictly positive.")
    if not np.all(np.diff(k) > 0):
        raise ValueError("k must be strictly increasing.")


def _base_params(cosmo: Cosmology, lmax: int) -> camb.model.CAMBparams:
    params = camb.CAMBparams()
    params.set_cosmology(H0=cosmo.H0, ombh2=cosmo.ombh2, omch2=cosmo.omch2, tau=cosmo.tau)
    params.NonLinear = camb.model.NonLinear_none
    params.Want_CMB_lensing = False
    params.set_for_lmax(lmax, lens_potential_accuracy=0)
    return params


def _extract_cls(cl: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # CAMB returns columns [TT, EE, BB, TE] for "total" raw_cl
    ell = np.arange(cl.shape[0])
    tt = cl[:, 0]
    ee = cl[:, 1]
    bb = cl[:, 2]
    te = cl[:, 3]
    return ell, tt, ee, te, bb




def run_camb_with_table(
    cosmo: Cosmology, lmax: int, k: np.ndarray, pk: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Run CAMB using a custom scalar power spectrum table (robust across CAMB builds)."""
    _validate_table(k, pk)
    logger.info("Running CAMB with custom scalar table (%d points, lmax=%d)", k.size, lmax)

    params = _base_params(cosmo, lmax)

    # Force a splined primordial spectrum object that supports set_scalar_table
    from camb.initialpower import SplinedInitialPower
    params.InitPower = SplinedInitialPower()
    params.InitPower.set_scalar_table(k, pk)
    logger.info("Using camb.initialpower.SplinedInitialPower.set_scalar_table")

    results = camb.get_results(params)
    powers = results.get_cmb_power_spectra(lmax=lmax, spectra=["total"], raw_cl=True)
    cl = powers["total"]
    return _extract_cls(cl)





def run_camb_baseline(
    cosmo: Cosmology, lmax: int, As: float, ns: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Run CAMB with a baseline power-law spectrum."""
    logger.info("Running CAMB baseline (As=%.3e, ns=%.3f, lmax=%d)", As, ns, lmax)

    params = _base_params(cosmo, lmax)
    params.InitPower.set_params(As=As, ns=ns)

    results = camb.get_results(params)
    powers = results.get_cmb_power_spectra(lmax=lmax, spectra=["total"], raw_cl=True)
    cl = powers["total"]
    return _extract_cls(cl)
