"""Metrics for comparing GenesisFT-derived spectra with a baseline."""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass


@dataclass
class ScoreConfig:
    ell_min: int = 30
    ell_max: int = 2000
    frac_error: float = 0.05


def fractional_residual(tt_genft: np.ndarray, tt_lcdm: np.ndarray) -> np.ndarray:
    """Compute (tt_genft - tt_lcdm) / tt_lcdm with safe handling."""
    tt_genft = np.asarray(tt_genft)
    tt_lcdm = np.asarray(tt_lcdm)
    with np.errstate(divide="ignore", invalid="ignore"):
        residual = np.divide(tt_genft - tt_lcdm, tt_lcdm, out=np.zeros_like(tt_genft), where=tt_lcdm != 0)
    return residual


def score_tt(ell: np.ndarray, tt_genft: np.ndarray, tt_lcdm: np.ndarray, config: ScoreConfig) -> dict[str, float]:
    """Compute RMS fractional residual and a simple chi-square-like metric."""
    ell = np.asarray(ell)
    tt_genft = np.asarray(tt_genft)
    tt_lcdm = np.asarray(tt_lcdm)

    mask = (ell >= config.ell_min) & (ell <= config.ell_max)
    if not np.any(mask):
        raise ValueError("No ell values fall within the scoring range.")

    res = fractional_residual(tt_genft, tt_lcdm)
    res_range = res[mask]

    rms = float(np.sqrt(np.mean(res_range ** 2)))
    sigma = config.frac_error * tt_lcdm[mask]
    chi2 = float(np.sum(((tt_genft[mask] - tt_lcdm[mask]) / sigma) ** 2))

    return {"rms_fractional_residual": rms, "chi2_like": chi2}
