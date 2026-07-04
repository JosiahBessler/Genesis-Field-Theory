"""Plotting utilities for GenesisFT CAMB pipeline."""
from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import numpy as np

from .io import ensure_output_dir


def plot_pzeta(outdir: pathlib.Path, k: np.ndarray, pzeta: np.ndarray) -> None:
    """Plot primordial spectrum in log-log space."""
    ensure_output_dir(outdir / "plots")
    fig, ax = plt.subplots()
    ax.loglog(k, pzeta)
    ax.set_xlabel(r"$k$ [Mpc$^{-1}$]")
    ax.set_ylabel(r"$P_\zeta(k)$")
    ax.set_title("Primordial curvature power spectrum")
    fig.tight_layout()
    fig.savefig(outdir / "plots" / "pzeta.png")
    plt.close(fig)


def plot_tt(outdir: pathlib.Path, ell: np.ndarray, tt_genft: np.ndarray, tt_lcdm: np.ndarray) -> None:
    """Plot TT spectra overlay."""
    ensure_output_dir(outdir / "plots")
    fig, ax = plt.subplots()
    ax.plot(ell, tt_genft, label="GenesisFT")
    ax.plot(ell, tt_lcdm, label=r"ΛCDM baseline", linestyle="--")
    ax.set_xlabel(r"$\ell$")
    ax.set_ylabel(r"$C_\ell^{TT}$")
    ax.set_title("CMB TT power spectra")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "plots" / "tt_overlay.png")
    plt.close(fig)


def plot_tt_residual(outdir: pathlib.Path, ell: np.ndarray, residual: np.ndarray, ell_min: int, ell_max: int) -> None:
    """Plot fractional TT residuals."""
    ensure_output_dir(outdir / "plots")
    fig, ax = plt.subplots()
    ax.plot(ell, residual)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xlim(ell_min, ell_max)
    ax.set_xlabel(r"$\ell$")
    ax.set_ylabel("Fractional residual")
    ax.set_title("TT fractional residuals")
    fig.tight_layout()
    fig.savefig(outdir / "plots" / "tt_residuals.png")
    plt.close(fig)
