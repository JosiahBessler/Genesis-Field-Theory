"""I/O utilities for GenesisFT spectral data and CAMB outputs.

This module provides helpers to read GenesisFT eigenvalues/weights from NPZ or CSV
formats and to write spectrum tables and CAMB outputs in both CSV and NPZ forms.
"""
from __future__ import annotations

import logging
import pathlib
import warnings
from typing import Iterable, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class InputFormatError(Exception):
    """Raised when an unsupported input format is provided."""


class DataValidationError(Exception):
    """Raised when input data cannot be validated for processing."""


def load_genft_input(path: pathlib.Path) -> Tuple[np.ndarray, np.ndarray]:
    """Load GenesisFT eigenvalues and weights from NPZ or CSV.

    Parameters
    ----------
    path:
        Path to the input file. Supported formats: `.npz` containing arrays
        ``lambdas`` and ``mus``, or `.csv` with columns ``lambda`` and ``mu``.

    Returns
    -------
    lambdas, mus : ndarray
        One-dimensional arrays of eigenvalues and weights.
    """

    path = pathlib.Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    suffix = path.suffix.lower()
    if suffix == ".npz":
        data = np.load(path)
        if "lambdas" not in data or "mus" not in data:
            raise DataValidationError("NPZ file must contain 'lambdas' and 'mus'.")
        lambdas = np.array(data["lambdas"], dtype=float)
        mus = np.array(data["mus"], dtype=float)
    elif suffix == ".csv":
        df = pd.read_csv(path)
        if not {"lambda", "mu"}.issubset(df.columns):
            raise DataValidationError("CSV input must have columns 'lambda' and 'mu'.")
        lambdas = df["lambda"].to_numpy(dtype=float)
        mus = df["mu"].to_numpy(dtype=float)
    else:
        raise InputFormatError(f"Unsupported input format: {suffix}")

    return lambdas, mus


def warn_and_filter_positive(values: np.ndarray, label: str) -> np.ndarray:
    """Filter nonpositive entries, issuing a warning when values are removed."""
    mask = values > 0
    removed = np.count_nonzero(~mask)
    if removed > 0:
        warnings.warn(f"Removed {removed} nonpositive {label} entries.")
        logger.warning("Removed %d nonpositive %s entries", removed, label)
    return values[mask]


def ensure_output_dir(outdir: pathlib.Path) -> None:
    """Create an output directory if it does not already exist."""
    outdir.mkdir(parents=True, exist_ok=True)


def save_pk(outdir: pathlib.Path, k: np.ndarray, pzeta: np.ndarray) -> None:
    """Save primordial spectrum table to CSV and NPZ."""
    ensure_output_dir(outdir)
    df = pd.DataFrame({"k": k, "Pzeta": pzeta})
    df.to_csv(outdir / "pk.csv", index=False)
    np.savez(outdir / "pk.npz", k=k, Pzeta=pzeta)


def save_cls(outdir: pathlib.Path, prefix: str, ell: np.ndarray, tt: np.ndarray, ee: np.ndarray, te: np.ndarray, bb: np.ndarray) -> None:
    """Save CAMB CMB power spectra."""
    ensure_output_dir(outdir)
    df = pd.DataFrame({"ell": ell, "TT": tt, "EE": ee, "TE": te, "BB": bb})
    df.to_csv(outdir / f"cls_{prefix}.csv", index=False)
    np.savez(outdir / f"cls_{prefix}.npz", ell=ell, TT=tt, EE=ee, TE=te, BB=bb)


def save_report(outdir: pathlib.Path, content: str) -> None:
    """Write a report markdown file."""
    ensure_output_dir(outdir)
    (outdir / "report.md").write_text(content)


def describe_parameters(params: Iterable[Tuple[str, object]]) -> str:
    """Format parameter key/value pairs for reporting."""
    lines = []
    for key, value in params:
        lines.append(f"- **{key}**: {value}")
    return "\n".join(lines)
