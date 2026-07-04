"""
analyze_ucp_relations.py

Runs several diagnostic tests on the GenFT / UCP halo fits:

1. RAR (Radial Acceleration Relation) scatter using:
   - observed accelerations g_obs
   - baryonic accelerations g_bar
   - model accelerations g_tot_model for each halo model

2. Per-galaxy residual statistics for each halo profile
   (RMS of V_obs - V_model).

3. Simple scaling summary (V_flat proxy and r_scale per model).

Outputs:
  - Results/rar_points.csv
  - Results/residuals_summary.csv
  - Results/scaling_summary.csv

This script expects to live in the same folder as:
  - halo_dm_baryons_results.csv  (preferred) or halo_fit_results.csv
  - halo_models.py
  - Rotmod_LTG/*.dat  (SPARC rotmod files)
"""

import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from halo_models import get_default_models


# ------------------------------------------------------------
# Paths and I/O helpers
# ------------------------------------------------------------

BASE_DIR = os.path.dirname(__file__)
RESULTS_DIR = os.path.join(BASE_DIR, "Results")
ROTMOD_DIR = os.path.join(BASE_DIR, "Rotmod_LTG")

os.makedirs(RESULTS_DIR, exist_ok=True)

RESULT_CANDIDATES = [
    "halo_dm_baryons_results.csv",
    "halo_fit_results.csv",
]


@dataclass
class RotmodGalaxy:
    name: str
    R_kpc: np.ndarray
    Vobs: np.ndarray
    eVobs: np.ndarray
    Vgas: np.ndarray
    Vdisk: np.ndarray
    Vbul: np.ndarray


def find_results_file() -> str:
    for fname in RESULT_CANDIDATES:
        path = os.path.join(BASE_DIR, fname)
        if os.path.exists(path):
            return path
    raise FileNotFoundError(
        f"Could not find any of {RESULT_CANDIDATES} in {BASE_DIR}"
    )


def load_halo_results() -> pd.DataFrame:
    path = find_results_file()
    print(f"Loading halo fit results from: {path}")
    df = pd.read_csv(path)
    required = {"galaxy", "model", "r_scale", "V0", "chi2", "dof"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Results file is missing required columns: {missing}")
    return df


def parse_rotmod_file(path: str, name: str) -> RotmodGalaxy:
    data: List[List[float]] = []
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#") or line.startswith(";"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                vals = [float(x) for x in parts[:6]]
            except ValueError:
                continue
            data.append(vals)

    if not data:
        raise ValueError(f"No numeric data in {path}")

    arr = np.array(data, dtype=float)
    R = arr[:, 0]
    Vobs = arr[:, 1]
    eVobs = arr[:, 2]
    Vgas = arr[:, 3]
    Vdisk = arr[:, 4]
    Vbul = arr[:, 5]

    return RotmodGalaxy(
        name=name,
        R_kpc=R,
        Vobs=Vobs,
        eVobs=eVobs,
        Vgas=Vgas,
        Vdisk=Vdisk,
        Vbul=Vbul,
    )


def load_rotmod_galaxies() -> Dict[str, RotmodGalaxy]:
    if not os.path.isdir(ROTMOD_DIR):
        raise FileNotFoundError(f"Rotmod directory not found: {ROTMOD_DIR}")

    galaxies: Dict[str, RotmodGalaxy] = {}
    for fname in os.listdir(ROTMOD_DIR):
        if not fname.lower().endswith("_rotmod.dat"):
            continue
        base = fname[:-len("_rotmod.dat")]  # strip suffix
        full = os.path.join(ROTMOD_DIR, fname)
        try:
            gal = parse_rotmod_file(full, base)
        except Exception as e:
            print(f"[WARN] Skipping {fname}: {e}")
            continue
        galaxies[base] = gal

    print(f"Loaded {len(galaxies)} rotmod galaxies from {ROTMOD_DIR}")
    return galaxies


# ------------------------------------------------------------
# Physics helpers
# ------------------------------------------------------------

def compute_accelerations(gal: RotmodGalaxy) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute observed and baryonic accelerations (up to a common unit factor).

    g_obs ~ Vobs^2 / R
    g_bar ~ Vbar^2 / R, where Vbar^2 = Vgas^2 + Vdisk^2 + Vbul^2

    We keep units in (km/s)^2 / kpc. For RAR shape and scatter,
    the overall constant factor does not matter.
    """
    R = gal.R_kpc
    Vobs = gal.Vobs
    Vbar = np.sqrt(
        np.maximum(gal.Vgas, 0.0) ** 2
        + np.maximum(gal.Vdisk, 0.0) ** 2
        + np.maximum(gal.Vbul, 0.0) ** 2
    )

    # Avoid division by zero
    mask = (R > 0) & np.isfinite(R) & np.isfinite(Vobs) & (Vobs >= 0)
    mask &= np.isfinite(Vbar)

    R = R[mask]
    Vobs = Vobs[mask]
    Vbar = Vbar[mask]

    g_obs = (Vobs ** 2) / R
    g_bar = (Vbar ** 2) / R

    return g_obs, g_bar


def compute_model_velocity(
    halo_vfunc,
    r_scale: float,
    V0: float,
    R_kpc: np.ndarray,
    Vgas: np.ndarray,
    Vdisk: np.ndarray,
    Vbul: np.ndarray,
) -> np.ndarray:
    """
    Total model velocity:
        V_tot^2 = V_bar^2 + V_halo^2
    using the given halo velocity function and parameters.
    """
    Vbar = np.sqrt(
        np.maximum(Vgas, 0.0) ** 2
        + np.maximum(Vdisk, 0.0) ** 2
        + np.maximum(Vbul, 0.0) ** 2
    )
    # halo_vfunc expects (r, params)
    Vhalo = halo_vfunc(R_kpc, (r_scale, V0))
    Vmodel = np.sqrt(Vbar ** 2 + np.maximum(Vhalo, 0.0) ** 2)
    return Vmodel


def estimate_v_flat(Vobs: np.ndarray, R_kpc: np.ndarray, frac_outer: float = 0.3) -> float:
    """
    Crude V_flat proxy: median Vobs over the outer fraction of the radius range.
    """
    if len(R_kpc) == 0:
        return np.nan
    # sort by radius just in case
    idx = np.argsort(R_kpc)
    R_sorted = R_kpc[idx]
    V_sorted = Vobs[idx]
    n = len(R_sorted)
    n_outer = max(3, int(np.ceil(frac_outer * n)))
    outer_vals = V_sorted[-n_outer:]
    return float(np.median(outer_vals))


# ------------------------------------------------------------
# Main analysis
# ------------------------------------------------------------

def main():
    # Load halo fit results
    df = load_halo_results()

    # Load rotmod data
    rotmod = load_rotmod_galaxies()

    # Set up halo models dict by name
    models = get_default_models(include_einasto=True)
    model_map = {m.name: m for m in models}

    # -------------------------
    # 1. RAR point cloud
    # -------------------------
    rar_rows = []

    # -------------------------
    # 2. Residuals per galaxy/model
    # -------------------------
    resid_rows = []

    # -------------------------
    # 3. Simple scaling per galaxy/model
    # -------------------------
    scale_rows = []

    # loop over galaxies that have rotmod and halo fits
    galaxies_with_data = sorted(set(df["galaxy"]) & set(rotmod.keys()))
    print(f"Galaxies with both halo fits and rotmod data: {len(galaxies_with_data)}")

    for gname in galaxies_with_data:
        gal = rotmod[gname]

        # observed vs baryonic accelerations (for RAR)
        g_obs, g_bar = compute_accelerations(gal)

        # store RAR points for this galaxy using direct SPARC data
        for go, gb in zip(g_obs, g_bar):
            rar_rows.append(
                {
                    "galaxy": gname,
                    "g_obs_data": go,
                    "g_bar": gb,
                }
            )

        # estimate V_flat proxy
        V_flat = estimate_v_flat(gal.Vobs, gal.R_kpc)

        # now handle each model separately
        sub = df[df["galaxy"] == gname]
        for _, row in sub.iterrows():
            model_name = row["model"]
            if model_name not in model_map:
                # skip models we don't know (e.g. alternate Einasto alphas)
                continue
            halo_model = model_map[model_name]

            r_scale = float(row["r_scale"])
            V0 = float(row["V0"])

            # compute model velocity curve
            V_model = compute_model_velocity(
                halo_model.vfunc,
                r_scale,
                V0,
                gal.R_kpc,
                gal.Vgas,
                gal.Vdisk,
                gal.Vbul,
            )

            # residual RMS (only where errors and velocities are valid)
            mask = (
                np.isfinite(gal.Vobs)
                & np.isfinite(gal.eVobs)
                & (gal.eVobs > 0)
                & np.isfinite(V_model)
            )
            if not np.any(mask):
                continue
            resid = gal.Vobs[mask] - V_model[mask]
            rms = float(np.sqrt(np.mean(resid ** 2)))

            resid_rows.append(
                {
                    "galaxy": gname,
                    "model": model_name,
                    "rms_resid_kms": rms,
                    "chi2": float(row["chi2"]),
                    "dof": int(row["dof"]),
                    "chi2_dof": float(row.get("chi2_dof", row["chi2"] / max(row["dof"], 1))),
                }
            )

            # add scaling info: r_scale, V0, V_flat
            scale_rows.append(
                {
                    "galaxy": gname,
                    "model": model_name,
                    "r_scale_kpc": r_scale,
                    "V0_kms": V0,
                    "V_flat_kms": V_flat,
                }
            )

            # For RAR: add model-based g_obs_model curve as well
            # g_tot_model ~ V_model^2 / R
            R = gal.R_kpc
            mask2 = (R > 0) & np.isfinite(R) & np.isfinite(V_model)
            Rm = R[mask2]
            Vm = V_model[mask2]
            Vbar = np.sqrt(
                np.maximum(gal.Vgas[mask2], 0.0) ** 2
                + np.maximum(gal.Vdisk[mask2], 0.0) ** 2
                + np.maximum(gal.Vbul[mask2], 0.0) ** 2
            )
            g_model_tot = (Vm ** 2) / Rm
            g_model_bar = (Vbar ** 2) / Rm

            for go, gb in zip(g_model_tot, g_model_bar):
                rar_rows.append(
                    {
                        "galaxy": gname,
                        "g_obs_data": np.nan,  # indicates model curve
                        "g_bar": gb,
                        "g_tot_model": go,
                        "model": model_name,
                    }
                )

    # -------------------------
    # Save outputs
    # -------------------------
    rar_df = pd.DataFrame(rar_rows)
    rar_path = os.path.join(RESULTS_DIR, "rar_points.csv")
    rar_df.to_csv(rar_path, index=False)
    print(f"Saved RAR point cloud to: {rar_path}")

    resid_df = pd.DataFrame(resid_rows)
    resid_path = os.path.join(RESULTS_DIR, "residuals_summary.csv")
    resid_df.to_csv(resid_path, index=False)
    print(f"Saved residual stats to: {resid_path}")

    scale_df = pd.DataFrame(scale_rows)
    scale_path = os.path.join(RESULTS_DIR, "scaling_summary.csv")
    scale_df.to_csv(scale_path, index=False)
    print(f"Saved scaling summary to: {scale_path}")

    # Some quick printed summaries
    if not resid_df.empty:
        print("\nMedian chi2/dof by model:")
        print(resid_df.groupby("model")["chi2_dof"].median().sort_values())

        print("\nMedian RMS residual (km/s) by model:")
        print(resid_df.groupby("model")["rms_resid_kms"].median().sort_values())

    if not scale_df.empty:
        print("\nMedian r_scale (kpc) by model:")
        print(scale_df.groupby("model")["r_scale_kpc"].median().sort_values())


if __name__ == "__main__":
    main()
