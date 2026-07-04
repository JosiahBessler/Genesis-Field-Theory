"""Verify the December 2025 MID-vs-UCP eigenmode comparison numbers.

Recomputes, from the stored mode-shape table, the r^2-weighted L2 norms
reported in GenesisFT_Eigenmode_Comparison_Report (2025-12-03) and checks
them against the published values:

    ||phi_mid||          = 3.595280
    ||phi_ucp||          = 1.223896
    ||phi_mid - phi_ucp||= 2.587374
    relative L2 error    = 0.7196587
    lambda_free          = 0.1506
    lambda_eff           = 0.7114

Usage:  python verify_eigenmode_numbers.py
Exits 0 if all values reproduce within tolerance, 1 otherwise.
"""
import sys
import numpy as np
import pandas as pd

EXPECTED = {
    "norm_mid": 3.595280,
    "norm_ucp": 1.223896,
    "norm_diff": 2.587374,
    "rel_l2": 0.7196587,
    "lambda_free": 0.15061,
    "lambda_eff": 0.71135,
}
TOL = 1e-4  # relative


def main() -> int:
    df = pd.read_csv("data/mid_vs_ucp_mode_shapes.csv")
    r = df["r_kpc"].values
    mid = df["phi_mid_norm"].values
    ucp = df["phi_ucp_norm"].values

    w = r ** 2  # 3-D radial volume measure
    def l2(v):
        return float(np.sqrt(np.trapz(w * v ** 2, r)))

    got = {
        "norm_mid": l2(mid),
        "norm_ucp": l2(ucp),
        "norm_diff": l2(mid - ucp),
    }
    got["rel_l2"] = got["norm_diff"] / got["norm_mid"]

    ev = pd.read_csv("data/optimized_mid_mode_eigenvalues.csv")
    got["lambda_free"] = float(ev["lambda_free"].iloc[0])
    got["lambda_eff"] = float(ev["lambda_eff"].iloc[0])

    ok = True
    for key, expect in EXPECTED.items():
        val = got[key]
        rel = abs(val - expect) / abs(expect)
        status = "OK  " if rel < TOL else "FAIL"
        if rel >= TOL:
            ok = False
        print(f"{status} {key:12s} expected {expect:<12.7g} got {val:<12.7g} (rel diff {rel:.2e})")

    print("\nRESULT:", "all values reproduce" if ok else "MISMATCH — investigate")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
