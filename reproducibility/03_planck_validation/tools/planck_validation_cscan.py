"""Planck TT validation c-scan at fixed a2=8000 (see planck_validation.py).

Quantifies the exclusion bound on the lambda-coupling c: how small must c be
for GenFT Model A to remain consistent with Planck 2018 binned TT?

Usage: python tools/planck_validation_cscan.py --outdir outputs/planck_validation_cscan_20260704
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import camb
from camb.initialpower import SplinedInitialPower
from genftcamb.build_pk import SpectrumParams, build_primordial_spectrum
from genftcamb.io import load_genft_input
from planck_validation import (
    DATA_DIR, PLANCK_FILE, load_planck_binned, camb_params, get_lensed_dl,
    bin_theory, chi2, pk_deviation_stats,
)

C_FILES = {
    0.0: "model_a_kvac_dense_ell100_n12_a2_8000_c_0.npz",
    0.05: "model_a_kvac_dense_ell100_n12_a2_8000_c_0p05.npz",
    0.1: "model_a_kvac_dense_ell100_n12_a2_8000_c_0p1.npz",
    0.15: "model_a_kvac_dense_ell100_n12_a2_8000_c_0p15.npz",
    0.2: "model_a_kvac_dense_ell100_n12_a2_8000_c_0p2.npz",
    0.25: "model_a_kvac_dense_ell100_n12_a2_8000_c_0p25.npz",
    0.3: "model_a_kvac_dense_ell100_n12_a2_8000_c_0p3.npz",
    0.35: "model_a_kvac_dense_ell100_n12_a2_8000_c_0p35.npz",
    0.4: "model_a_kvac_dense_ell100_n12_a2_8000_c_0p4.npz",
    0.5: "model_a_kvac_dense_ell100_n12_a2_8000_c_0p5.npz",
    0.7: "model_a_kvac_dense_ell100_n12_a2_8000_c_0p7.npz",
    1.0: "model_a_kvac_dense_ell100_n12_a2_8000_c_1.npz",
    1.5: "model_a_kvac_dense_ell100_n12_a2_8000_c_1p5.npz",
    2.0: "model_a_kvac_dense_ell100_n12_a2_8000_c_2.npz",
    3.0: "model_a_kvac_dense_ell100_n12_a2_8000_c_3.npz",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", type=float, default=0.00118)
    ap.add_argument("--As", type=float, default=2.1e-9)
    ap.add_argument("--ns", type=float, default=0.965)
    ap.add_argument("--k0", type=float, default=0.05)
    ap.add_argument("--lmax", type=int, default=2600)
    ap.add_argument("--ell-min", type=float, default=30)
    ap.add_argument("--ell-max", type=float, default=2000)
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    stamp = time.strftime("%Y%m%d_%H%M%S")
    outdir = pathlib.Path(args.outdir or f"outputs/planck_validation_cscan_{stamp}")
    outdir.mkdir(parents=True, exist_ok=True)

    ell_b, dl_b, err_b, edges = load_planck_binned(PLANCK_FILE, args.ell_min, args.ell_max)
    nbins = ell_b.size

    p = camb_params(args.lmax)
    p.InitPower.set_params(As=args.As, ns=args.ns)
    dl_lcdm = get_lensed_dl(p, args.lmax)
    chi2_lcdm = chi2(dl_b, err_b, bin_theory(dl_lcdm, edges))
    print(f"[lcdm] chi2 = {chi2_lcdm:.1f} / {nbins} bins")

    sp = SpectrumParams(kappa=args.kappa, As=args.As, ns=args.ns,
                        k0=args.k0, mu_ref="median")
    rows = []
    for c, fname in sorted(C_FILES.items()):
        path = DATA_DIR / fname
        if not path.exists():
            print(f"[skip] {fname}")
            continue
        t0 = time.time()
        lambdas, mus = load_genft_input(path)
        k, pz, meta = build_primordial_spectrum(lambdas, mus, sp)
        pparams = camb_params(args.lmax)
        pparams.InitPower = SplinedInitialPower()
        pparams.InitPower.effective_ns_for_nonlinear = args.ns
        pparams.InitPower.set_scalar_table(k, pz)
        dl_g = get_lensed_dl(pparams, args.lmax)
        c2 = chi2(dl_b, err_b, bin_theory(dl_g, edges))
        lo, hi = int(args.ell_min), int(args.ell_max)
        frac = dl_g[lo:hi + 1] / dl_lcdm[lo:hi + 1] - 1.0
        rms = float(np.sqrt(np.mean(frac ** 2)))
        dev = pk_deviation_stats(lambdas, mus, sp)
        mu_rel_std = float(np.std(mus) / np.mean(mus))
        rows.append(dict(c=c, chi2_planck=c2, dchi2_vs_lcdm=c2 - chi2_lcdm,
                         rms_frac_vs_lcdm_theory=rms,
                         pk_dev_rms=dev["dev_rms"],
                         mu_rel_std=mu_rel_std,
                         runtime_s=round(time.time() - t0, 1)))
        print(f"[c={c:g}] chi2={c2:.1f} dchi2={c2 - chi2_lcdm:+.2f} "
              f"rms={rms:.3e} mu_rel_std={mu_rel_std:.3e}")

    with open(outdir / "planck_cscan_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(outdir / "summary.json", "w") as f:
        json.dump(dict(a2=8000, kappa=args.kappa, lmax=args.lmax,
                       nbins=nbins, chi2_lcdm=chi2_lcdm,
                       camb_version=camb.__version__), f, indent=2)
    print(f"[done] {outdir}")


if __name__ == "__main__":
    main()
