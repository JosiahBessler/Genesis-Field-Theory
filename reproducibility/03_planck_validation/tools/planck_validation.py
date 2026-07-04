"""Planck 2018 TT validation scan for GenesisFT Model A.

Scores GenesisFT-derived TT spectra (lensed, for fair comparison with the
observed sky) against the Planck 2018 binned TT power spectrum
(COM_PowerSpect_CMB-TT-binned_R3.01.txt, ell >= 30), as a function of the
kernel-vacuum stabilization parameter a2.

For each a2 input file the script:
  1. builds the primordial spectrum Pzeta(k) with the standard pipeline
     mapping k = kappa*sqrt(lambda), Pzeta = As*(mu/mu_ref)^2*(k/k0)^(ns-1);
  2. runs CAMB with lensing (lens_potential_accuracy=1) for the GenesisFT
     table and once for the LCDM power-law baseline (same cosmology);
  3. bin-averages the theory D_ell over the Planck bins and computes a
     diagonal chi^2 against the binned data (symmetrized errors);
  4. records the fractional TT residual of GenesisFT vs the LCDM theory
     curve, and the primordial-spectrum deviation amplitude
     (rms of the detrended (mu/mu_ref)^2 ratio in the observable window).

Caveats (stated in the report): diagonal chi^2 with symmetrized errors and
no bin-bin covariance, no foreground/nuisance marginalization, fixed
cosmology (H0=67.4, ombh2=0.0224, omch2=0.12, tau=0.054), TT only.
This is a physics-level discriminator, not a publication-grade likelihood.

Usage (from the pipeline root):
  python tools/planck_validation.py --kappa 0.00118 --lmax 2600
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import camb
from genftcamb.build_pk import SpectrumParams, build_primordial_spectrum
from genftcamb.io import load_genft_input

DATA_DIR = pathlib.Path(__file__).resolve().parents[1] / "data"
PLANCK_FILE = DATA_DIR / "COM_PowerSpect_CMB-TT-binned_R3.01.txt"

A2_FILES = {
    0.01: "model_a_kvac_dense_ell100_n12_a2_0p01.npz",
    0.1: "model_a_kvac_dense_ell100_n12_a2_0p1.npz",
    1: "model_a_kvac_dense_ell100_n12_a2_1.npz",
    2: "model_a_kvac_dense_ell100_n12_a2_2.npz",
    3: "model_a_kvac_dense_ell100_n12_a2_3.npz",
    4: "model_a_kvac_dense_ell100_n12_a2_4.npz",
    5: "model_a_kvac_dense_ell100_n12_a2_5.npz",
    6: "model_a_kvac_dense_ell100_n12_a2_6.npz",
    7: "model_a_kvac_dense_ell100_n12_a2_7.npz",
    8: "model_a_kvac_dense_ell100_n12_a2_8.npz",
    9: "model_a_kvac_dense_ell100_n12_a2_9.npz",
    12: "model_a_kvac_dense_ell100_n12_a2_12.npz",
    20: "model_a_kvac_dense_ell100_n12_a2_20.npz",
    30: "model_a_kvac_dense_ell100_n12_a2_30.npz",
    50: "model_a_kvac_dense_ell100_n12_a2_50.npz",
    80: "model_a_kvac_dense_ell100_n12_a2_80.npz",
    120: "model_a_kvac_dense_ell100_n12_a2_120.npz",
    200: "model_a_kvac_dense_ell100_n12_a2_200.npz",
    400: "model_a_kvac_dense_ell100_n12_a2_400.npz",
    600: "model_a_kvac_dense_ell100_n12_a2_600.npz",
    800: "model_a_kvac_dense_ell100_n12_a2_800.npz",
    1000: "model_a_kvac_dense_ell100_n12_a2_1000.npz",
    1500: "model_a_kvac_dense_ell100_n12_a2_1500.npz",
    2000: "model_a_kvac_dense_ell100_n12_a2_2000.npz",
    3000: "model_a_kvac_dense_ell100_n12_a2_3000.npz",
    5000: "model_a_kvac_dense_ell100_n12_a2_5000.npz",
    8000: "model_a_kvac_dense_ell100_n12_a2_8000.npz",
    12000: "model_a_kvac_dense_ell100_n12_a2_12000.npz",
}
BEST_HI_RES = ("best_ell120_n14_a2_8000", "BEST_dense_ell120_n14_a2_8000_c_0.npz")


def load_planck_binned(path: pathlib.Path, ell_min: float, ell_max: float):
    raw = np.loadtxt(path)
    ell, dl, err_lo, err_hi = raw[:, 0], raw[:, 1], raw[:, 2], raw[:, 3]
    err = 0.5 * (np.abs(err_lo) + np.abs(err_hi))
    keep = (ell >= ell_min) & (ell <= ell_max)
    ell, dl, err = ell[keep], dl[keep], err[keep]
    # reconstruct bin edges from midpoints between consecutive centers
    edges = np.empty(ell.size + 1)
    edges[1:-1] = 0.5 * (ell[1:] + ell[:-1])
    edges[0] = ell[0] - (edges[1] - ell[0])
    edges[-1] = ell[-1] + (ell[-1] - edges[-2])
    return ell, dl, err, edges


def camb_params(lmax: int) -> camb.model.CAMBparams:
    p = camb.CAMBparams()
    p.set_cosmology(H0=67.4, ombh2=0.0224, omch2=0.12, tau=0.054)
    p.NonLinear = camb.model.NonLinear_both
    p.set_for_lmax(lmax + 300, lens_potential_accuracy=1)
    return p


def get_lensed_dl(params: camb.model.CAMBparams, lmax: int) -> np.ndarray:
    results = camb.get_results(params)
    powers = results.get_cmb_power_spectra(
        params, lmax=lmax, spectra=["total"], CMB_unit="muK", raw_cl=False
    )
    return powers["total"][:, 0]  # D_ell TT in muK^2, index = ell


def bin_theory(dl_theory: np.ndarray, edges: np.ndarray) -> np.ndarray:
    ells = np.arange(dl_theory.size)
    out = np.empty(edges.size - 1)
    for i in range(edges.size - 1):
        m = (ells >= edges[i]) & (ells < edges[i + 1])
        out[i] = dl_theory[m].mean()
    return out


def chi2(dl_data, err, dl_model_binned) -> float:
    return float(np.sum(((dl_data - dl_model_binned) / err) ** 2))


def pk_deviation_stats(lambdas, mus, params: SpectrumParams,
                       k_lo=0.005, k_hi=0.2, window=101):
    """rms of the detrended Pzeta ratio (mu/mu_ref)^2 in the observable window."""
    k, pz, meta = build_primordial_spectrum(lambdas, mus, params)
    smooth = params.As * (k / params.k0) ** (params.ns - 1)
    ratio = pz / smooth  # equals (mu/mu_ref)^2 after pipeline filtering
    m = (k >= k_lo) & (k <= k_hi)
    r = ratio[m]
    if r.size < window + 1:
        return dict(dev_rms=np.nan, dev_detrended_rms=np.nan, n=int(r.size))
    # moving-average trend in index space (k is dense and near-uniformly packed)
    kernel = np.ones(window) / window
    trend = np.convolve(r, kernel, mode="same")
    resid = (r - trend) / trend
    return dict(
        dev_rms=float(np.sqrt(np.mean((r / np.median(r) - 1.0) ** 2))),
        dev_detrended_rms=float(np.sqrt(np.mean(resid ** 2))),
        n=int(r.size),
    )


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
    outdir = pathlib.Path(args.outdir or f"outputs/planck_validation_{stamp}")
    outdir.mkdir(parents=True, exist_ok=True)

    ell_b, dl_b, err_b, edges = load_planck_binned(
        PLANCK_FILE, args.ell_min, args.ell_max
    )
    nbins = ell_b.size
    print(f"[planck] {nbins} bins in {args.ell_min} <= ell <= {args.ell_max}")

    # LCDM baseline (lensed power law)
    p = camb_params(args.lmax)
    p.InitPower.set_params(As=args.As, ns=args.ns)
    dl_lcdm = get_lensed_dl(p, args.lmax)
    dl_lcdm_binned = bin_theory(dl_lcdm, edges)
    chi2_lcdm = chi2(dl_b, err_b, dl_lcdm_binned)
    print(f"[lcdm] chi2 = {chi2_lcdm:.1f} / {nbins} bins "
          f"(chi2/bin = {chi2_lcdm / nbins:.3f})")
    np.savetxt(outdir / "dl_lcdm_lensed.csv",
               np.column_stack([np.arange(dl_lcdm.size), dl_lcdm]),
               delimiter=",", header="ell,DlTT_muK2", comments="")

    sp = SpectrumParams(kappa=args.kappa, As=args.As, ns=args.ns,
                        k0=args.k0, mu_ref="median")

    rows = []
    cases = [(f"a2_{a2:g}", fn, a2) for a2, fn in sorted(A2_FILES.items())]
    cases.append((BEST_HI_RES[0], BEST_HI_RES[1], 8000))
    for label, fname, a2 in cases:
        path = DATA_DIR / fname
        if not path.exists():
            print(f"[skip] {fname} missing")
            continue
        t0 = time.time()
        lambdas, mus = load_genft_input(path)
        k, pz, meta = build_primordial_spectrum(lambdas, mus, sp)
        pparams = camb_params(args.lmax)
        from camb.initialpower import SplinedInitialPower
        pparams.InitPower = SplinedInitialPower()
        pparams.InitPower.effective_ns_for_nonlinear = args.ns
        pparams.InitPower.set_scalar_table(k, pz)
        dl_g = get_lensed_dl(pparams, args.lmax)
        dl_g_binned = bin_theory(dl_g, edges)
        c2 = chi2(dl_b, err_b, dl_g_binned)

        # fractional residual vs LCDM theory (continuity metric)
        lo, hi = int(args.ell_min), int(args.ell_max)
        frac = dl_g[lo:hi + 1] / dl_lcdm[lo:hi + 1] - 1.0
        rms_vs_lcdm = float(np.sqrt(np.mean(frac ** 2)))

        dev = pk_deviation_stats(lambdas, mus, sp)
        rows.append(dict(label=label, a2=a2, chi2_planck=c2,
                         chi2_per_bin=c2 / nbins,
                         dchi2_vs_lcdm=c2 - chi2_lcdm,
                         rms_frac_vs_lcdm_theory=rms_vs_lcdm,
                         pk_dev_rms=dev["dev_rms"],
                         pk_dev_detrended_rms=dev["dev_detrended_rms"],
                         n_modes=int(lambdas.size),
                         runtime_s=round(time.time() - t0, 1)))
        print(f"[{label}] chi2={c2:.1f} dchi2={c2 - chi2_lcdm:+.1f} "
              f"rms_vs_lcdm={rms_vs_lcdm:.3e} "
              f"pk_dev_rms={dev['dev_rms']:.3e} ({rows[-1]['runtime_s']}s)")
        np.savetxt(outdir / f"dl_genft_{label}.csv",
                   np.column_stack([np.arange(dl_g.size), dl_g]),
                   delimiter=",", header="ell,DlTT_muK2", comments="")

    import csv
    keys = list(rows[0].keys())
    with open(outdir / "planck_scan_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    summary = dict(
        planck_file=str(PLANCK_FILE.name),
        nbins=nbins, ell_min=args.ell_min, ell_max=args.ell_max,
        kappa=args.kappa, As=args.As, ns=args.ns, k0=args.k0,
        lmax=args.lmax, mu_ref="median",
        cosmology=dict(H0=67.4, ombh2=0.0224, omch2=0.12, tau=0.054),
        lensed=True, chi2_lcdm=chi2_lcdm,
        camb_version=camb.__version__,
        note="diagonal chi2, symmetrized errors, no covariance/foregrounds",
    )
    with open(outdir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    np.savetxt(outdir / "planck_binned_used.csv",
               np.column_stack([ell_b, dl_b, err_b]),
               delimiter=",", header="ell_eff,Dl_muK2,err_muK2", comments="")
    print(f"[done] results in {outdir}")


if __name__ == "__main__":
    main()
