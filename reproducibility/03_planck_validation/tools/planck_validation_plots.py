"""Figures for the Planck TT validation scan (run tools/planck_validation.py first).

Usage: python tools/planck_validation_plots.py --outdir outputs/planck_validation_20260704
"""
from __future__ import annotations

import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cscan-dir", default=None,
                    help="optional planck_validation_cscan output dir for fig5")
    args = ap.parse_args()
    out = pathlib.Path(args.outdir)
    plots = out / "plots"
    plots.mkdir(exist_ok=True)

    res = pd.read_csv(out / "planck_scan_results.csv")
    res = res[~res["label"].str.startswith("best_")].sort_values("a2")
    binned = np.loadtxt(out / "planck_binned_used.csv", delimiter=",", skiprows=1)
    ell_b, dl_b, err_b = binned[:, 0], binned[:, 1], binned[:, 2]
    lcdm = np.loadtxt(out / "dl_lcdm_lensed.csv", delimiter=",", skiprows=1)
    ell_t, dl_lcdm = lcdm[:, 0], lcdm[:, 1]
    nbins = ell_b.size
    chi2_lcdm = res["chi2_planck"].iloc[-1] - res["dchi2_vs_lcdm"].iloc[-1]

    # --- Fig 1: chi2 vs a2 ---
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(res["a2"], res["chi2_planck"], "o-", label="GenFT Model A (TT vs Planck 2018)")
    ax.axhline(chi2_lcdm, color="k", ls="--",
               label=f"$\\Lambda$CDM baseline ($\\chi^2$={chi2_lcdm:.1f}/{nbins} bins)")
    ax.axhline(chi2_lcdm + 25, color="r", ls=":", label="$\\Lambda$CDM + 25 ($\\sim$5$\\sigma$)")
    ax.set_xlabel("$a_2$ (kernel-vacuum stabilization)")
    ax.set_ylabel("$\\chi^2$ (66 Planck TT bins, $30\\leq\\ell\\leq2000$)")
    ax.set_title("GenFT Model A vs Planck 2018 binned TT")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(plots / "fig1_chi2_vs_a2.png", dpi=150)
    plt.close(fig)

    # --- Fig 2: spectra for selected a2 with Planck data ---
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(ell_b, dl_b, yerr=err_b, fmt="k.", ms=4, capsize=2,
                label="Planck 2018 binned TT", zorder=5)
    ax.plot(ell_t, dl_lcdm, "k-", lw=1, label="$\\Lambda$CDM (lensed)")
    for a2, color in [(20, "tab:red"), (200, "tab:orange"), (2000, "tab:green"), (8000, "tab:blue")]:
        f = out / f"dl_genft_a2_{a2:g}.csv"
        if f.exists():
            d = np.loadtxt(f, delimiter=",", skiprows=1)
            ax.plot(d[:, 0], d[:, 1], lw=1, color=color, label=f"GenFT $a_2$={a2:g}")
    ax.set_xlim(30, 2000)
    ax.set_ylim(0, None)
    ax.set_xlabel("$\\ell$")
    ax.set_ylabel("$D_\\ell^{TT}$ [$\\mu K^2$]")
    ax.set_title("TT power spectra vs Planck 2018")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(plots / "fig2_tt_spectra.png", dpi=150)
    plt.close(fig)

    # --- Fig 3: fractional residual vs LCDM theory, selected a2 ---
    fig, ax = plt.subplots(figsize=(8, 5))
    for a2, color in [(20, "tab:red"), (200, "tab:orange"), (2000, "tab:green"), (8000, "tab:blue")]:
        f = out / f"dl_genft_a2_{a2:g}.csv"
        if f.exists():
            d = np.loadtxt(f, delimiter=",", skiprows=1)
            ell = d[:, 0].astype(int)
            m = (ell >= 30) & (ell <= 2000)
            frac = d[m, 1] / dl_lcdm[ell[m]] - 1.0
            ax.plot(ell[m], frac, lw=0.8, color=color, label=f"$a_2$={a2:g}")
    # Planck data relative to LCDM, for scale
    interp_lcdm = np.interp(ell_b, ell_t, dl_lcdm)
    ax.errorbar(ell_b, dl_b / interp_lcdm - 1, yerr=err_b / interp_lcdm,
                fmt="k.", ms=4, capsize=2, label="Planck data", zorder=5)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("$\\ell$")
    ax.set_ylabel("$D_\\ell / D_\\ell^{\\Lambda CDM} - 1$")
    ax.set_title("Fractional TT residual vs $\\Lambda$CDM (high-$\\ell$ deficit structure)")
    ax.set_ylim(-0.5, 0.5)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(plots / "fig3_residuals.png", dpi=150)
    plt.close(fig)

    # --- Fig 4: rms residual and P(k) deviation vs a2 ---
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(res["a2"], res["rms_frac_vs_lcdm_theory"], "o-",
              label="rms fractional TT residual vs $\\Lambda$CDM")
    ax.loglog(res["a2"], res["pk_dev_rms"], "s-",
              label="$P_\\zeta(k)$ deviation rms (vs smooth power law)")
    ax.loglog(res["a2"], res["pk_dev_detrended_rms"], "^-",
              label="$P_\\zeta(k)$ detrended (oscillatory) rms")
    ax.set_xlabel("$a_2$")
    ax.set_ylabel("deviation amplitude")
    ax.set_title("Deviation amplitudes vs stabilization $a_2$")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(plots / "fig4_deviation_amplitudes.png", dpi=150)
    plt.close(fig)

    # --- Fig 5: c-scan exclusion (optional) ---
    if args.cscan_dir:
        cs = pd.read_csv(pathlib.Path(args.cscan_dir) / "planck_cscan_results.csv")
        m = cs["c"] > 0
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.loglog(cs.loc[m, "c"], cs.loc[m, "dchi2_vs_lcdm"], "o-",
                  color="tab:blue", label="GenFT Model A ($a_2$=8000)")
        ax.axhline(25, color="r", ls=":", label="$\\Delta\\chi^2$ = 25 exclusion threshold")
        cc = np.geomspace(0.01, 3, 100)
        ax.loglog(cc, 28.0 * (cc / 0.05) ** 2, "k--", lw=0.8,
                  label="$\\Delta\\chi^2 \\propto c^2$ scaling")
        ax.axvline(0.047, color="gray", ls="-.", lw=0.8, label="$c \\approx 0.047$ bound")
        ax.set_xlabel("$c$ (eigenvalue coupling in kernel equation)")
        ax.set_ylabel("$\\Delta\\chi^2$ vs $\\Lambda$CDM (66 Planck TT bins)")
        ax.set_title("Planck 2018 TT exclusion of the GenFT eigenvalue coupling")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3, which="both")
        fig.tight_layout()
        fig.savefig(plots / "fig5_cscan_exclusion.png", dpi=150)
        plt.close(fig)

    print(f"figures written to {plots}")


if __name__ == "__main__":
    main()
