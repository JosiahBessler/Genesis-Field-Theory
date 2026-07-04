"""CLI entry point for GenesisFT CAMB pipeline."""
from __future__ import annotations

import argparse
import datetime as dt
import logging
import pathlib
from typing import Optional

from .build_pk import SpectrumParams, build_primordial_spectrum
from .io import describe_parameters, load_genft_input, save_cls, save_pk, save_report
from .plots import plot_pzeta, plot_tt, plot_tt_residual
from .run_camb import Cosmology, run_camb_baseline, run_camb_with_table
from .score import ScoreConfig, fractional_residual, score_tt

logger = logging.getLogger(__name__)


def _default_outdir() -> pathlib.Path:
    stamp = dt.datetime.now().strftime("run_%Y%m%d_%H%M%S")
    return pathlib.Path("outputs") / stamp


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="GenesisFT CAMB pipeline")
    parser.add_argument("--input", required=True, help="Path to GenesisFT NPZ or CSV input.")
    parser.add_argument("--outdir", default=None, help="Output directory (default timestamped under outputs/)")
    parser.add_argument("--kappa", type=float, default=1.0, help="Scaling applied to sqrt(lambda) to obtain k.")
    parser.add_argument("--As", type=float, default=2.1e-9, help="Scalar amplitude A_s")
    parser.add_argument("--ns", type=float, default=0.965, help="Scalar spectral index n_s")
    parser.add_argument("--k0", type=float, default=0.05, help="Pivot scale k0 (Mpc^-1)")
    parser.add_argument("--mu-ref", default="median", help="Reference mu value; 'median' or numeric")
    parser.add_argument("--lmax", type=int, default=1200, help="Maximum multipole for CAMB outputs")
    parser.add_argument("--ell-min", type=int, default=30, help="Minimum ell for scoring")
    parser.add_argument("--ell-max", type=int, default=2000, help="Maximum ell for scoring")
    parser.add_argument("--frac-error", type=float, default=0.05, help="Fractional error model for chi2 metric")
    parser.add_argument("--baseline-only", action="store_true", help="Only run baseline ΛCDM")
    parser.add_argument("--genft-only", action="store_true", help="Only run GenesisFT table")
    return parser.parse_args(argv)


def _setup_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")


def _parse_mu_ref(value: str) -> str | float:
    try:
        return float(value)
    except ValueError:
        return value


def main(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)
    _setup_logging()

    outdir = pathlib.Path(args.outdir) if args.outdir else _default_outdir()
    outdir.mkdir(parents=True, exist_ok=True)
    logger.info("Using output directory %s", outdir)

    lambdas, mus = load_genft_input(pathlib.Path(args.input))

    spectrum_params = SpectrumParams(
        k0=args.k0,
        As=args.As,
        ns=args.ns,
        kappa=args.kappa,
        mu_ref=_parse_mu_ref(args.mu_ref),
    )
    ell_max_eff = min(args.ell_max, args.lmax)
    score_config = ScoreConfig(ell_min=args.ell_min, ell_max=ell_max_eff, frac_error=args.frac_error)
    cosmo = Cosmology()

    ell_genft = tt_genft = ee_genft = te_genft = bb_genft = None
    if not args.baseline_only:
        k, pzeta, meta = build_primordial_spectrum(lambdas, mus, spectrum_params)
        save_pk(outdir, k, pzeta)
        plot_pzeta(outdir, k, pzeta)
        ell_genft, tt_genft, ee_genft, te_genft, bb_genft = run_camb_with_table(cosmo, args.lmax, k, pzeta)
        save_cls(outdir, "genft", ell_genft, tt_genft, ee_genft, te_genft, bb_genft)
    else:
        meta = {}

    ell_lcdm = tt_lcdm = ee_lcdm = te_lcdm = bb_lcdm = None
    if not args.genft_only:
        ell_lcdm, tt_lcdm, ee_lcdm, te_lcdm, bb_lcdm = run_camb_baseline(cosmo, args.lmax, args.As, args.ns)
        save_cls(outdir, "lcdm", ell_lcdm, tt_lcdm, ee_lcdm, te_lcdm, bb_lcdm)

    report_lines = ["# GenesisFT CAMB run", ""]
    report_lines.append("## Parameters")
    report_lines.append(
        describe_parameters(
            [
                ("input", args.input),
                ("outdir", str(outdir)),
                ("kappa", args.kappa),
                ("As", args.As),
                ("ns", args.ns),
                ("k0", args.k0),
                ("mu_ref", args.mu_ref),
                ("lmax", args.lmax),
                ("ell_min", args.ell_min),
                ("ell_max", ell_max_eff),
                ("frac_error", args.frac_error),
                ("meta", meta),
            ]
        )
    )

    if ell_genft is not None and ell_lcdm is not None:
        plot_tt(outdir, ell_genft, tt_genft, tt_lcdm)
        residual = fractional_residual(tt_genft, tt_lcdm)
        plot_tt_residual(outdir, ell_genft, residual, args.ell_min, ell_max_eff)
        metrics = score_tt(ell_genft, tt_genft, tt_lcdm, score_config)
        report_lines.append("")
        report_lines.append("## Metrics (TT)")
        report_lines.append(describe_parameters(metrics.items()))
    else:
        report_lines.append("")
        report_lines.append("## Metrics")
        report_lines.append("Not available (only one spectrum computed).")

    save_report(outdir, "\n".join(report_lines))
    logger.info("Run complete. Outputs stored in %s", outdir)


if __name__ == "__main__":
    main()
