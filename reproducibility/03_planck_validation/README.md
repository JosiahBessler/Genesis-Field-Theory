# GenesisFT CAMB Pipeline

An end-to-end Python pipeline to ingest GenesisFT spectral outputs (eigenvalues `lambda` and weights `mu`), build a primordial curvature power spectrum, run CAMB, and compare against a baseline ΛCDM spectrum.

## Installation

```bash
pip install -e .
```

Dependencies: `numpy`, `pandas`, `matplotlib`, and `camb`.

The default run uses a linear, unlensed CAMB configuration with `lmax=1200` for stability. Once you verify the pipeline, you can increase `lmax` or enable lensing as needed.

## Generate dummy inputs

Create plausible test inputs in `data/`:

```bash
python make_dummy_input.py
```

This writes `data/dummy_genft.npz` and `data/dummy_genft.csv` with sorted positive eigenvalues and weights.

## Run the pipeline

Example using the generated NPZ input:

```bash
python -m genftcamb.pipeline --input data/dummy_genft.npz
```

You can also provide the CSV version:

```bash
python -m genftcamb.pipeline --input data/dummy_genft.csv
```

Outputs are written to `outputs/run_YYYYMMDD_HHMMSS/` by default and include:

- `pk.csv` and `pk.npz`: primordial spectrum table (k, Pzeta)
- `cls_genft.*`: CAMB CMB power spectra using the GenesisFT-derived table
- `cls_lcdm.*`: baseline ΛCDM spectra
- `report.md`: summary of parameters and comparison metrics
- `plots/`: PNG diagnostics for the primordial spectrum and TT spectra/residuals

## CLI options

Run `python -m genftcamb.pipeline --help` to see all options. Key overrides include `--kappa`, `--As`, `--ns`, `--k0`, `--mu-ref`, `--lmax`, `--ell-min`, `--ell-max`, and `--frac-error`. You can skip either CAMB run with `--baseline-only` or `--genft-only`.

## Notes

- Nonpositive `lambda` or `mu` values are removed with warnings.
- The pipeline ensures the `k` grid is strictly increasing for CAMB and clips extreme values to avoid numerical issues.
- Cosmology defaults: H0=67.4, ombh2=0.0224, omch2=0.12, tau=0.054, lmax=1200 (linear, unlensed default run).
