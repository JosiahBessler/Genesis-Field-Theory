# Planck 2018 TT Validation of GenFT Model A

This folder contains everything needed to reproduce the Planck validation scan
(see `../docs/GenesisFT_ModelA_Planck2018_TT_Validation_Report.docx` for the full
report, and `reference_results/` for the expected outputs and figures).

**Result being reproduced:** every Model A configuration in which the operator
eigenvalues influence the primordial spectrum (coupling c = 1, stabilization
a₂ = 0.01 → 12000) is excluded by Planck 2018 binned TT (minimum Δχ² = +15,189
vs ΛCDM over 66 bins; ΛCDM itself scores χ² = 64.0). The c = 0 configuration
matches Planck exactly because its mode weights are constant, i.e., it *is* the
ΛCDM power law. The c-scan bounds the coupling at c ≲ 0.047 (a₂ = 8000).

## Steps

1. Install the pipeline package (from this folder):

   ```
   pip install -e .
   ```

2. Regenerate the model input spectra (not shipped — each file is ~25 MB).
   The generator is deterministic; ~1.58M modes per file.

   a₂ scan (c = 1), one command per value of A2 in
   {0.01, 0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 20, 30, 50, 80, 120, 200, 400,
   600, 800, 1000, 1500, 2000, 3000, 5000, 8000, 12000}:

   ```
   python make_model_a_input.py --ell-max 100 --n-max 12 --a2 <A2> --c 1.0 --out data/model_a_kvac_dense_ell100_n12_a2_<A2>.npz
   ```

   (File-name convention: 0.01 → `0p01`, 0.1 → `0p1`; see the `A2_FILES` map at the
   top of `tools/planck_validation.py` for the exact expected names.)

   c scan at a₂ = 8000, for C in {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
   0.5, 0.7, 1.0, 1.5, 2.0, 3.0}:

   ```
   python make_model_a_input.py --ell-max 100 --n-max 12 --a2 8000 --c <C> --out data/model_a_kvac_dense_ell100_n12_a2_8000_c_<C>.npz
   ```

   (Naming: 0.05 → `0p05`, …, 1.0 → `1`; see `C_FILES` in
   `tools/planck_validation_cscan.py`.)

   High-resolution reference input (used for the December 2025 "best run"
   replication):

   ```
   python make_model_a_input.py --ell-max 120 --n-max 14 --a2 8000 --c 0 --out data/BEST_dense_ell120_n14_a2_8000_c_0.npz
   ```

3. Run the scans (each CAMB case ~3 s; full a₂ scan ~2 min):

   ```
   python tools/planck_validation.py --kappa 0.00118 --lmax 2600
   python tools/planck_validation_cscan.py
   python tools/planck_validation_plots.py --outdir outputs/<scan_dir> --cscan-dir outputs/<cscan_dir>
   ```

4. Compare your `planck_scan_results.csv` / `planck_cscan_results.csv` against
   `reference_results/`. Key check values:

   | Case | χ² (66 bins) | Δχ² vs ΛCDM |
   |---|---|---|
   | ΛCDM baseline | 64.0 | — |
   | a₂ = 12000, c = 1 | 15,252.9 | +15,188.9 |
   | a₂ = 8000, c = 1 | 16,902.1 | +16,838.1 |
   | a₂ = 8000, c = 0.05 | 92.0 | +28.0 |
   | a₂ = 8000, c = 0 | 64.0 | −0.0 |

5. Optional: replicate the original (theory-vs-theory, unlensed) December 2025
   reference run and check rms fractional residual = 1.5405718×10⁻⁵:

   ```
   python -m genftcamb.pipeline --input data/BEST_dense_ell120_n14_a2_8000_c_0.npz --kappa 0.00118 --lmax 3000 --ell-min 30 --ell-max 2000
   ```

## Method caveats (stated in full in the report)

Diagonal χ² with symmetrized errors, no bin-bin covariance or foreground/nuisance
parameters, TT only, cosmology fixed at H0 = 67.4, ω_b = 0.0224, ω_c = 0.12,
τ = 0.054. The ΛCDM calibration (χ²/bin = 0.970) shows this approximation is
adequate at the effect sizes involved (excluded deviations are 10–300% against
~1–3% data errors).
