# Genesis Field Theory (GenFT) — Reproducibility Package

**Author:** Josiah Bessler
**Package version:** 1.0 (2026-07-04)
**Status:** Independent research. Not peer reviewed. Everything in this repository is
offered for scrutiny, replication, and refutation.

This package contains the code, data, and reference outputs for the three main
computational results of the GenFT research program. Each was independently
re-executed on 2026-07-04 and reproduced exactly (details below and in each
subfolder). Both positive and negative results are included — the Planck CMB
test in `03_planck_validation` is a **falsification** of the tested model branch,
published here deliberately.

## Contents and headline results

| Folder | What it tests | Result |
|---|---|---|
| `01_halo_comparison` | GenFT/UCP halo template vs NFW, Burkert, Pseudo-Isothermal, Einasto on 175 SPARC rotation curves (DM + baryons) | UCP is the most frequent best-fit model (62/175 galaxies, 35.4%); Burkert has a slightly better median χ²/dof (0.472 vs 0.521). Competitive, not universally best. |
| `02_eigenmode_verification` | The optimized MID-band eigenmode vs the UCP template (Dec 2025 analysis) | The verified eigenmode differs from UCP by relative L2 = 0.7197 — UCP is an approximation, not the exact operator eigenmode. |
| `03_planck_validation` | GenFT Model A TT spectra vs the real Planck 2018 binned TT power spectrum | **Falsification of the coupled branch:** all configurations where the operator spectrum influences P(k) are excluded (min Δχ² = +15,189 over 66 bins); the only Planck-consistent configuration (c = 0) is exactly ΛCDM by construction. Coupling bound c ≲ 0.047. |
| `docs` | Full research report for the Planck validation | Methods, figures, exclusion bounds, replication appendix. |

## Reproduction status (2026-07-04)

- **Halo comparison:** re-run from scratch → output CSV is **byte-identical** to the
  archived reference (SHA-256
  `72ECB22339BF0EC2D3A4249BE10D444F7DDC7B802EA1EC595E9ECC2316F058C5`).
  All 875 fits and all 175 per-galaxy winners identical.
- **Eigenmode numbers:** recomputed from the stored mode-shape table —
  digit-for-digit agreement (`verify_eigenmode_numbers.py` automates this).
- **CAMB reference run:** re-executed → rms fractional TT residual
  1.5405718×10⁻⁵, numerically identical to the December 2025 archived run.

Note on scope: reproduction confirms the results are deterministic properties of the
code and data — it does not by itself establish physical correctness of the theory.

## Setup

Python 3.11+ recommended.

```
pip install -r requirements.txt
```

Exact-hash reproduction was verified on Windows 11 / Python 3.11.9 with the pinned
versions. Other platforms should agree to high numerical precision (float
last-digit differences may change file hashes but not conclusions).

## Quick start — reproduce everything

```
# 1. Halo comparison (~ minutes; writes halo_dm_baryons_results.csv + plots)
cd 01_halo_comparison
python fit_all_models_dm_baryons.py
#   compare to reference_results/halo_dm_baryons_results.csv

# 2. Eigenmode verification (instant; checks stored numbers)
cd ../02_eigenmode_verification
python verify_eigenmode_numbers.py

# 3. Planck validation (needs camb; inputs regenerated first — see folder README)
cd ../03_planck_validation
pip install -e .
python make_model_a_input.py --ell-max 100 --n-max 12 --a2 8000 --c 1.0 --out data/model_a_kvac_dense_ell100_n12_a2_8000.npz
#   ... (see 03_planck_validation/README_VALIDATION.md for the full input grid)
python tools/planck_validation.py --kappa 0.00118 --lmax 2600
python tools/planck_validation_cscan.py
python tools/planck_validation_plots.py --outdir outputs/<scan dir> --cscan-dir outputs/<cscan dir>
#   compare to reference_results/planck_scan_results.csv and planck_cscan_results.csv
```

## Data provenance and attribution

- **SPARC rotation curves** (`01_halo_comparison/Rotmod_LTG/`): Lelli, F., McGaugh,
  S. S., & Schombert, J. M. (2016), AJ, 152, 157. Publicly available at
  http://astroweb.cwru.edu/SPARC/ — redistributed here for reproducibility with
  attribution; please cite the SPARC paper in any derived work.
- **Planck 2018 binned TT spectrum**
  (`03_planck_validation/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt`): ESA Planck
  Legacy Archive, mirrored by IRSA
  (https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/).
  © ESA, Planck Collaboration; cite Planck 2018 results (A&A 641, A6, 2020).
- **CAMB**: Lewis, A., Challinor, A., & Lasenby, A. (2000), ApJ, 538, 473;
  https://camb.info.

## Honest-disclosure notes

- The GenFT vs MOND comparison sometimes referenced in early project material was
  retracted the day it was run (mismatched fit conditions) and is **not** included
  here; its numbers should not be cited.
- The December 2025 CAMB "best run" matches ΛCDM because its configuration (c = 0)
  makes the mode weights exactly constant, reducing the primordial spectrum to the
  ΛCDM power law by construction. See the report in `docs/` for the full analysis.
- All statistical caveats (diagonal χ², TT-only, fixed cosmology) are stated in the
  report and in the tool docstrings.

## License

Code: (to be chosen by the author before publication — MIT suggested).
SPARC and Planck data remain under their original terms (see above).
