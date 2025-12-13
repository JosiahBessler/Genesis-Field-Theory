# 📘 Genesis Field Theory (GenesisFT) – CAMB Numerical Validation (Exploratory)
--------------------------------------------------------------------------------------------------------------------------------------------

⚠️ **This folder contains exploratory, non-peer-reviewed numerical validation work.**  
Nothing here should be interpreted as established cosmology or accepted physical theory.

This directory documents **preliminary numerical compatibility testing** between **Genesis Field Theory–derived primordial spectra** and the standard cosmological Boltzmann code **CAMB**.

The goal of this work is **not** to claim physical correctness, but to test whether GenesisFT-derived structures can be *numerically propagated, compared, and stress-tested* within an established cosmological pipeline.

---

## ⚠️ Scientific Disclaimer (Please Read Carefully)

This material represents **early-stage numerical validation only** and was developed through a combination of:

- **Human conceptual development**  
  - All physical motivation, theoretical framing, and interpretation originate from the author.
- **AI-assisted formalization and workflow construction**  
  - Documentation, numerical testing structure, and analysis scripts were AI-assisted.
- **Theory-derived primordial inputs**  
  - GenesisFT produces a primordial scalar power spectrum through a stability-based mode construction.
  - No observational CMB data were used to *construct* the GenesisFT spectrum.
- **Independent numerical propagation using CAMB**
  - CAMB is treated strictly as a numerical engine.
  - No modifications were made to CAMB’s internal physics solvers.
- **Comparative numerical diagnostics**
  - Residuals, convergence, and χ² statistics are used *only* as numerical consistency measures.

This work is **not peer-reviewed**, **not likelihood-fitted**, and **not a cosmological parameter inference**.

---

### Because of this exploratory methodology:

🚫 This is **not** a claim that GenesisFT replaces inflation or ΛCDM  
🚫 This is **not** a best-fit comparison to Planck or other CMB datasets  
🚫 χ² values here are **diagnostic**, not inferential  
🚫 No cosmological parameters are claimed to be measured  
🚫 Numerical agreement does **not** imply physical correctness  

All results should be considered **experimental and provisional**.

---

## 🧭 Purpose of This Folder

This folder exists to answer a **single, narrow question**:

> *Can a GenesisFT-derived primordial spectrum be numerically evolved through CAMB in a stable, convergent, and internally consistent manner?*

It does **not** attempt to answer whether GenesisFT is *true*.

---

## 📂 Contents

### **GenesisFT_CMB_Numerical_Validation_Paper_v11.docx**

A detailed technical report documenting:

- Construction of GenesisFT primordial spectra
- Injection into CAMB via `SplinedInitialPower`
- Multi-resolution convergence tests (L1 / L2 / L3)
- Residual diagnostics (TT / TE / EE)
- Cosmic-variance–weighted χ² comparisons
- Baseline ΛCDM parameter sensitivity checks
- Numerical stability and reproducibility analysis

This document is intended as a **working technical record**, not a finalized journal submission.

---

## 🧪 What Was Tested (High-Level)

1. **GenesisFT primordial spectra construction**
   - Stability-derived mode structure
   - Fixed shape; no observational tuning

2. **CAMB numerical propagation**
   - Standard CAMB release
   - No code modifications
   - Identical cosmological parameters for baseline and GenesisFT runs

3. **Resolution convergence**
   - Increasing spectral resolution (L1 → L2 → L3)
   - Tracking residual decay and χ²/dof convergence

4. **Baseline robustness**
   - Variation of \( n_s \) and \( A_s \)
   - Comparison against fixed GenesisFT reference spectrum

5. **Diagnostics**
   - Fractional residuals
   - χ²/dof (cosmic-variance weighted)
   - Pull statistics

---

## ⚠️ What This Does *Not* Show

❌ That GenesisFT is physically correct  
❌ That GenesisFT explains inflation  
❌ That GenesisFT fits real CMB data  
❌ That CAMB “confirms” GenesisFT  
❌ That ΛCDM is invalid  

This is a **numerical compatibility test only**.

---

## 🧠 Relationship to the UCP Work

- **UCP (galactic scale):**  
  Tests a theory-derived dark-matter halo profile against galaxy rotation curves.

- **CAMB validation (cosmological scale):**  
  Tests whether a GenesisFT-derived primordial spectrum behaves consistently when evolved through a Boltzmann solver.

These are **independent tests** of the same underlying stability-based framework, operating at *very different physical scales*.

Neither validates the other.

---

## 🔍 Why This Is Public (But Not a Preprint)

This material is released publicly in order to:

### **1. Preserve provenance and timestamp**
Establish a clear public record of the numerical work performed.

### **2. Enable scrutiny**
Allow others to examine assumptions, scripts, and diagnostics.

### **3. Allow independent reproduction**
If the results cannot be reproduced, the methodology is flawed.

### **4. Allow time for consolidation**
The author intends to fully understand all terminology, diagnostics, and implications before any formal submission.

---

## 🧠 Author Note

This CAMB validation was performed rapidly following the release of the UCP preprint.  
The author is intentionally **pausing formal publication** to:

- consolidate understanding,
- verify terminology and interpretation,
- and ensure all claims are properly scoped.

This folder should be viewed as a **technical notebook**, not a finished theory paper.

---

## 📌 Citation & Use

If you reference or discuss this CAMB validation work, please cite it informally as:

> **Bessler, J. (2025).**  
> *Genesis Field Theory: CAMB Numerical Validation (Exploratory).*  
> GitHub repository (working document).

Formal citation guidance will be provided **only after** peer-reviewed submission.

---

## 🧭 Final Summary

This folder provides:

✔ transparent numerical testing  
✔ reproducible CAMB workflows  
✔ convergence and robustness diagnostics  
✔ conservative scientific framing  

It does **not** claim correctness — only **numerical viability**.

