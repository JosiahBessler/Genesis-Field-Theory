# 📘 Genesis Field Theory (GenesisFT) – Preliminary Research Notes
--------------------------------------------------------------------------------------------------------------------------------------------

⚠️ **This repository contains exploratory, unverified, non-peer-reviewed research.**  
Nothing here should be interpreted as established scientific fact.

This project presents a theory-derived dark-matter halo shape — the **Universal Cored Profile (UCP)** — constructed from the stability equation of a proposed K-operator in the **Genesis Field Theory (GenesisFT)** framework.

The purpose of releasing this material is to invite criticism, correction, and independent verification, *not* to assert correctness.

---

## ⚠️ Scientific Disclaimer (Please Read)

This repository contains research developed through a combination of:

- **Human conceptual development** (all core ideas, motivation, and direction)  
- **AI-assisted formalization** (mathematical exposition, derivations, documentation)  
- **Theory-based derivations**  
  - To the best of my knowledge, the UCP profile used here was **not tuned** using SPARC rotation-curve data.  
  - During development I explored several **lensing-motivated band-mix simulations** to understand what spectral-band combinations could plausibly produce realistic dark-field lensing behavior.  
  - These exploratory checks were **purely theoretical** and did not involve any SPARC galaxy data.  
  - Once the CSV profile was fixed, it was applied uniformly across all galaxies, with only **amplitude (A)** and **scale radius (s)** fitted per galaxy.  
- **Independent numerical testing** (rotation-curve fitting and model comparison)

This work is **preliminary, unverified, and not peer-reviewed**.  
All results should be considered experimental until independently reproduced and validated.

---

### Because of this hybrid methodology:

🚫 The work has **not** been validated by the cosmology or astrophysics community.  
🚫 The mathematics may contain **mistakes**.  
🚫 The code may contain **implementation errors**.  
🚫 Results may reflect misunderstandings by a **non-expert user**.  
🚫 Rotation-curve fits are **not** evidence of correctness.

No scientific claims are made.  
No physical conclusions should be drawn.  
Nothing in this repository should be used as a substitute for established cosmological models.

This is released for transparency and open scientific discussion only.

---

## 🧭 Why This Repository Is Public

This project is released publicly for four reasons:

### **1. To allow experts to find mistakes.**

If errors exist in:

- the K-operator formulation  
- the stability equation  
- the numerical implementation  
- the rotation-curve pipeline  

…releasing the full code allows others to discover and correct them.

---

### **2. To allow independent reproduction.**

If others cannot reproduce the results, something in the methodology is incorrect.

---

### **3. To maintain research integrity.**

Science progresses only when methods and results are openly accessible.

---

### **4. To evaluate whether a stability-derived DM profile can compete with empirical profiles.**

This is a **testable question**, not a claim.

---

## 🧩 How the UCP Was Constructed

The UCP is **not empirically tuned**.  
It is derived from theory as follows:

1. **Propose a K-operator**  
   Based on GenesisFT stability principles, representing how modes evolve.

2. **Solve the stability equation**

       K φ = λ φ

   where:  
   - **positive λ → unstable → decay**  
   - **λ = 0 → stable → surviving mode**

3. **Extract the radial equilibrium solution**  
   The dark field’s lowest-energy stable configuration.

4. **Generate a halo density profile numerically**  
   Producing the universal, fixed-shape “band-derived” CSV halo.

5. **Use the CSV template directly**  
   No galaxy data is used in constructing the UCP profile.

6. **Fit only two parameters per galaxy:**

   - **Amplitude (A)**  
   - **Scale radius (s)**  

Thus the GenesisFT UCP is **purely theory-derived** and globally fixed.

---

## 📊 How Rotation-Curve Tests Work

For each SPARC galaxy:

1. Observed rotation curve is loaded.  
2. Baryonic contribution (disk + bulge) is computed.  
3. Each dark-matter profile is fitted:

   - **GenesisFT UCP** (theory-derived)  
   - Burkert  
   - NFW  
   - Einasto  
   - DC14  
   - Pseudo-Isothermal  
   - Stable-Attractor (SIDM-inspired)

4. For each model we compute:

   - Reduced χ²  
   - AIC / BIC  
   - RAR residuals  
   - Slope and curvature tests  
   - ΔV scatter  
   - *(Optional)* lensing α(R)

---

## 🧪 Replication Instructions

**Install Python + dependencies:**

    pip install numpy pandas scipy lmfit matplotlib astropy tqdm

**Run the full pipeline:**

    python fit_all_multiDM.py

**Results output to:**

- `Output/Tables/MultiDM_all.csv`  
- `Output/Plots/`

All files required for replication are included.

---

## ⚠️ Key Transparency Points

✔ The UCP used here is the **original theory-derived CSV**, not an analytic approximation.  
✔ Re-running the script produces **identical results** — deterministic fit.  
✔ No hidden parameters or galaxy-specific tuning exist.  
✔ UCP is **not designed to outperform** empirical halo models.  
✔ This is a **test**, not a claim.

---

## 🧠 Final Summary

This project provides:

✔ a purely theory-derived universal halo profile  
✔ a complete open-source rotation-curve testing pipeline  
✔ transparent disclaimers and documentation  
✔ an invitation for expert critique and verification

The UCP is:

- derivable  
- reproducible  
- testable  
- **not yet validated**

This project exists **for open scientific analysis**, not advocacy of correctness.

---

## 📖 How to Cite

If you use the Universal Cored Profile (UCP), the accompanying theory paper, or results derived from this repository, please cite the Zenodo preprint:

**Bessler, J. (2025).**  
*The Universal Cored Profile (UCP): Theory, Derivation, and Validation.*  
Zenodo.  
https://doi.org/10.5281/zenodo.17904509

If you also make use of the code or pipeline in this repository, please additionally cite the repository itself.

---

## 📚 SPARC Dataset & Pipeline Acknowledgment

This project uses the publicly available **SPARC rotation-curve dataset**:

**Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016)**  
*AJ, 152, 157 — doi:10.3847/0004-6256/152/6/157*

Elements of the rotation-curve pipeline were adapted from scripts provided on the SPARC website.  
All modifications, extensions, and interpretations are my own.

SPARC website:  
https://astroweb.cwru.edu/SPARC/

