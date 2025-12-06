📘 Genesis Field Theory (GenesisFT) – Preliminary Research Notes
--------------------------------------------------------------------------------------------------------------------------------------------
⚠️ This repository contains exploratory, unverified, non-peer-reviewed research.
Nothing here should be interpreted as established scientific fact.

This project presents a theory-derived dark-matter halo shape — the Universal Cored Profile (UCP) — constructed from the stability equation of a proposed K-operator in the Genesis Field Theory framework.

The purpose of releasing this material is to invite criticism, correction, and independent verification, not to assert correctness.
________________________________________
⚠️ Scientific Disclaimer (Please Read)
--------------------------------------------------------------------------------------------------------------------------------------------
This repository contains research developed through a combination of:

Human conceptual development (all core ideas, motivation, and direction)

AI-assisted formalization (mathematical exposition, derivations, documentation)

Theory-based derivations (no known empirical tuning — To the best of my knowledge, the UCP profile used here was not tuned using SPARC rotation-curve data. During development, I explored several lensing-motivated band-mix simulations to understand what spectral-band combinations could plausibly produce realistic dark-field lensing behavior. These exploratory checks were theoretical and did not involve any SPARC galaxy data. Once the CSV profile was fixed, it was applied uniformly across all galaxies, with only amplitude and scale radius fitted per galaxy.)

Independent numerical testing (rotation-curve fitting and model comparison)

This work is preliminary, unverified, and not peer-reviewed.
All results should be considered experimental until independently reproduced and validated by experts.
________________________________________
Because of this hybrid methodology:

🚫 The work has not been validated by the cosmology or astrophysics community.

🚫 The mathematics may contain mistakes.

🚫 The code may contain implementation errors.

🚫 Results may reflect misunderstandings by a non-expert user.

🚫 Rotation curve fits are not evidence of correctness.

No scientific claims are made.
No physical conclusions should be drawn.
Nothing in this repository should be used as a substitute for established cosmological models.

This is released for transparency and open scientific discussion only.
________________________________________
🧭 Why This Repository Is Public

This project is released publicly for four reasons:

1. To allow experts to find mistakes.

If errors exist in:

the K-operator formulation  
the stability equation  
the numerical implementation  
the rotation-curve pipeline  

…then releasing the full code enables others to identify and correct them.

2. To allow independent reproduction.

If others cannot reproduce the results, that signals an error in the methodology.

3. To maintain research integrity.

Science progresses only when methods and results are openly accessible.

4. To evaluate whether a stability-derived DM profile can compete with empirical profiles.

This is a testable question, not a claim.

________________________________________

🧩 How the UCP Was Constructed

The UCP in this repository is not empirically tuned.
It is derived from theory as follows:

Propose a K-operator  
Based on GenesisFT stability principles, representing how modes evolve.

Solve the stability equation  

K φ = λ φ

with  

positive λ = unstable → decay  
λ = 0 = stable (“surviving modes”)

Extract the radial equilibrium solution  
The dark field’s lowest-energy stable configuration.

Generate a halo density profile numerically  
This produces the “band-derived” CSV halo shape included here.

Use that shape directly  

No galaxy data is used in the K-operator solution.  
No parameters are tuned per galaxy.  

Only two free parameters vary across galaxies:

amplitude A  
scale radius s  

Thus the GenesisFT UCP is purely theory-derived and globally fixed.
________________________________________

📊 How Rotation-Curve Tests Work

For each SPARC galaxy:

Observed rotation curve is loaded.

Baryonic contribution (disk + bulge) is computed.

Each dark-matter profile is fitted:

GenesisFT UCP (theory-derived)  
Burkert  
NFW  
Einasto  
DC14  
Pseudo-Isothermal  
Stable-Attractor (SIDM)

For each model we compute:

reduced χ²  
AIC / BIC  
RAR residuals  
slope and curvature tests  
ΔV scatter  
(optional) lensing α(R) — experimental
________________________________________

🧪 Replication Instructions

Install Python + dependencies:  
pip install numpy pandas scipy lmfit matplotlib astropy tqdm

Run:  
python fit_all_multiDM.py

Inspect:  
Output/Tables/MultiDM_all.csv

All files required for replication are provided.
________________________________________

The main test script is:  
python fit_all_multiDM.py

Results appear in:  
Output/Tables/MultiDM_all.csv

Plots appear in:  
Output/Plots/
________________________________________

⚠️ Key Transparency Points:

✔ The UCP used here is the original theory-derived CSV, not an analytic approximation  
(The analytic version is still being developed and must perfectly match the CSV to be used.)

✔ Re-running the script produces identical results  
No hidden parameters or galaxy-specific tuning exist.

✔ The GenesisFT UCP is not designed to win  
If it performs poorly compared to NFW, Burkert, or Einasto, that is acceptable and expected.

This is a test, not a claim.
________________________________________

🧠 Final Summary

This project presents:

✔ a purely theory-derived universal halo profile  
✔ a complete open-source testing pipeline  
✔ transparent disclaimers and documentation  
✔ an invitation for expert review, critique, and correction

The UCP is:

derivable  
reproducible  
testable  
not yet validated

The goal is not to assert correctness, but to enable open scientific analysis.

📚 **SPARC Dataset & Pipeline Acknowledgment**

This project uses the publicly available **SPARC rotation-curve dataset** (Lelli, McGaugh & Schombert 2016).  
Elements of the rotation-curve fitting pipeline were adapted from the analysis scripts provided on the SPARC website.  
All modifications, extensions, and interpretations are my own.

Please cite SPARC if using these data in scientific work. Website: https://astroweb.cwru.edu/SPARC/
