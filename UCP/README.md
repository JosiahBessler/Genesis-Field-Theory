📘 Genesis Field Theory (GenFT) – Preliminary Research Notes
________________________________________
⚠️ This repository contains exploratory, unverified, non-peer-reviewed research.
Nothing here should be interpreted as established scientific fact.

This project presents a theory-derived dark-matter halo shape — the Universal Cored Profile (UCP) — constructed from the stability equation of a proposed K-operator in the Genesis Field Theory framework.

The purpose of releasing this material is to invite criticism, correction, and independent verification, not to assert correctness.
________________________________________
⚠️ Scientific Disclaimer (Please Read)

This repository includes research created through a combination of:

Human conceptual development (all core ideas and direction)

AI-assisted formalization (mathematical exposition, documentation)

Theory-based derivations (no empirical tuning)

Independent numerical testing (rotation curve fits)
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
Based on GenFT stability principles, representing how modes evolve.

Solve the stability equation

Modes satisfying

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

Thus the GenFT UCP is purely theory-derived and globally fixed.
________________________________________

📊 How Rotation-Curve Tests Work

For each SPARC galaxy:

Observed rotation curve is loaded.

Baryonic contribution (disk + bulge) is computed.

Each dark-matter profile is fitted:

GenFT UCP (theory-derived)

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

✔ The GenFT UCP is not designed to win

If it performs poorly compared to NFW, Burkert, or Einasto, that is acceptable and expected.

This is a test, not a claim.
________________________________________

🧪 Replication Instructions

Install Python + dependencies

Run:

python fit_all_multiDM.py


Inspect:

Output/Tables/MultiDM_all.csv


All files required for replication are provided.

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
