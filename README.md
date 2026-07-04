# Genesis Field Theory (GenesisFT)
### A Human–AI Assisted Exploration of Emergent Cosmology

This repository presents **Genesis Field Theory (GenesisFT)** — a speculative emergent cosmology framework developed through a **human-guided, AI-assisted research process**.  
The core *ideas*, *conceptual architecture*, *research direction*, and *experimental goals* were created by me.  
The *mathematical exposition*, *formal derivations*, and *scientific writing* were generated with the assistance of large-language-model AI tools under my direction and curation.

---

## 📄 UCP Theory Paper (Preprint)

The complete theoretical derivation and validation of the **Universal Cored Profile (UCP)** — including the stability equation, construction of the UCP template, rotation-curve methodology, and SPARC multi-model comparison — is available here:

👉 **[UCP Theory Paper (PDF)](UCP/UCP_Theory_Paper.pdf)** · **[Revised v2 (July 2026)](UCP/UCP_Theory_Paper_v2.docx)**

[![Paper DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17904509.svg)](https://doi.org/10.5281/zenodo.17904509)  
[![Software DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17914774.svg)](https://doi.org/10.5281/zenodo.17914774)

This document is the **primary scientific reference** for the UCP portion of this repository and is intended for **transparency, reproducibility, and open review**.

> **v2 revision note (July 2026):** the v2 paper corrects the characterization of UCP relative to the GenesisFT operator eigenmodes (the verified lensing-consistent MID-band eigenmode differs from the UCP template by relative L2 ≈ 0.72), reports the fit statistics precisely (UCP is the most frequent best-fit model on SPARC — 62/175 galaxies — while Burkert attains a slightly better median χ²/dof, 0.472 vs 0.521), and fixes reference details. The v1 PDF is retained for the DOI record.

---

## 🔬 Reproducibility Package & Validation Status (July 2026)

The [`reproducibility/`](reproducibility/) folder contains code, data, reference outputs, and
step-by-step instructions to independently reproduce the program's main computational
results — including a **negative result published deliberately**:

| Test | Result | Reproduction status |
|---|---|---|
| [SPARC halo comparison](reproducibility/01_halo_comparison/) (175 galaxies, 5 models) | UCP most frequent best-fit (62/175); competitive, not universally best | Re-run 2026-07-04: output CSV **byte-identical** to archive (SHA-256 verified) |
| [Eigenmode verification](reproducibility/02_eigenmode_verification/) | Verified MID-band eigenmode differs from UCP (rel. L2 = 0.7197) | Automated check script reproduces all values |
| [Planck 2018 TT validation](reproducibility/03_planck_validation/) | **Falsification:** every configuration where the operator spectrum influences P(k) is excluded (min Δχ² = +15,189 / 66 bins); the Planck-consistent c = 0 limit is exactly ΛCDM by construction; coupling bound c ≲ 0.047 | Deterministic pipeline; inputs regenerate bit-identically |

Full analysis: [Planck validation report](reproducibility/docs/GenesisFT_ModelA_Planck2018_TT_Validation_Report.docx).
See [`reproducibility/README.md`](reproducibility/README.md) for setup and exact commands.

**What this means:** the CMB temperature spectrum now acts as a *constraint* on GenesisFT
rather than a validation of it — any distinctive spectral signature must lie below ≈0.6%
in the primordial spectrum. The framework's empirically competitive results remain at
galaxy scale (rotation curves, halo eigenmode structure). An earlier GenFT-vs-MOND
comparison was retracted on methodological grounds the day it was run and its numbers
should not be cited.

---

## 🔷 What This Project Is

GenesisFT is an exploration of:

- emergent spacetime from spectral geometry  
- a generative operator whose eigenmodes produce visible matter and dark fields  
- an arrow of time emerging from kernel condensation  
- dark-field halos built from uncondensed spectral modes  
- oscillatory features in the matter power spectrum P(k)  
- dark-field time-dilation effects  
- gravitational-wave echo phenomenology  
- rotation-curve predictions across the SPARC dataset  

These items represent **areas of exploration and hypothesis generation**, not established results.

**This is not a validated physical theory**, but a **conceptual experiment in using AI to help formalize a new cosmological framework**.

---

## 🔷 Human Contribution

My primary role in this project has been **system design, numerical execution, documentation, and orchestration of AI-assisted formalization**. I am actively learning the underlying **cosmological and astrophysical formalism**, and this repository reflects an **early-stage integration of theory, numerics, and tooling**.

I contributed:



- all **core ideas**, conceptual structures, and research motivations  
- the **framework** for GenesisFT’s operator, kernel, and sector decomposition  
- the **halo model design**, including cored-template behavior  
- the **rotation-curve test design** (DM-only, DM+baryons, extreme DM, etc.)  
- the **rotation-curve modeling and comparisons** across the SPARC dataset using AI-generated fitting tools  
- the **curation, organization, and integration** of all documents  

I directed the AI step-by-step to build the mathematical formalism and scientific language around these ideas.

---

## 🔷 AI Contribution

Large-language-model AI tools were used to:

- generate formal mathematical exposition  
- write proofs, derivations, and operator definitions  
- translate conceptual ideas into scientific formatting  
- draft sections of text, explanations, and terminology  
- help structure long documents and ensure internal consistency  

I do **not** claim to have independently derived every equation or mathematical detail.  
All content has been reviewed, curated, and shaped by me, but not rigorously verified.

---

## ⚠️ Disclaimer

This project is a **speculative, exploratory exercise** in theoretical cosmology and AI-assisted scientific writing.  
Nothing here should be assumed to be:

- experimentally confirmed  
- mathematically validated  
- scientifically correct  
- a replacement for standard cosmological models  

GenesisFT is presented as an **open-ended conceptual framework** rather than an established scientific theory.  
None of the models, equations, or results in this repository constitute validated physical laws.  
All rotation curve fits, halo models, and numerical results are experimental and should not be interpreted as empirical confirmation.  
This work is presented for transparency, reproducibility, and open scientific discussion.

---

## 🔷 License

This project is licensed under **Creative Commons Attribution 4.0 International (CC BY 4.0)**.  
You are free to share and adapt the material for any purpose, even commercially, provided that proper attribution is given.

See the `LICENSE` file for full details.

---

## 🔷 Contact / Contribution

Feedback, critiques, and extensions are welcome.  
Contact Email: **GenesisFieldTheory@outlook.com**

The intent of publicizing this repository is to allow individuals with expertise to review, test, and expand on the ideas presented here.

If you use this material, please cite the repository and credit this project.

---

## 🔷 Notes

In the process of rebranding due to existing uses of the acronym GenFT. Converting to **GenesisFT**.  
This theoretical framework is still in development.

---

## 🌌 Conceptual Narrative Overview (Non-Technical)

The following section provides an **intuitive, non-technical narrative description** of Genesis Field Theory.  
It is intended to convey conceptual motivation and interpretive structure, **not** to serve as a scientific derivation or empirical claim.

---

### What Is Genesis Field Theory?

Genesis Field Theory (GenesisFT) is built on one central idea:

> “The universe always evolves toward the most stable state available.”

GenesisFT extends this principle beyond ordinary physics by proposing that, prior to the emergence of spacetime, matter, or light, the underlying field existed in a state of perfect symmetry:

- no waves  
- no motion  
- no time  
- absolute stability  

Because the field is assumed to possess infinite structure, even extraordinarily rare fluctuations are permitted.  
Eventually, one such fluctuation occurred — an extremely improbable ripple that broke the symmetry.

This rupture triggered a release of dynamics (qualitatively similar to the Big Bang).  
Once broken, the system could not return to perfect neutrality and instead evolved toward the next most stable configuration available.

Everything that exists today — space, time, dark matter, galaxies — is interpreted as emerging from this relaxation toward stability.

---

### 🌐 Why Stability Matters

Stability-driven behavior appears throughout nature:

- A marble rolls to the bottom of a bowl  
- A stretched rubber band snaps back  
- A bubble forms a sphere  

GenesisFT proposes that the universe behaves analogously, but on a cosmic scale.

The initial fluctuation produced many modes in the field:

- **Unstable modes** decayed  
- **Stable modes** persisted  

The surviving modes are interpreted as giving rise to:

- dark matter  
- large-scale structure  
- spacetime geometry  
- effective physical laws  

---

### 🌩️ The Big Fluctuation → Chaos → Stability

GenesisFT frames cosmic evolution in three conceptual phases:

1. **Perfect Stability** — the symmetric, motionless field  
2. **The Fluctuation** — a rare symmetry-breaking event  
3. **Return Toward Stability** —  
   - unstable modes fade  
   - stable modes persist  
   - structure emerges  

Within this narrative framework, **dark matter is interpreted as one class of surviving stable modes**.
