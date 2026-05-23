# ECSoC Framework — Emergent Criticality and Self-Organized Collapse

> **Hypothesis-Generating Framework for Cardiovascular Dynamical Analysis**
>
>　Okabe H. Heterogeneous Dynamical Pathways Preceding Ventricular Arrhythmia: 
A Multi-Cohort Phase-Space Analysis of Cardiovascular Criticality. 2026.　 

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language: R](https://img.shields.io/badge/Language-R%20%3E%3D%204.0-blue)](https://www.r-project.org/)
[![Data: PhysioNet](https://img.shields.io/badge/Data-PhysioNet%20Open%20Access-green)](https://physionet.org/)

---

## ⚠️ Important Disclaimer

All analyses in this repository are **exploratory and hypothesis-generating**.  
No prospective validation has been performed. Results **must not** be used for clinical decision-making.

---

## Overview

ECSoC (Emergent Criticality and Self-Organised Collapse) is a 
multi-domain dynamical framework for characterising physiological 
systems approaching catastrophic failure.

The framework proposes that pre-collapse dynamics share a 
common three-layer mathematical structure across organ systems 
— independent of the specific physiological mechanism involved.

### Published Evidence

| Paper | Domain | Journal | Status |
|-------|--------|---------|--------|
| Paper 1 | Cardiac (9 cohorts, N>1,500) | Frontiers in Network Physiology | Under Review |
| Paper 2 | Theoretical (Langevin/order parameter) | Frontiers in Network Physiology | Accepted 2026 |
| Paper 3 | Neural EEG (5 cohorts) | [Journal] | In Preparation |

Cross-domain replication: the same three-layer dissociation 
(CHI scaffold / β traversal / Phase V structural marker) was 
observed independently in cardiac and neural datasets.
---

## Repository Structure

```
ECSoC-Framework/
├── R/
│   ├── 01_dfa.R                  # DFA implementation (α₁, α₂, R²)
│   ├── 02_chi_calculation.R      # CHI computation and phase classification
│   ├── 03_phase_classification.R # Phase I–V assignment
│   ├── 04_trajectory_analysis.R  # Layer 2: ΔP(III), Stage 1/2 detection
│   ├── 05_cohort_analysis.R      # Multi-cohort summary statistics
│   └── 06_figures.R              # Figure generation (ggplot2)
├── matlab/
│   ├── dfa.m                     # DFA implementation
│   ├── chi_calculation.m         # CHI computation
│   ├── phase_classification.m    # Phase assignment
│   └── run_ecsoc.m               # Main analysis script
├── data/
│   └── sample/
│       ├── sample_rr_normal.csv  # Simulated healthy RR intervals (N=500)
│       ├── sample_rr_chf.csv     # Simulated CHF-like RR intervals (N=500)
│       └── README_data.md        # Data format specification
├── docs/
│   ├── methods.md                # Detailed methods description
│   ├── phase_definitions.md      # Phase I–V classification criteria
│   └── chi_reference.md         # CHI reference frame (revised)
├── figures/                      # Output directory for generated figures
├── tests/
│   └── test_dfa.R                # Unit tests for DFA implementation
└── README.md
```

---

## Quick Start (R)

### Installation

```r
# Required packages
install.packages(c("ggplot2", "dplyr", "tidyr", "purrr", "pracma", "scales"))
```

### Run on sample data

```r
source("R/01_dfa.R")
source("R/02_chi_calculation.R")
source("R/03_phase_classification.R")

# Load sample RR interval data (in milliseconds)
rr <- read.csv("data/sample/sample_rr_normal.csv")$rr_ms

# Compute DFA exponents
result <- compute_dfa(rr, scale_min = 4, scale_max_short = 16, scale_max_long = 64)
cat("α₁ =", result$alpha1, "\n")
cat("α₂ =", result$alpha2, "\n")
cat("R² =", result$r_squared, "\n")

# Compute CHI
chi <- compute_chi(result$alpha1, result$alpha2)
cat("CHI =", chi, "\n")

# Classify phase
phase <- classify_phase(result$alpha1, result$r_squared)
cat("Phase:", phase$label, "\n")
```

---

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Short-scale DFA | 4–16 beats | α₁ computation window |
| Long-scale DFA | 16–64 beats | α₂ computation window |
| Phase V threshold | R² < 0.93 | Structural breakdown criterion |
| CHI supercritical | CHI > 0 | Dominant short-scale fluctuations |
| CHI subcritical | CHI < 0 | Suppressed short-scale dynamics |
| ΔP(III) Late-rise | ≥ +0.15 | Stage 2 trajectory criterion |
| ΔP(III) Suppress | < −0.10 | Stage 1 trajectory criterion |

---

## Datasets Used in the Paper

All datasets are publicly available via [PhysioNet](https://physionet.org/):

| Cohort | N | PhysioNet ID | Role |
|--------|---|--------------|------|
| MUSIC | 100 | — | Primary derivation cohort |
| CUDB | 35 | `cudb` | Acute peri-event |
| MVEDB | 22 | `mvedb` | Long-term ambulatory |
| SVTDB | 76 valid | `svtdb` | Trajectory analysis |
| BIDMC CHF | 15 | `chfdb` | External replication |
| CAST | 1,543 | `cast-rr` | Large-scale validation |
| NSR Healthy | 54 | `nsrdb` | Healthy reference |
| CHF2 | 29 | `chf2db` | Intermediate CHF control |

See `docs/methods.md` for data download and preprocessing instructions.

---

## Citation

**Paper 2 (Accepted):**
Okabe H. Empirically Constrained Order Parameter Dynamics in 
Cardiovascular Criticality: A Synergetic Langevin Framework for 
Arrhythmic Transitions with Cross-Cohort Parameter Estimation 
and Kramers Escape-Time Validation.
*Frontiers in Network Physiology.* Accepted 2026.

**Paper 1 (Under Review):**
Okabe H. Heterogeneous Dynamical Pathways Preceding Ventricular 
Arrhythmia: A Multi-Cohort Phase-Space Analysis of Cardiovascular 
Criticality. *Frontiers in Network Physiology.* Under Review 2026.


```

---

## License

MIT License. See [LICENSE](LICENSE).

Data are from PhysioNet (open access). Please cite PhysioNet appropriately:  
Goldberger AL, et al. PhysioBank, PhysioToolkit, and PhysioNet. *Circulation*. 2000;101(23):e215–e220.
