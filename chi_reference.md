# CHI Reference Frame — ECSoC v2 (Revised)

> **This document describes the revised CHI reference frame introduced in the paper.**  
> Prior ECSoC formulations placed the healthy operating point near CHI ≈ 0.  
> New NSR Healthy Database data (N=54) necessitate a systematic revision.

---

## Formula

```
CHI = 2 × (α₁ − α₂)
```

where α₁ = DFA short-scale exponent (4–16 beats) and α₂ = long-scale exponent (16–64 beats).

---

## Revised Reference Values

| Population | N | Mean CHI ± SD | Median | % Supercritical (CHI > 0) |
|------------|---|--------------|--------|--------------------------|
| **NSR Healthy** | 54 | **+0.440 ± 0.441** | +0.449 | **87%** |
| CAST baseline | 734 | +0.113 ± 0.461 | — | 65.7% |
| CHF2 NYHA 1–3 | 29 | −0.243 ± 0.671 | −0.280 | 27.6% |
| MUSIC SCD | 30 | (subcritical) | — | — |
| SVTDB pre-arrhythmic | 86 | −0.453 ± 0.707 | −0.483 | 27% |

---

## Domain Classification

```
CHI >  0.05  →  Supercritical   (dominant short-scale fluctuations)
|CHI| ≤ 0.05 →  Near-critical   (critical manifold transition zone)
CHI < -0.05  →  Subcritical     (suppressed short-scale dynamics)
```

---

## Key Revision Points

### 1. Healthy baseline is supercritical, not critical
Prior formulations: healthy operating point ≈ CHI = 0  
**Revised**: healthy operating point ≈ CHI = +0.440 (87% supercritical in NSR, N=54)

### 2. Disease direction is subcritical
CHF (NYHA 1–3) shows mean CHI = −0.243 — a shift *toward* subcritical, not supercritical.  
This is consistent with established HRV literature (reduced variability in heart failure).

### 3. Supercritical CHI alone does not indicate risk
The Metastable Supercritical State (CHI > 0, R² high, ΔP(III) ≈ 0) is the **healthy norm**.  
Supercritical CHI requires Layer 2 (trajectory) and Layer 3 (R²) for risk context.

### 4. Phase Vb is pathologically supercritical
Phase Vb (CHI > 0 + R² < 0.93) is distinct from healthy supercritical:
- Healthy: CHI > 0, R² > 0.99 (mean 0.9952), no trajectory activation
- Vb: CHI > 0, R² < 0.93, structural scaling breakdown

---

## CHI Gradient Across Disease Severity

```
NSR Healthy  →  CAST  →  CHF2 NYHA 1-3  →  Advanced disease
  +0.440         +0.113       −0.243            < −0.5
   87% SC         66% SC       28% SC            ~30% SC
```

SC = supercritical. This four-point gradient provides a quantitative framework  
for tracking autonomic disease progression (cross-sectional; exploratory).

---

## Drug Effect (CAST, N=734 paired subjects)

Antiarrhythmic drug administration:
- **ΔCHI = −0.138 ± 0.484** (Wilcoxon p < 0.0001) — significant CHI reduction
- **ΔR² = −0.0004** (p = 0.592, NS) — no significant change in structural integrity

This supports **Layer 1 / Layer 3 independence**: drugs affect dynamical position (CHI)  
without altering scaling manifold integrity (R²). *[Exploratory]*

---

## Coding Reference

```r
# R
chi <- compute_chi(alpha1, alpha2)       # from 02_chi_calculation.R
domain <- classify_chi_domain(chi)
interp <- interpret_chi(chi)             # relative to NSR healthy reference
```

```matlab
% MATLAB
[chi, domain, interpretation] = chi_calculation(alpha1, alpha2);
```

---

*All values exploratory. Prospective validation required.*
