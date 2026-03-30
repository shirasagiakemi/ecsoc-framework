# ECSoC Methods — Detailed Documentation

> All methods are **exploratory**. No prospective validation has been performed.

---

## 1. Data Acquisition

All primary datasets are available from [PhysioNet](https://physionet.org/) (open access).

```bash
# Install WFDB tools for PhysioNet data access
pip install wfdb

# Or download directly:
# https://physionet.org/content/nsrdb/      — NSR Healthy
# https://physionet.org/content/chf2db/     — CHF2
# https://physionet.org/content/chfdb/      — BIDMC CHF
# https://physionet.org/content/cudb/       — CUDB
# https://physionet.org/content/svtdb/      — SVTDB
# https://physionet.org/content/mvedb/      — MVEDB
# https://physionet.org/content/cast-rr/    — CAST
```

---

## 2. RR Interval Extraction

```r
# From WFDB annotation files, extract RR intervals in milliseconds
# Recommended: use the wfdb R package or Python wfdb library

# Example (Python):
# import wfdb
# record = wfdb.rdann('nsrdb/16265', 'atr')
# rr_ms = np.diff(record.sample) / record.fs * 1000
```

**Preprocessing:**
- Remove ectopic beats (RR < 300 ms or > 2000 ms)
- Remove artifacts (consecutive RR change > 20%)
- Minimum series length: 256 beats (recommended: ≥ 512)

---

## 3. DFA Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Scale min | 4 beats | Minimum for reliable regression |
| Short-scale max | 16 beats | Standard α₁ upper limit |
| Long-scale max | 64 beats | Standard α₂ upper limit |
| Number of scales | 20 | Log-spaced |
| R² threshold (Phase V) | 0.93 | Derived from MUSIC cohort |

---

## 4. CHI Computation

```
CHI = 2 × (α₁ − α₂)
```

**Revised reference frame (ECSoC v2 — this paper):**

| Population | Mean CHI | % Supercritical |
|------------|----------|-----------------|
| NSR Healthy (N=54) | +0.440 ± 0.441 | 87% |
| CAST baseline (N=734) | +0.113 ± 0.461 | 65.7% |
| CHF2 NYHA 1–3 (N=29) | −0.243 ± 0.671 | 27.6% |

---

## 5. Phase Classification

```
Phase I    : α₁ < 0.85,  R² ≥ 0.93  → Subcritical stable
Phase II   : 0.85 ≤ α₁ < 1.15, R² ≥ 0.93 → Near-critical (critical manifold)
Phase III  : α₁ ≥ 1.15, R² ≥ 0.93  → Supercritical stable
Phase Va   : R² < 0.93, CHI < 0     → Structural breakdown, subcritical
Phase Vb   : R² < 0.93, CHI > 0     → Structural breakdown, supercritical
Phase Vc   : R² < 0.93, CHI ≈ 0     → Structural breakdown, near-critical
```

---

## 6. Layer 2 Trajectory Analysis

**P(III)_local(t):** Rolling proportion of windows with α₁ ≥ 1.15

**ΔP(III):** Change in P(III)_local from first 20% to last 20% of recording

| Stage | ΔP(III) | Prevalence (SVTDB) |
|-------|---------|---------------------|
| Late-rise (Stage 2) | ≥ +0.148 | 1.6–6.7% |
| Stable | \|ΔP(III)\| < 0.05 | ~60–66% |
| Suppress (Stage 1) | < −0.05 | 23–43% |
| Marginal+ | 0.05–0.148 | ~11% |

---

## 7. Statistical Methods

- Wilcoxon rank-sum / Mann-Whitney U for non-parametric group comparisons
- 95% CI for proportions: exact binomial (Clopper-Pearson)
- AUC reported for internal consistency only (no cross-validation)
- Spearman ρ for correlation with static HRV markers
