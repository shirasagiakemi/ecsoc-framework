# =============================================================================
# ECSoC Framework — CHI Calculation & Phase Classification
# 02_chi_calculation.R  /  03_phase_classification.R
# =============================================================================

# --- CHI Calculation ---------------------------------------------------------

#' Compute Criticality Heterogeneity Index (CHI)
#'
#' CHI = 2(α₁ − α₂)
#'   CHI > 0  → Supercritical (dominant short-scale fluctuations)
#'   CHI ≈ 0  → Critical manifold
#'   CHI < 0  → Subcritical (suppressed short-scale dynamics)
#'
#' Revised reference frame (from NSR Healthy Database, N=54):
#'   Healthy norm: mean CHI ≈ +0.440 (87% supercritical)
#'   CHF (NYHA 1–3): mean CHI ≈ −0.243 (28% supercritical)
#'
#' @param alpha1  Short-scale DFA exponent
#' @param alpha2  Long-scale DFA exponent
#' @return CHI value (numeric)
compute_chi <- function(alpha1, alpha2) {
  2 * (alpha1 - alpha2)
}

#' Classify CHI domain
#'
#' @param chi  CHI value
#' @return Character: "supercritical", "near-critical", or "subcritical"
classify_chi_domain <- function(chi) {
  dplyr::case_when(
    chi >  0.05  ~ "supercritical",
    chi < -0.05  ~ "subcritical",
    TRUE         ~ "near-critical"
  )
}

#' Interpret CHI relative to healthy reference (Revised ECSoC CHI frame)
#'
#' @param chi  CHI value
#' @return List with domain, deviation from healthy mean, and interpretation
interpret_chi <- function(chi) {
  healthy_mean <- 0.440   # NSR Healthy Database (N=54)
  chf_mean     <- -0.243  # CHF2 NYHA 1–3 (N=29)

  deviation <- chi - healthy_mean
  domain    <- classify_chi_domain(chi)

  interpretation <- dplyr::case_when(
    chi >= healthy_mean - 0.441        ~ "Within healthy supercritical range",
    chi >= 0                           ~ "Supercritical but below healthy norm",
    chi >= chf_mean                    ~ "Subcritical shift; consistent with early-moderate CHF",
    TRUE                               ~ "Markedly subcritical; consistent with advanced disease"
  )

  list(
    chi            = chi,
    domain         = domain,
    deviation_from_healthy = deviation,
    interpretation = interpretation
  )
}


# =============================================================================
# Phase Classification (Layer 3: structural breakdown via R²)
# =============================================================================

#' Classify ECSoC Phase based on α₁ and R²
#'
#' Phase definitions (exploratory; thresholds derived from MUSIC cohort):
#'
#'   Phase I   — Subcritical stable    : α₁ < 0.85, R² ≥ 0.93
#'   Phase II  — Near-critical stable  : 0.85 ≤ α₁ < 1.15, R² ≥ 0.93
#'   Phase III — Supercritical stable  : α₁ ≥ 1.15, R² ≥ 0.93
#'   Phase IV  — Mixed / transitional  : intermediate R²
#'   Phase Va  — Structural breakdown, subcritical  : R² < 0.93, CHI < 0
#'   Phase Vb  — Structural breakdown, supercritical: R² < 0.93, CHI > 0
#'   Phase Vc  — Structural breakdown, near-critical: R² < 0.93, CHI ≈ 0
#'
#' @param alpha1      Short-scale DFA exponent
#' @param r_squared   DFA goodness-of-fit R²
#' @param chi         CHI value (required for Phase V subtype)
#' @param r2_threshold  R² threshold for Phase V. Default: 0.93
#'
#' @return List with phase label, subtype, and description
classify_phase <- function(alpha1, r_squared, chi = NULL, r2_threshold = 0.93) {

  # --- Phase V: Structural breakdown ---
  if (r_squared < r2_threshold) {
    if (!is.null(chi)) {
      subtype <- dplyr::case_when(
        chi < -0.05  ~ "Va",   # subcritical + breakdown
        chi >  0.05  ~ "Vb",   # supercritical + breakdown
        TRUE         ~ "Vc"    # near-critical + breakdown
      )
      desc <- dplyr::case_when(
        subtype == "Va" ~ "Structural breakdown, subcritical (suppressed dynamics)",
        subtype == "Vb" ~ "Structural breakdown, supercritical (Path 2 candidate)",
        subtype == "Vc" ~ "Structural breakdown, near-critical manifold"
      )
    } else {
      subtype <- "V (subtype unknown)"
      desc    <- "Structural breakdown (CHI not provided for subtype)"
    }
    return(list(
      phase       = "V",
      subtype     = subtype,
      label       = paste0("Phase ", subtype),
      description = desc,
      r_squared   = r_squared,
      alpha1      = alpha1
    ))
  }

  # --- Phases I–III: Stable scaling ---
  phase <- dplyr::case_when(
    alpha1 < 0.85  ~ "I",
    alpha1 < 1.15  ~ "II",
    TRUE           ~ "III"
  )
  desc <- dplyr::case_when(
    phase == "I"   ~ "Subcritical: rigid, low-variability dynamics",
    phase == "II"  ~ "Near-critical: proximity to critical manifold (α₁ ≈ 1.0)",
    phase == "III" ~ "Supercritical: dominant short-scale fluctuations"
  )

  list(
    phase       = phase,
    subtype     = NA,
    label       = paste0("Phase ", phase),
    description = desc,
    r_squared   = r_squared,
    alpha1      = alpha1
  )
}


#' Run full ECSoC classification on a single RR series
#'
#' Convenience wrapper: DFA → CHI → Phase
#'
#' @param rr  Numeric vector of RR intervals (ms)
#' @param ...  Additional arguments passed to compute_dfa()
#' @return List with all ECSoC metrics
run_ecsoc <- function(rr, ...) {
  source("R/01_dfa.R", local = TRUE)

  dfa    <- compute_dfa(rr, ...)
  chi    <- compute_chi(dfa$alpha1, dfa$alpha2)
  phase  <- classify_phase(dfa$alpha1, dfa$r_squared, chi)
  interp <- interpret_chi(chi)

  list(
    # DFA
    alpha1    = dfa$alpha1,
    alpha2    = dfa$alpha2,
    r_squared = dfa$r_squared,
    # CHI
    chi       = chi,
    chi_domain      = interp$domain,
    chi_interpretation = interp$interpretation,
    # Phase
    phase     = phase$phase,
    phase_label = phase$label,
    phase_description = phase$description,
    # Raw
    dfa_result = dfa
  )
}
