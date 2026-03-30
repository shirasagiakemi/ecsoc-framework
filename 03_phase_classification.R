# =============================================================================
# ECSoC Framework — Phase Classification (Layer 3)
# 03_phase_classification.R
#
# Assigns ECSoC phase (I–V) based on α₁ and DFA goodness-of-fit (R²).
# Phase V indicates structural scaling breakdown (R² < threshold).
# Phase V subtypes (Va/Vb/Vc) require CHI from 02_chi_calculation.R.
#
# Thresholds derived from MUSIC cohort (N=100) — exploratory only.
# Sensitivity analysis at R²= 0.90 / 0.93 / 0.95 in paper Section 6.
# =============================================================================

#' Classify ECSoC Phase
#'
#' @param alpha1       Short-scale DFA exponent (4–16 beats)
#' @param r_squared    DFA goodness-of-fit R² (full scale range)
#' @param chi          CHI value; required for Phase V subtype. Default: NULL
#' @param r2_threshold R² cutoff for structural breakdown. Default: 0.93
#'
#' @return List:
#'   $phase       — "I", "II", "III", "V"
#'   $subtype     — "Va", "Vb", "Vc", or NA
#'   $label       — Human-readable label (e.g. "Phase Vb")
#'   $description — Mechanistic interpretation
#'
#' @examples
#' classify_phase(1.3, 0.97, chi = 0.4)   # Phase III
#' classify_phase(1.2, 0.88, chi = 0.3)   # Phase Vb
classify_phase <- function(alpha1, r_squared, chi = NULL, r2_threshold = 0.93) {

  # ---- Phase V: Structural scaling breakdown --------------------------------
  if (r_squared < r2_threshold) {

    if (!is.null(chi)) {
      subtype <- dplyr::case_when(
        chi < -0.05 ~ "Va",
        chi >  0.05 ~ "Vb",
        TRUE        ~ "Vc"
      )
      description <- dplyr::case_when(
        subtype == "Va" ~ "Structural breakdown, subcritical: rigid low-variability dynamics with scaling loss",
        subtype == "Vb" ~ "Structural breakdown, supercritical: Path 2 collapse candidate (CHI > 0 with R² loss)",
        subtype == "Vc" ~ "Structural breakdown, near-critical manifold: transitional"
      )
    } else {
      subtype     <- "V"
      description <- "Structural breakdown (R² < threshold); CHI not supplied for subtype"
    }

    return(list(
      phase       = "V",
      subtype     = subtype,
      label       = paste0("Phase ", subtype),
      description = description,
      alpha1      = alpha1,
      r_squared   = r_squared,
      chi         = chi
    ))
  }

  # ---- Phases I–III: Stable scaling ----------------------------------------
  phase <- dplyr::case_when(
    alpha1 <  0.85 ~ "I",
    alpha1 <  1.15 ~ "II",
    TRUE           ~ "III"
  )

  description <- dplyr::case_when(
    phase == "I"   ~ "Subcritical stable: rigid, low-variability dynamics; suppressed autonomic modulation",
    phase == "II"  ~ "Near-critical stable: proximity to critical manifold (α₁ ≈ 1.0); optimal adaptability zone",
    phase == "III" ~ "Supercritical stable: dominant short-scale fluctuations; includes Metastable Supercritical State (healthy norm)"
  )

  list(
    phase       = phase,
    subtype     = NA_character_,
    label       = paste0("Phase ", phase),
    description = description,
    alpha1      = alpha1,
    r_squared   = r_squared,
    chi         = chi
  )
}


#' Classify phase with R² sensitivity analysis
#'
#' Re-runs classify_phase at three standard R² thresholds (0.90, 0.93, 0.95).
#' Used to assess robustness of Phase V calls (paper Section 6 / Supp. Fig. S2).
#'
#' @param alpha1   Short-scale DFA exponent
#' @param r_squared DFA R²
#' @param chi      CHI value
#'
#' @return Data frame with columns: threshold, phase, subtype, label
sensitivity_phase <- function(alpha1, r_squared, chi = NULL) {
  thresholds <- c(0.90, 0.93, 0.95)
  results <- lapply(thresholds, function(thr) {
    p <- classify_phase(alpha1, r_squared, chi, r2_threshold = thr)
    data.frame(
      threshold = thr,
      phase     = p$phase,
      subtype   = ifelse(is.na(p$subtype), "", p$subtype),
      label     = p$label,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}


#' Summarise phase distribution for a cohort
#'
#' @param df  Data frame with columns: alpha1, r_squared, chi (optional)
#' @param r2_threshold  R² threshold. Default: 0.93
#'
#' @return Data frame with phase counts and percentages
summarise_phases <- function(df, r2_threshold = 0.93) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' required.")

  chi_col <- if ("chi" %in% names(df)) df$chi else NULL

  phases <- mapply(
    function(a, r, c) classify_phase(a, r, c, r2_threshold)$label,
    df$alpha1, df$r_squared,
    if (is.null(chi_col)) rep(list(NULL), nrow(df)) else chi_col,
    SIMPLIFY = TRUE
  )

  phase_counts <- table(phases)
  data.frame(
    phase   = names(phase_counts),
    n       = as.integer(phase_counts),
    percent = round(as.numeric(phase_counts) / nrow(df) * 100, 1),
    stringsAsFactors = FALSE
  )
}
