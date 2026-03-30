# =============================================================================
# ECSoC Framework — Multi-Cohort Analysis
# 05_cohort_analysis.R
#
# Functions for:
#   - Running ECSoC pipeline across a full cohort
#   - Generating Table 9 (cross-cohort Phase V temporal gradient)
#   - Computing CHI group comparisons (Mann-Whitney U)
#   - NYHA dose-response gradient analysis
#   - Va/Vb directional alignment verification
#
# All results exploratory. No cross-validation performed.
# =============================================================================

library(dplyr)

# Source dependencies
source("R/01_dfa.R")
source("R/02_chi_calculation.R")
source("R/03_phase_classification.R")
source("R/04_trajectory_analysis.R")


# =============================================================================
# 1. Batch ECSoC Analysis
# =============================================================================

#' Run full ECSoC pipeline on a list of RR series
#'
#' @param rr_list    Named list of numeric vectors (RR intervals in ms)
#' @param r2_threshold  Phase V threshold. Default: 0.93
#' @param verbose    Print progress. Default: TRUE
#'
#' @return Data frame with one row per subject:
#'   id, alpha1, alpha2, r_squared, chi, chi_domain, phase, label, n_beats
run_cohort <- function(rr_list, r2_threshold = 0.93, verbose = TRUE) {

  results <- lapply(names(rr_list), function(id) {
    if (verbose) message("  Processing: ", id)

    rr <- rr_list[[id]]

    tryCatch({
      dfa   <- compute_dfa(rr)
      chi   <- compute_chi(dfa$alpha1, dfa$alpha2)
      phase <- classify_phase(dfa$alpha1, dfa$r_squared, chi, r2_threshold)

      data.frame(
        id          = id,
        n_beats     = length(rr),
        alpha1      = dfa$alpha1,
        alpha2      = dfa$alpha2,
        r_squared   = dfa$r_squared,
        chi         = chi,
        chi_domain  = classify_chi_domain(chi),
        phase       = phase$phase,
        subtype     = ifelse(is.na(phase$subtype), "", phase$subtype),
        label       = phase$label,
        phase_v     = phase$phase == "V",
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      warning("Failed on subject ", id, ": ", conditionMessage(e))
      data.frame(
        id = id, n_beats = length(rr),
        alpha1 = NA, alpha2 = NA, r_squared = NA,
        chi = NA, chi_domain = NA, phase = NA,
        subtype = NA, label = NA, phase_v = NA,
        stringsAsFactors = FALSE
      )
    })
  })

  do.call(rbind, results)
}


# =============================================================================
# 2. Cohort Summary Statistics
# =============================================================================

#' Summarise a cohort data frame (output of run_cohort)
#'
#' @param df       Data frame from run_cohort()
#' @param cohort_name  Label for this cohort
#' @param r2_threshold  Phase V threshold used. Default: 0.93
#'
#' @return Named list of summary statistics
summarise_cohort <- function(df, cohort_name = "Cohort", r2_threshold = 0.93) {

  df_valid <- df[!is.na(df$alpha1), ]
  n_total  <- nrow(df)
  n_valid  <- nrow(df_valid)

  phase_v_n   <- sum(df_valid$phase_v, na.rm = TRUE)
  phase_v_pct <- phase_v_n / n_valid * 100

  # Exact binomial 95% CI for Phase V proportion
  ci <- binom.test(phase_v_n, n_valid)$conf.int * 100

  super_n   <- sum(df_valid$chi > 0, na.rm = TRUE)
  super_pct <- super_n / n_valid * 100

  list(
    cohort        = cohort_name,
    n_total       = n_total,
    n_valid       = n_valid,
    alpha1_mean   = mean(df_valid$alpha1, na.rm = TRUE),
    alpha1_sd     = sd(df_valid$alpha1,   na.rm = TRUE),
    alpha2_mean   = mean(df_valid$alpha2, na.rm = TRUE),
    r2_mean       = mean(df_valid$r_squared, na.rm = TRUE),
    chi_mean      = mean(df_valid$chi, na.rm = TRUE),
    chi_sd        = sd(df_valid$chi,   na.rm = TRUE),
    chi_median    = median(df_valid$chi, na.rm = TRUE),
    pct_super     = super_pct,
    phase_v_n     = phase_v_n,
    phase_v_pct   = phase_v_pct,
    phase_v_ci_lo = ci[1],
    phase_v_ci_hi = ci[2],
    r2_threshold  = r2_threshold
  )
}


#' Print cohort summary as formatted table
print_cohort_summary <- function(summary_list) {
  cat(sprintf("\n=== %s (N=%d valid / %d total) ===\n",
              summary_list$cohort, summary_list$n_valid, summary_list$n_total))
  cat(sprintf("  α₁            : %.3f ± %.3f\n", summary_list$alpha1_mean, summary_list$alpha1_sd))
  cat(sprintf("  α₂            : %.3f\n", summary_list$alpha2_mean))
  cat(sprintf("  R² (mean)     : %.4f\n", summary_list$r2_mean))
  cat(sprintf("  CHI           : %.3f ± %.3f (median %.3f)\n",
              summary_list$chi_mean, summary_list$chi_sd, summary_list$chi_median))
  cat(sprintf("  %% Supercritical: %.1f%%\n", summary_list$pct_super))
  cat(sprintf("  Phase V       : %d/%d (%.1f%%; 95%% CI: %.1f–%.1f%%)\n",
              summary_list$phase_v_n, summary_list$n_valid,
              summary_list$phase_v_pct,
              summary_list$phase_v_ci_lo, summary_list$phase_v_ci_hi))
  cat(sprintf("  [R² threshold : %.2f; all values exploratory]\n", summary_list$r2_threshold))
}


# =============================================================================
# 3. Cross-Cohort CHI Comparison (Mann-Whitney U)
# =============================================================================

#' Compare CHI between two cohorts
#'
#' @param df1, df2  Data frames from run_cohort()
#' @param name1, name2  Cohort labels
#'
#' @return List with statistic, p-value, and effect size (rank-biserial r)
compare_chi <- function(df1, df2, name1 = "Group 1", name2 = "Group 2") {

  chi1 <- df1$chi[!is.na(df1$chi)]
  chi2 <- df2$chi[!is.na(df2$chi)]

  test <- wilcox.test(chi1, chi2, exact = FALSE)

  # Rank-biserial correlation (effect size)
  n1 <- length(chi1); n2 <- length(chi2)
  r_rb <- 1 - 2 * test$statistic / (n1 * n2)

  cat(sprintf("\n--- CHI Comparison: %s vs %s ---\n", name1, name2))
  cat(sprintf("  %s: mean=%.3f ± %.3f, median=%.3f (N=%d)\n",
              name1, mean(chi1), sd(chi1), median(chi1), n1))
  cat(sprintf("  %s: mean=%.3f ± %.3f, median=%.3f (N=%d)\n",
              name2, mean(chi2), sd(chi2), median(chi2), n2))
  cat(sprintf("  Mann-Whitney U = %.0f, p = %.4g\n", test$statistic, test$p.value))
  cat(sprintf("  Rank-biserial r = %.3f\n", r_rb))
  cat("  [Exploratory; no correction for multiple comparisons]\n")

  list(
    statistic   = test$statistic,
    p_value     = test$p.value,
    effect_size = r_rb,
    n1 = n1, n2 = n2,
    mean1 = mean(chi1), mean2 = mean(chi2)
  )
}


# =============================================================================
# 4. Cross-Cohort Phase V Gradient (Table 9 equivalent)
# =============================================================================

#' Build cross-cohort Phase V temporal gradient table
#'
#' Replicates Table 9 from the paper using a list of cohort summaries.
#'
#' @param summary_list  List of outputs from summarise_cohort()
#' @return Data frame formatted as Table 9
build_gradient_table <- function(summary_list) {
  rows <- lapply(summary_list, function(s) {
    data.frame(
      Cohort       = s$cohort,
      N            = s$n_valid,
      PhaseV_n     = s$phase_v_n,
      PhaseV_pct   = round(s$phase_v_pct, 1),
      CI_lo        = round(s$phase_v_ci_lo, 1),
      CI_hi        = round(s$phase_v_ci_hi, 1),
      CHI_mean     = round(s$chi_mean, 3),
      Pct_super    = round(s$pct_super, 1),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  df <- df[order(df$PhaseV_pct), ]
  rownames(df) <- NULL
  df
}


# =============================================================================
# 5. Va/Vb Directional Alignment Check
# =============================================================================

#' Verify Va/Vb CHI directional alignment
#'
#' Checks that all Phase Va cases have CHI < 0 and all Phase Vb cases have CHI > 0.
#' Paper reports 115/115 = 100% concordance.
#'
#' @param df  Data frame from run_cohort() — Phase V subjects only
#' @return List with concordance count and percentage
check_vavb_alignment <- function(df) {
  pv <- df[df$phase == "V" & !is.na(df$chi), ]

  if (nrow(pv) == 0) {
    message("No Phase V cases found.")
    return(NULL)
  }

  va_correct <- sum(pv$subtype == "Va" & pv$chi < 0, na.rm = TRUE)
  vb_correct <- sum(pv$subtype == "Vb" & pv$chi > 0, na.rm = TRUE)
  total_pv   <- nrow(pv)
  concordant <- va_correct + vb_correct
  pct        <- concordant / total_pv * 100
  ci         <- binom.test(concordant, total_pv)$conf.int * 100

  cat(sprintf("\n--- Va/Vb Directional Alignment ---\n"))
  cat(sprintf("  Phase V total : %d\n", total_pv))
  cat(sprintf("  Va (CHI < 0)  : %d correct\n", va_correct))
  cat(sprintf("  Vb (CHI > 0)  : %d correct\n", vb_correct))
  cat(sprintf("  Concordance   : %d/%d = %.1f%% (95%% CI: %.1f–%.1f%%)\n",
              concordant, total_pv, pct, ci[1], ci[2]))
  cat("  [Paper reports 115/115 = 100% across five datasets]\n")

  list(concordant = concordant, total = total_pv,
       pct = pct, ci_lo = ci[1], ci_hi = ci[2])
}


# =============================================================================
# 6. CAST Drug Effect Analysis (ΔCHI)
# =============================================================================

#' Paired CHI analysis: pre vs post antiarrhythmic drug (CAST protocol)
#'
#' @param df_pre   Data frame (run_cohort output) for baseline recordings
#' @param df_post  Data frame (run_cohort output) for post-drug recordings
#'        Must have matching 'id' column.
#'
#' @return List with ΔCHI statistics and Wilcoxon test result
analyse_drug_effect <- function(df_pre, df_post) {
  merged <- merge(df_pre[, c("id","chi","r_squared")],
                  df_post[, c("id","chi","r_squared")],
                  by = "id", suffixes = c("_pre","_post"))
  merged <- merged[complete.cases(merged), ]

  delta_chi <- merged$chi_post - merged$chi_pre
  delta_r2  <- merged$r_squared_post - merged$r_squared_pre

  test_chi <- wilcox.test(delta_chi, mu = 0, exact = FALSE)
  test_r2  <- wilcox.test(delta_r2,  mu = 0, exact = FALSE)

  cat(sprintf("\n--- CAST Drug Effect (N=%d paired subjects) ---\n", nrow(merged)))
  cat(sprintf("  ΔCHI = %.3f ± %.3f  (Wilcoxon p = %.4g)\n",
              mean(delta_chi), sd(delta_chi), test_chi$p.value))
  cat(sprintf("  ΔR²  = %.4f ± %.4f  (Wilcoxon p = %.4f)\n",
              mean(delta_r2), sd(delta_r2), test_r2$p.value))
  cat("  [Paper: ΔCHI = -0.138, p < 0.0001; ΔR² NS, p = 0.592]\n")
  cat("  [Supports Layer 1/Layer 3 independence — exploratory]\n")

  list(
    n          = nrow(merged),
    delta_chi  = mean(delta_chi),
    delta_chi_sd = sd(delta_chi),
    p_chi      = test_chi$p.value,
    delta_r2   = mean(delta_r2),
    p_r2       = test_r2$p.value
  )
}
