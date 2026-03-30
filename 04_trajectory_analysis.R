# =============================================================================
# ECSoC Framework — Layer 2 Trajectory Analysis
# 04_trajectory_analysis.R
#
# Computes ΔP(III) — the non-stationary trajectory of supercritical dynamics
# — and detects Stage 1 (Suppress) and Stage 2 (Late-rise) patterns.
#
# All results are exploratory. See paper Section 2.2.
# =============================================================================

#' Compute P(III)_local(t) — rolling proportion of supercritical windows
#'
#' Divides the RR series into overlapping windows, computes α₁ for each,
#' and tracks the proportion of windows with α₁ ≥ 1.15 (Phase III territory).
#'
#' @param rr          Numeric vector of RR intervals (ms)
#' @param window_size Number of beats per window. Default: 200
#' @param step_size   Step between windows (beats). Default: 50
#' @param alpha1_threshold  α₁ threshold for supercritical classification. Default: 1.15
#'
#' @return Data frame with columns: time_index, p_local, alpha1_local
compute_p_local <- function(rr,
                            window_size     = 200,
                            step_size       = 50,
                            alpha1_threshold = 1.15) {

  source("R/01_dfa.R", local = TRUE)

  N       <- length(rr)
  starts  <- seq(1, N - window_size + 1, by = step_size)

  results <- lapply(starts, function(s) {
    window <- rr[s:(s + window_size - 1)]
    tryCatch({
      dfa <- compute_dfa(window)
      data.frame(
        time_index  = s + window_size / 2,
        p_local     = as.numeric(dfa$alpha1 >= alpha1_threshold),
        alpha1_local = dfa$alpha1
      )
    }, error = function(e) {
      data.frame(time_index = s + window_size / 2, p_local = NA, alpha1_local = NA)
    })
  })

  df <- do.call(rbind, results)

  # Rolling mean of p_local to smooth
  k <- max(3, round(300 / step_size))
  df$p_local_smooth <- stats::filter(df$p_local, rep(1/k, k), sides = 2)

  df
}


#' Compute ΔP(III) — global trajectory metric
#'
#' ΔP(III) = P(III)_local(end) − P(III)_local(start)
#' Captures whether supercritical proportion rises or falls over the recording.
#'
#' @param p_local_df  Output from compute_p_local()
#' @param frac        Fraction of recording used for start/end windows. Default: 0.2
#'
#' @return List with delta_p3, stage, and interpretation
compute_delta_p3 <- function(p_local_df, frac = 0.2) {

  n       <- nrow(p_local_df)
  n_frac  <- max(1, round(n * frac))

  p_start <- mean(p_local_df$p_local[1:n_frac], na.rm = TRUE)
  p_end   <- mean(p_local_df$p_local[(n - n_frac + 1):n], na.rm = TRUE)

  delta_p3 <- p_end - p_start

  # Stage classification (thresholds from MUSIC cohort derivation — exploratory)
  stage <- dplyr::case_when(
    delta_p3 >=  0.148        ~ "Late-rise",   # Stage 2 candidate
    delta_p3 < -0.05          ~ "Suppress",    # Stage 1
    abs(delta_p3) < 0.05      ~ "Stable",
    TRUE                       ~ "Marginal+"
  )

  interpretation <- dplyr::case_when(
    stage == "Late-rise"  ~ "Stage 2 (Late-rise): putative SCD-precursor trajectory [RARE: ~2-7% of events]",
    stage == "Suppress"   ~ "Stage 1 (Suppress): subcritical shift or post-event withdrawal",
    stage == "Stable"     ~ "Stable: ergodic Layer 1 dominance; no trajectory activation",
    stage == "Marginal+"  ~ "Marginal+: near-critical boundary; monitoring may be warranted"
  )

  full_pattern <- stage == "Late-rise"  # requires prior Suppress in multi-episode analysis

  list(
    delta_p3      = delta_p3,
    p_start       = p_start,
    p_end         = p_end,
    stage         = stage,
    interpretation = interpretation,
    full_pattern_candidate = full_pattern
  )
}


#' Full Layer 2 trajectory analysis
#'
#' @param rr  Numeric vector of RR intervals (ms)
#' @param ...  Arguments passed to compute_p_local()
#' @return List with trajectory data and ΔP(III) result
analyze_trajectory <- function(rr, ...) {
  p_local <- compute_p_local(rr, ...)
  delta   <- compute_delta_p3(p_local)

  list(
    p_local_timeseries = p_local,
    delta_p3_result    = delta,
    stage              = delta$stage,
    interpretation     = delta$interpretation
  )
}


#' Plot P(III)_local(t) trajectory
#'
#' @param traj_result  Output from analyze_trajectory()
#' @param title        Plot title
plot_trajectory <- function(traj_result, title = "P(III)_local(t) Trajectory") {
  library(ggplot2)

  df    <- traj_result$p_local_timeseries
  stage <- traj_result$stage
  dp3   <- round(traj_result$delta_p3_result$delta_p3, 3)

  df_plot <- df[!is.na(df$p_local_smooth), ]

  ggplot(df_plot, aes(x = time_index)) +
    geom_line(aes(y = p_local), color = "lightblue", alpha = 0.5, linewidth = 0.4) +
    geom_line(aes(y = p_local_smooth), color = "#2C5F8A", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    annotate("text", x = max(df_plot$time_index),
             y = max(df_plot$p_local_smooth, na.rm = TRUE),
             label = sprintf("ΔP(III) = %.3f\nStage: %s", dp3, stage),
             hjust = 1, vjust = 1, size = 3.5, color = "#E05A3A") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(title = title,
         subtitle = "[Exploratory — Layer 2 trajectory analysis]",
         x = "Beat index (window midpoint)",
         y = "P(III)_local") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(color = "gray50", size = 9))
}
