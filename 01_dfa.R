# =============================================================================
# ECSoC Framework — DFA Implementation
# 01_dfa.R
#
# Detrended Fluctuation Analysis (DFA) for RR interval time series.
# Computes short-scale (α₁) and long-scale (α₂) scaling exponents,
# along with DFA goodness-of-fit (R²) used as a proxy for manifold integrity.
#
# Reference: Peng CK, et al. Chaos. 1995;5(1):82–87.
#            Chen Z, et al. Phys Rev E. 2002;65:041107.
# =============================================================================

#' Detrended Fluctuation Analysis
#'
#' @param rr         Numeric vector of RR intervals (milliseconds)
#' @param scale_min  Minimum box size (beats). Default: 4
#' @param scale_max_short  Upper limit of short-scale regime (beats). Default: 16
#' @param scale_max_long   Upper limit of long-scale regime (beats). Default: 64
#' @param n_scales   Number of logarithmically spaced scales. Default: 20
#'
#' @return List with:
#'   alpha1     — Short-scale DFA exponent (4–16 beats)
#'   alpha2     — Long-scale DFA exponent (16–64 beats)
#'   r_squared  — Goodness-of-fit R² of the log-log DFA scaling (full range)
#'   scales     — All box sizes used
#'   F_n        — Corresponding DFA fluctuation values
#'   fit_short  — Linear fit object for short-scale regime
#'   fit_long   — Linear fit object for long-scale regime
#'
#' @examples
#' rr <- read.csv("data/sample/sample_rr_normal.csv")$rr_ms
#' result <- compute_dfa(rr)
#' cat("alpha1:", result$alpha1, "alpha2:", result$alpha2, "R2:", result$r_squared)
compute_dfa <- function(rr,
                        scale_min       = 4,
                        scale_max_short = 16,
                        scale_max_long  = 64,
                        n_scales        = 20) {

  # --- Input validation ---
  if (!is.numeric(rr) || length(rr) < scale_max_long * 4) {
    stop(sprintf(
      "RR series too short. Need at least %d beats; got %d.",
      scale_max_long * 4, length(rr)
    ))
  }
  if (any(is.na(rr))) {
    warning("NA values detected in RR series. Removing before DFA.")
    rr <- rr[!is.na(rr)]
  }
  if (any(rr <= 0)) {
    stop("RR intervals must be positive.")
  }

  N <- length(rr)

  # --- Step 1: Integrate (cumulative sum of mean-subtracted series) ---
  y <- cumsum(rr - mean(rr))

  # --- Step 2: Define logarithmically spaced box sizes ---
  scales <- unique(round(exp(seq(
    log(scale_min),
    log(scale_max_long),
    length.out = n_scales
  ))))
  scales <- scales[scales >= scale_min & scales <= scale_max_long]
  scales <- scales[scales <= floor(N / 4)]  # each box needs at least 4 points

  if (length(scales) < 4) {
    stop("Too few valid scales. Consider a longer RR series.")
  }

  # --- Step 3: Compute DFA fluctuation F(n) for each box size ---
  F_n <- sapply(scales, function(s) {
    n_boxes <- floor(N / s)
    residuals_sq <- numeric(n_boxes)
    for (b in seq_len(n_boxes)) {
      idx   <- ((b - 1) * s + 1):(b * s)
      t_box <- seq_along(idx)
      fit   <- lm(y[idx] ~ t_box)
      residuals_sq[b] <- mean(residuals(fit)^2)
    }
    sqrt(mean(residuals_sq))
  })

  # --- Step 4: Log-log regression for α₁ and α₂ ---
  log_s  <- log10(scales)
  log_Fn <- log10(F_n)

  # Short-scale regime: scale_min to scale_max_short
  idx_short <- scales >= scale_min & scales <= scale_max_short
  # Long-scale regime: scale_max_short to scale_max_long
  idx_long  <- scales >= scale_max_short & scales <= scale_max_long

  if (sum(idx_short) < 2 || sum(idx_long) < 2) {
    stop("Insufficient scales in one regime. Adjust scale parameters.")
  }

  fit_short <- lm(log_Fn[idx_short] ~ log_s[idx_short])
  fit_long  <- lm(log_Fn[idx_long]  ~ log_s[idx_long])

  alpha1 <- coef(fit_short)[2]
  alpha2 <- coef(fit_long)[2]

  # --- Step 5: Global R² over the full scale range ---
  fit_global <- lm(log_Fn ~ log_s)
  ss_res     <- sum(residuals(fit_global)^2)
  ss_tot     <- sum((log_Fn - mean(log_Fn))^2)
  r_squared  <- 1 - ss_res / ss_tot

  # --- Return results ---
  list(
    alpha1    = unname(alpha1),
    alpha2    = unname(alpha2),
    r_squared = r_squared,
    scales    = scales,
    F_n       = F_n,
    fit_short = fit_short,
    fit_long  = fit_long,
    N         = N
  )
}


#' Plot DFA log-log scaling plot
#'
#' @param dfa_result  Output from compute_dfa()
#' @param title       Plot title
#' @param save_path   If provided, saves plot to this path
plot_dfa <- function(dfa_result, title = "DFA Scaling Plot", save_path = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
  }
  library(ggplot2)

  df <- data.frame(
    log_scale = log10(dfa_result$scales),
    log_Fn    = log10(dfa_result$F_n)
  )

  alpha1 <- round(dfa_result$alpha1, 3)
  alpha2 <- round(dfa_result$alpha2, 3)
  r2     <- round(dfa_result$r_squared, 4)

  p <- ggplot(df, aes(x = log_scale, y = log_Fn)) +
    geom_point(size = 2, color = "#2C5F8A") +
    geom_smooth(method = "lm", se = FALSE, color = "#E05A3A", linewidth = 0.8) +
    annotate("text", x = min(df$log_scale) + 0.05, y = max(df$log_Fn),
             label = sprintf("α₁ = %.3f | α₂ = %.3f | R² = %.4f", alpha1, alpha2, r2),
             hjust = 0, size = 3.5, color = "gray30") +
    labs(title = title,
         x = "log₁₀(Box size n)",
         y = "log₁₀(F(n))") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 7, height = 5, dpi = 150)
    message("Plot saved to: ", save_path)
  }

  p
}
