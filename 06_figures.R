# =============================================================================
# ECSoC Framework — Figure Generation
# 06_figures.R
#
# Reproduces key figures from the paper using ggplot2.
# All figures are exploratory; data distributions simulated from
# reported cohort statistics (mean ± SD) where originals are unavailable.
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# --- Helper: Create output directory -----------------------------------------
ensure_fig_dir <- function(path = "figures") {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  path
}


# =============================================================================
# Figure 1: Phase Space Scatter Plot (α₁ × R²)
# =============================================================================

#' Plot ECSoC Phase Space
#'
#' @param df  Data frame with columns: alpha1, r_squared, group (factor), chi
#' @param r2_threshold  R² threshold for Phase V line. Default: 0.93
plot_phase_space <- function(df, r2_threshold = 0.93, save = FALSE) {
  fig_dir <- ensure_fig_dir()

  # Assign phase colors
  df <- df %>%
    mutate(
      phase_v = r_squared < r2_threshold,
      fill_col = case_when(
        phase_v & chi >  0.05  ~ "#E05A3A",  # Vb — supercritical breakdown
        phase_v & chi < -0.05  ~ "#8B2FC9",  # Va — subcritical breakdown
        phase_v                ~ "#E0A33A",  # Vc
        chi > 0                ~ "#2C8AF5",  # supercritical stable
        TRUE                   ~ "#4CAF50"   # subcritical/near-critical stable
      )
    )

  p <- ggplot(df, aes(x = alpha1, y = r_squared, color = fill_col, shape = group)) +
    geom_hline(yintercept = r2_threshold, linetype = "dashed",
               color = "red", linewidth = 0.7, alpha = 0.7) +
    geom_vline(xintercept = 1.0, linetype = "dotted",
               color = "gray40", linewidth = 0.5) +
    geom_point(size = 2.5, alpha = 0.75) +
    scale_color_identity() +
    annotate("text", x = 1.0, y = r2_threshold + 0.003,
             label = sprintf("Phase V threshold (R² = %.2f)", r2_threshold),
             hjust = -0.05, size = 3, color = "red") +
    annotate("text", x = 1.0, y = min(df$r_squared, na.rm = TRUE),
             label = "α₁ = 1.0\n(critical point)",
             hjust = -0.1, vjust = 0, size = 2.8, color = "gray40") +
    labs(
      title    = "ECSoC Phase Space: α₁ × R²",
      subtitle = "[Exploratory] Phase V (R² < 0.93) = structural scaling breakdown",
      x        = "DFA Short-scale Exponent (α₁)",
      y        = "DFA Goodness-of-Fit (R²)",
      shape    = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray50", size = 9),
      legend.position = "bottom"
    )

  if (save) {
    out <- file.path(fig_dir, "fig1_phase_space.png")
    ggsave(out, p, width = 8, height = 6, dpi = 150)
    message("Saved: ", out)
  }
  p
}


# =============================================================================
# Figure 2: CHI Distribution Across Cohorts (Violin Plot)
# Reproduces Figure 8 from the paper
# =============================================================================

#' CHI Distribution Violin Plot (Figure 8 equivalent)
#'
#' Uses reported cohort statistics (mean ± SD) to simulate distributions.
#' NOTE: Actual individual-level data from MUSIC/SVTDB are not included here.
plot_chi_distribution <- function(save = FALSE) {
  fig_dir <- ensure_fig_dir()

  set.seed(42)

  # Reported cohort statistics (Table 13 and related tables in paper)
  cohorts <- list(
    list(name = "NSR Healthy\n(N=54)",     mean_chi =  0.440, sd_chi = 0.441, n = 54,   nyha = "Healthy"),
    list(name = "CAST Baseline\n(N=734)",  mean_chi =  0.113, sd_chi = 0.461, n = 734,  nyha = "Ectopy"),
    list(name = "CHF2 NYHA 1-3\n(N=29)",  mean_chi = -0.243, sd_chi = 0.671, n = 29,   nyha = "CHF mild"),
    list(name = "MUSIC SCD\n(N=30)",       mean_chi = -0.350, sd_chi = 0.750, n = 30,   nyha = "SCD"),
    list(name = "MUSIC PFD\n(N=38)",       mean_chi = -0.450, sd_chi = 0.800, n = 38,   nyha = "PFD"),
    list(name = "SVTDB VT/VF\n(N=76)",    mean_chi = -0.530, sd_chi = 0.920, n = 76,   nyha = "Arrhythmia")
  )

  df <- do.call(rbind, lapply(cohorts, function(co) {
    data.frame(
      cohort = co$name,
      chi    = rnorm(co$n, co$mean_chi, co$sd_chi),
      group  = co$nyha,
      stringsAsFactors = FALSE
    )
  }))

  df$cohort <- factor(df$cohort, levels = sapply(cohorts, `[[`, "name"))

  # Percent supercritical per cohort
  pct_super <- df %>%
    group_by(cohort) %>%
    summarise(pct = mean(chi > 0) * 100, .groups = "drop") %>%
    mutate(label = sprintf("%.0f%%", pct))

  p <- ggplot(df, aes(x = cohort, y = chi, fill = cohort)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +
    geom_violin(alpha = 0.7, trim = TRUE, width = 0.8) +
    geom_boxplot(width = 0.15, outlier.size = 0.8, fill = "white", alpha = 0.9) +
    geom_text(data = pct_super,
              aes(x = cohort, y = 2.2, label = label),
              size = 3.2, fontface = "bold", inherit.aes = FALSE) +
    scale_fill_brewer(palette = "Blues", guide = "none") +
    annotate("text", x = 0.55, y = 2.3, label = "% supercritical (CHI > 0)",
             hjust = 0, size = 3, color = "gray40") +
    labs(
      title    = "CHI Distribution Across Disease Severity",
      subtitle = "[Exploratory] Distributions simulated from reported mean ± SD; original N values used",
      x        = NULL,
      y        = "CHI = 2(α₁ − α₂)"
    ) +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray50", size = 8),
      axis.text.x   = element_text(size = 8)
    )

  if (save) {
    out <- file.path(fig_dir, "fig8_chi_distribution.png")
    ggsave(out, p, width = 10, height = 6, dpi = 150)
    message("Saved: ", out)
  }
  p
}


# =============================================================================
# Figure 3: NYHA Phase V Dose-Response Gradient
# =============================================================================

plot_nyha_gradient <- function(save = FALSE) {
  fig_dir <- ensure_fig_dir()

  df <- data.frame(
    cohort   = c("NSR Healthy\n(N=54)", "CHF2 NYHA 1-3\n(N=29)",
                 "BIDMC NYHA 3-4\n(N=15)", "MVEDB\n(N=22)"),
    phase_v  = c(0, 0, 53.3, 90.9),
    ci_lower = c(0, 0, 26.6, 69.9),
    ci_upper = c(6.6, 11.7, 78.7, 97.2),
    severity = 1:4
  )
  df$cohort <- factor(df$cohort, levels = df$cohort)

  p <- ggplot(df, aes(x = cohort, y = phase_v)) +
    geom_col(fill = "#2C5F8A", alpha = 0.8, width = 0.6) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  width = 0.2, color = "#E05A3A", linewidth = 0.8) +
    geom_text(aes(label = ifelse(phase_v == 0,
                                  "0%",
                                  sprintf("%.1f%%", phase_v))),
              vjust = -0.8, size = 4, fontface = "bold", color = "#1A3A5C") +
    scale_y_continuous(limits = c(0, 105), labels = function(x) paste0(x, "%")) +
    labs(
      title    = "NYHA Dose-Response Gradient: Phase V Prevalence",
      subtitle = "[Exploratory] Error bars = 95% CI; thresholds from MUSIC cohort",
      x        = "Cohort (ordered by disease severity)",
      y        = "Phase V Prevalence (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray50", size = 9)
    )

  if (save) {
    out <- file.path(fig_dir, "fig_nyha_gradient.png")
    ggsave(out, p, width = 8, height = 5, dpi = 150)
    message("Saved: ", out)
  }
  p
}


# =============================================================================
# Run all figures
# =============================================================================
if (!interactive()) {
  message("Generating all ECSoC figures...")
  plot_chi_distribution(save = TRUE)
  plot_nyha_gradient(save = TRUE)
  message("Done. Figures saved to figures/")
}
