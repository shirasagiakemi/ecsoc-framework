#!/usr/bin/env Rscript
# =============================================================================
# Generate Sample RR Interval Data for ECSoC Framework
# data/sample/generate_sample_data.R
#
# Creates two synthetic RR interval series with physiologically plausible
# properties for testing the ECSoC pipeline.
# NOTE: These are SIMULATED data, not real patient recordings.
# =============================================================================

set.seed(2025)

# --- Helper: Simulate fractal-like RR series via AR + noise -----------------
simulate_rr <- function(n, mean_rr, ar_coef, noise_sd, label) {
  x <- numeric(n)
  x[1] <- mean_rr
  for (i in 2:n) {
    x[i] <- ar_coef * x[i-1] + (1 - ar_coef) * mean_rr + rnorm(1, 0, noise_sd)
  }
  # Clip to physiologically plausible range (400–1600 ms)
  x <- pmax(400, pmin(1600, x))
  data.frame(rr_ms = round(x, 1), source = label)
}

# Healthy-like: higher variability, α₁ closer to 1.0
normal_rr <- simulate_rr(
  n = 500, mean_rr = 850, ar_coef = 0.85, noise_sd = 45,
  label = "simulated_healthy"
)

# CHF-like: lower variability, more correlated
chf_rr <- simulate_rr(
  n = 500, mean_rr = 780, ar_coef = 0.97, noise_sd = 18,
  label = "simulated_chf_like"
)

# Save
dir.create("data/sample", recursive = TRUE, showWarnings = FALSE)
write.csv(normal_rr, "data/sample/sample_rr_normal.csv", row.names = FALSE)
write.csv(chf_rr,    "data/sample/sample_rr_chf.csv",    row.names = FALSE)

cat("Sample data saved:\n")
cat("  data/sample/sample_rr_normal.csv  (N=500, simulated healthy-like)\n")
cat("  data/sample/sample_rr_chf.csv     (N=500, simulated CHF-like)\n")
cat("\nData format:\n")
cat("  rr_ms  — RR interval in milliseconds\n")
cat("  source — Data origin label\n")
