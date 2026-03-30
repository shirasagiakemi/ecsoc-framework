# =============================================================================
# ECSoC Framework — Unit Tests
# tests/test_dfa.R
#
# Run with: source("tests/test_dfa.R")
# =============================================================================

source("R/01_dfa.R")
source("R/02_chi_calculation.R")

cat("=== ECSoC Unit Tests ===\n\n")

passed <- 0
failed <- 0

assert <- function(condition, name) {
  if (isTRUE(condition)) {
    cat(sprintf("  ✓ PASS: %s\n", name))
    passed <<- passed + 1
  } else {
    cat(sprintf("  ✗ FAIL: %s\n", name))
    failed <<- failed + 1
  }
}

# --- Test 1: DFA runs on valid input -----------------------------------------
cat("--- DFA Tests ---\n")
set.seed(42)
rr_test <- rnorm(1000, mean = 850, sd = 50)
rr_test <- pmax(400, pmin(1600, rr_test))

result <- tryCatch(compute_dfa(rr_test), error = function(e) NULL)
assert(!is.null(result), "compute_dfa() runs without error on valid input")

if (!is.null(result)) {
  assert(is.numeric(result$alpha1), "alpha1 is numeric")
  assert(is.numeric(result$alpha2), "alpha2 is numeric")
  assert(is.numeric(result$r_squared), "r_squared is numeric")
  assert(result$r_squared >= 0 && result$r_squared <= 1, "R² is in [0,1]")
  assert(result$alpha1 > 0 && result$alpha1 < 3, "alpha1 in plausible range (0–3)")
  assert(result$alpha2 > 0 && result$alpha2 < 3, "alpha2 in plausible range (0–3)")
}

# --- Test 2: DFA fails gracefully on bad input --------------------------------
too_short <- rnorm(50, 850, 50)
err <- tryCatch({ compute_dfa(too_short); NULL }, error = function(e) e)
assert(inherits(err, "error"), "compute_dfa() errors on series too short")

# --- Test 3: CHI calculation --------------------------------------------------
cat("\n--- CHI Tests ---\n")
chi_pos <- compute_chi(1.2, 0.9)
assert(abs(chi_pos - 0.6) < 1e-10, "CHI = 2*(1.2 - 0.9) = 0.6")

chi_neg <- compute_chi(0.7, 1.0)
assert(abs(chi_neg - (-0.6)) < 1e-10, "CHI = 2*(0.7 - 1.0) = -0.6")

chi_zero <- compute_chi(1.0, 1.0)
assert(abs(chi_zero) < 1e-10, "CHI = 0 when alpha1 == alpha2")

# --- Test 4: CHI domain classification ----------------------------------------
assert(classify_chi_domain(0.5) == "supercritical", "CHI=0.5 → supercritical")
assert(classify_chi_domain(-0.5) == "subcritical",  "CHI=-0.5 → subcritical")
assert(classify_chi_domain(0.0) == "near-critical", "CHI=0.0 → near-critical")

# --- Test 5: Phase classification ---------------------------------------------
cat("\n--- Phase Classification Tests ---\n")
source("R/02_chi_calculation.R")

p1 <- classify_phase(0.7, 0.98, -0.3)
assert(p1$phase == "I",   "alpha1=0.7, R²=0.98 → Phase I")

p2 <- classify_phase(1.0, 0.98, 0.0)
assert(p2$phase == "II",  "alpha1=1.0, R²=0.98 → Phase II")

p3 <- classify_phase(1.3, 0.98, 0.4)
assert(p3$phase == "III", "alpha1=1.3, R²=0.98 → Phase III")

pVa <- classify_phase(0.7, 0.90, -0.4)
assert(pVa$phase == "V" && pVa$subtype == "Va",
       "alpha1=0.7, R²=0.90, CHI=-0.4 → Phase Va")

pVb <- classify_phase(1.3, 0.90, 0.4)
assert(pVb$phase == "V" && pVb$subtype == "Vb",
       "alpha1=1.3, R²=0.90, CHI=+0.4 → Phase Vb")

# --- Test 6: Healthy reference frame ------------------------------------------
cat("\n--- CHI Reference Frame Tests ---\n")
interp_healthy <- interpret_chi(0.44)
assert(grepl("healthy", interp_healthy$interpretation, ignore.case = TRUE),
       "CHI=+0.44 interpreted as within healthy range")

interp_chf <- interpret_chi(-0.5)
assert(grepl("CHF|subcritical|disease", interp_chf$interpretation, ignore.case = TRUE),
       "CHI=-0.5 interpreted as CHF/subcritical")

# --- Summary ------------------------------------------------------------------
cat(sprintf("\n=== Results: %d passed, %d failed ===\n", passed, failed))
if (failed == 0) {
  cat("All tests passed.\n")
} else {
  cat("Some tests FAILED. Please review the implementation.\n")
}
