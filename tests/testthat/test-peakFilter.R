library(testthat)
library(dplyr)
library(igraph)

# ==============================================================================
# 1. SETUP SHARED MOCK DATA FOR THE TESTS
# ==============================================================================
# This setup creates 5 mock peak features across 4 samples:
# - Peak1: Clean, perfect data.
# - Peak2: High missingness (3 out of 4 are NA) -> Should fail fNA.
# - Peak3: High variance / high noise -> Should fail fCV.
# - Peak4 & Peak5: Share almost identical RT (1.20 vs 1.21) and are perfectly
#   correlated (1.0). Peak5 has a larger m.z, so Peak4 should be filtered out.

mock_peaks <- data.frame(
  peak_id = c("Peak1", "Peak2", "Peak3", "Peak4", "Peak5"),
  SampleA = c(100,  NA,  10,  500,  500),
  SampleB = c(105,  NA, 900,  510,  510),
  SampleC = c(95,   NA,  15,  490,  490),
  SampleD = c(102, 500, 850,  520,  520),
  `m.z`   = c(150.1, 200.2, 250.3, 300.4, 350.5), # Peak5 has higher m.z than Peak4
  RT      = c(0.50,  0.80,  2.10,  1.20,  1.21),  # Peak4 and Peak5 have close RTs
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mock_treatment <- data.frame(
  sample = c("SampleA", "SampleB", "SampleC", "SampleD"),
  Treatment = c("Control", "Control", "Treated", "Treated"),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2. THE FORMAL TEST SUITE
# ==============================================================================

test_that("peakFilter correctly filters out high NA rows", {
  # Run only the NA filter (Threshold = 0.50)
  # Peak2 has 75% NAs across all samples, and does NOT have a complete group.
  result <- peakFilter(
    x = mock_peaks,
    y = mock_treatment,
    fNA = c('T', 0.5),
    fCV = FALSE,
    fRT = FALSE
  )

  # Assertions
  expect_false("Peak2" %in% result$peak_id)
  expect_true("Peak1" %in% result$peak_id)
  expect_equal(nrow(result), 4) # 4 out of 5 peaks should remain
})


test_that("peakFilter correctly filters out high CV rows", {
  # Run only the CV filter (Threshold = 0.50)
  # Peak3 has massive intensity swings (10, 900, 15, 850) -> CV will be > 1.0
  result <- peakFilter(
    x = mock_peaks,
    y = mock_treatment,
    fNA = FALSE,
    fCV = c('T', 0.5),
    fRT = FALSE
  )

  # Assertions
  expect_false("Peak3" %in% result$peak_id)
  expect_true("Peak1" %in% result$peak_id)
})


test_that("peakFilter keeps features with high NA/CV if they have a complete treatment group", {
  # Create a specialized peak that has NAs in the Control group,
  # but is 100% complete in the Treated group (SampleC and SampleD).
  complete_group_peak <- data.frame(
    peak_id = "Peak_Marker",
    SampleA = NA, SampleB = NA, SampleC = 600, SampleD = 610,
    `m.z` = 180.1, RT = 3.5, check.names = FALSE
  )
  test_data <- rbind(mock_peaks, complete_group_peak)

  # Run NA filter at a strict 20% allowance.
  # "Peak_Marker" has 50% total NAs, but should survive due to the complete Treated group rule.
  result <- peakFilter(test_data, mock_treatment, fNA = c('T', 0.2), fCV = FALSE, fRT = FALSE)

  # Assertions
  expect_true("Peak_Marker" %in% result$peak_id)
})


test_that("peakFilter network de-duplicates correlated RT peaks and keeps the max m.z", {
  # Run only the RT filter
  # Window = 0.02, Correlation Threshold = 0.95
  # Peak4 and Peak5 are within 0.01 min of each other and track perfectly.
  result <- peakFilter(
    x = mock_peaks,
    y = mock_treatment,
    fNA = FALSE,
    fCV = FALSE,
    fRT = c('T', 2, 0.02, 0.95)
  )

  # Assertions
  # The cluster containing Peak4 and Peak5 should collapse into just one row.
  # Because Peak5 has an m.z of 350.5 (higher than Peak4's 300.4), Peak5 must be kept.
  expect_true("Peak5" %in% result$peak_id)
  expect_false("Peak4" %in% result$peak_id)
})


test_that("peakFilter parallel execution matches sequential execution perfectly", {
  # Skip this test if the user doesn't have the parallel libraries loaded globally
  skip_if_not_installed("future.apply")

  # Run sequentially
  res_seq <- peakFilter(mock_peaks, mock_treatment, fNA = c('T', 0.5), fCV = c('T', 0.8), fRT = c('T', 2, 0.02, 0.95), parallel = FALSE)

  # Run in parallel
  res_par <- peakFilter(mock_peaks, mock_treatment, fNA = c('T', 0.5), fCV = c('T', 0.8), fRT = c('T', 2, 0.02, 0.95), parallel = TRUE, n_cores = 2)

  # Assertions
  expect_equal(res_seq, res_par)
})
