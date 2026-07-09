library(testthat)

# ==============================================================================
# 1. SETUP SHARED MOCK DATA FOR THE TESTS
# ==============================================================================
# This setup creates 5 mock peak features across 4 samples:
# - peak_id: Name of the peaks
# - S1-S3: Samples with their abundances
# - m.z: mass to charge ratio
# - RT: retention time

mock_peaks <- data.frame(
  peak_id = c("Peak1", "Peak2", "Peak3_Isobaric", "Peak4", "Peak5"),
  S1      = c(100,     NA,     500,              15,    100),
  S2      = c(950,     50,      NA,              15,    200),
  S3      = c(300,    100,      20,              15,      0), # Contains a zero
  `m.z`   = c(447.128, 200.200, 447.130,         400.4, 500.5),
  RT      = c(5.25,    2.20,    12.38,           4.40,  5.50),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mock_treatment <- data.frame(
  sample_name = c("S1", "S2", "S3"),
  Treatment   = c("Control", "Control", "Treated"),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2. THE FORMAL TEST SUITE
# ==============================================================================

test_that("treatCases identifies correct treatment groups with complete cases", {
  # Peak1 has valid values in S1 (100), S2 (950), and S3 (300).
  # S1 and S2 belong to "Control". S3 belongs to "Treated".
  # Both groups should be complete.
  res <- treatCases(mock_peaks, mock_treatment, mass = 447.128, RT = 5.25, n = 3)

  # Assertions
  expect_true("Control" %in% res)
  expect_true("Treated" %in% res)
  expect_equal(length(res), 2)
})


test_that("treatCases filters out treatment groups containing NA values", {
  # Peak3_Isobaric has an NA value in S2.
  # S1 and S2 are "Control", so "Control" is incomplete.
  # S3 is "Treated" and contains a valid value (20), so "Treated" should be complete.
  res <- treatCases(mock_peaks, mock_treatment, mass = 447.130, RT = 12.38, n = 3)

  # Assertions
  expect_true("Treated" %in% res)
  expect_false("Control" %in% res)
  expect_equal(length(res), 1)
})


test_that("treatCases treats legacy 0 values as incomplete cases", {
  # Peak5 has a 0 value in S3.
  # S3 belongs to "Treated", so "Treated" should be flagged as incomplete and omitted.
  # S1 (100) and S2 (200) belong to "Control" and are valid, so "Control" passes.
  res <- treatCases(mock_peaks, mock_treatment, mass = 500.5, RT = 5.5, n = 3)

  # Assertions
  expect_true("Control" %in% res)
  expect_false("Treated" %in% res)
  expect_equal(length(res), 1)
})


test_that("treatCases errors out cleanly if the requested peak is missing", {
  # Attempt to run a coordinate tracking match on values completely outside the table boundaries
  expect_error(
    treatCases(mock_peaks, mock_treatment, mass = 999.9, RT = 99.9, n = 3),
    regexp = "No peak found matching m/z"
  )
})
