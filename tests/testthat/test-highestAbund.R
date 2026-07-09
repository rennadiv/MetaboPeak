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
  S1      = c(10,      NA,     500,              15,    100),
  S2      = c(950,     50,      10,              15,    200),
  S3      = c(30,     100,      20,              15,    300),
  `m.z`   = c(593.128, 200.200, 593.130,         400.4, 500.5), # Peak1 and Peak3 are isobaric
  RT      = c(4.2,     2.2,     8.5,             4.4,   5.5),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ==============================================================================
# 2. THE FORMAL TEST SUITE
# ==============================================================================

test_that("highestAbund returns the correct sample column name for a unique peak", {
  # Peak1 has its max value (950) in column S2.
  # Using a small tolerance (0.001) isolates Peak1 perfectly.
  res <- highestAbund(mock_peaks, mass = 593.128, n = 3, tolerance = 0.001)

  # Assertions
  expect_equal(res, "S2")
})


test_that("highestAbund safely resolves multiple matching peaks without crashing", {
  # Searching for 593.13 with a 0.02 tolerance catches BOTH Peak1 and Peak3_Isobaric.
  # Peak1 max value is 950 (in S2). Peak3_Isobaric max value is 500 (in S1).
  # The overall absolute max is 950, so it must return "S2" cleanly.
  res <- highestAbund(mock_peaks, mass = 593.13, n = 3, tolerance = 0.02)

  # Assertions
  expect_equal(res, "S2")
})


test_that("highestAbund safely ignores legacy 0 values", {
  # Modify Peak5 to have a 0 value to make sure it doesn't cause index math issues
  zero_mock_peaks <- mock_peaks
  zero_mock_peaks[5, "S1"] <- 0

  # Peak5 values are now: 0 (ignored/NA), 200, 300. Max is 300 in S3.
  res <- highestAbund(zero_mock_peaks, mass = 500.5, n = 3)

  # Assertions
  expect_equal(res, "S3")
})


test_that("highestAbund errors out cleanly if the requested mass is completely missing", {
  # Attempt to look up a mass window that does not exist in our data schema
  expect_error(
    highestAbund(mock_peaks, mass = 999.99, n = 3),
    regexp = "was not found within a tolerance"
  )
})

