# ==============================================================================
# 1. SETUP SHARED MOCK DATA FOR THE TESTS
# ==============================================================================
# This setup creates 5 mock peak features across 3 samples:
# - peak_id: Name of the peaks
# - S1-S3: Samples with their abundances
# - m.z: mass to charge ratio
# - RT: retention time (varying from noisy edges to biological region)

mock_peaks <- data.frame(
  peak_id = c("Peak1_Noise", "Peak2", "Peak3", "Peak4", "Peak5_Noise"),
  S1      = c(10,           NA,    500,     15,    100),
  S2      = c(950,          50,    10,      15,    200),
  S3      = c(30,           100,   20,      15,    300),
  `m.z`   = c(593.128,      200.2, 593.130, 400.4, 500.5),
  RT      = c(0.4,          2.2,   4.5,     12.1,  28.7),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ==============================================================================
# 2. THE FORMAL TEST SUITE
# ==============================================================================

test_that("interestRegion filters retention time correctly within inclusive boundaries", {
  # Isolating biological window between 1.5 and 15.0 minutes.
  # This should keep Peak2, Peak3, and Peak4 while removing edge noise.
  res <- interestRegion(mock_peaks, start = 1.5, end = 15.0)

  # Assertions
  expect_equal(nrow(res), 3)
  expect_true(all(res$RT >= 1.5 & res$RT <= 15.0))
  expect_false("Peak1_Noise" %in% res$peak_id)
  expect_false("Peak5_Noise" %in% res$peak_id)
})


test_that("interestRegion throws a clear error if the mandatory RT column is missing", {
  # Create bad data framework missing the 'RT' column entirely
  bad_peaks <- mock_peaks
  bad_peaks$RT <- NULL

  # Assertions
  expect_error(
    interestRegion(bad_peaks, start = 1.5, end = 15.0),
    "The input data frame does not have a column named 'RT'"
  )
})


test_that("interestRegion stops execution if boundary inputs are non-numeric", {
  # Assertions passing character values to numeric thresholds
  expect_error(
    interestRegion(mock_peaks, start = "one", end = 15.0),
    "Both 'start' and 'end' must be numeric values"
  )
  expect_error(
    interestRegion(mock_peaks, start = 1.5, end = "fifteen"),
    "Both 'start' and 'end' must be numeric values"
  )
})
