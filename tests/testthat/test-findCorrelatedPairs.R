library(testthat)

# ==============================================================================
# 1. SETUP SHARED MOCK DATA FOR THE TESTS
# ==============================================================================
# This setup creates 5 mock peak features across 4 samples:
# - peak_id: Name of the peaks
# - S1-S3: Samples with their abundances
# - m.z: mass to charge ration
# - RT: retention time

mock_data <- data.frame(
  peak_id = c("PeakA", "PeakB", "PeakC"),
  S1 = c(10, 11, 500),
  S2 = c(20, 21,  10),
  S3 = c(30, 31, 900),
  `m.z` = c(150.1, 152.1, 400.5),
  RT = c(5.00, 5.02, 5.04), # All close in time
  stringsAsFactors = FALSE, check.names = FALSE
)

# ==============================================================================
# 2. THE FORMAL TEST SUITE
# ==============================================================================

test_that("findCorrelatedPairs discovers and ranks correct pairs internally", {
  # PeakA and PeakB increase perfectly together (10->30 vs 11->31). Correlation = 1.0
  # PeakC bounces around randomly.
  res <- findCorrelatedPairs(mock_data, n = 3, rt_window = 0.05, cor_threshold = 0.90)

  expect_equal(nrow(res), 1)
  expect_equal(res$Peak_X[1], "PeakA")
  expect_equal(res$Peak_Y[1], "PeakB")
  expect_equal(res$Correlation[1], 1.0)
})

test_that("findCorrelatedPairs returns empty dataframe cleanly if no pairs pass criteria", {
  # Set a correlation threshold of 1.1 (impossible)
  res <- suppressMessages(findCorrelatedPairs(mock_data, n = 3, rt_window = 0.05, cor_threshold = 1.1))
  expect_equal(nrow(res), 0)
  expect_true(is.data.frame(res))
})

# ==============================================================================
# 3. THE VISUALIZATION TESTS
# ==============================================================================

test_that("plotPeakCorrelation runs smoothly and throws appropriate missing ID errors", {
  # Should run silently without errors
  expect_silent(
    plotPeakCorrelation(mock_data, peak_id_x = "PeakA", peak_id_y = "PeakB", n = 3)
  )

  # Should crash cleanly if a user passes a wrong name tag
  expect_error(
    plotPeakCorrelation(mock_data, peak_id_x = "PeakA", peak_id_y = "Fake_Name", n = 3)
  )
})

