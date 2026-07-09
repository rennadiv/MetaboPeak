library(testthat)
library(dplyr)

# ==============================================================================
# 1. SETUP SHARED MOCK DATA FOR THE TESTS
# ==============================================================================
# This setup creates 5 mock peak features across 4 samples:
# - peak_id: Name of the peaks
# - S1-S3: Samples with their abundances
# - m.z: mass to charge ratio
# - RT: retention time
mock_peaks <- data.frame(
  peak_id = c("Peak1", "Peak2", "Peak3", "Peak4", "Peak5"),
  S1      = c(100,     200,     500,    15,    100),
  S2      = c(950,     250,     600,    15,    200),
  S3      = c(300,     220,     700,    15,    300),
  `m.z`   = c(150.005, 250.500, 320.115, 400.4, 500.5),
  RT      = c(5.20,    2.20,    8.50,    4.40,  5.50),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Reference library definitions
mock_library <- data.frame(
  name = c("Glucose-Derivative", "Tryptophan-Isomer"),
  `m.z` = c(150.000, 320.120),
  RT   = c(5.25, 8.40),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ==============================================================================
# 2. THE FORMAL TEST SUITE
# ==============================================================================

test_that("libraryHits correctly annotates true matching features", {
  # Peak1 and Peak3 match the target parameters inside the mock library data frame
  res <- libraryHits(mock_peaks, mock_library, RT.range = 0.25, mz.range = 0.01, unique = TRUE)

  # Assertions
  expect_true("Glucose-Derivative" %in% res$unique_names)
  expect_true("Tryptophan-Isomer" %in% res$unique_names)
  expect_equal(length(res$unique_names), 2)
})


test_that("libraryHits marks unassigned features correctly as unknown", {
  res <- libraryHits(mock_peaks, mock_library, RT.range = 0.25, mz.range = 0.01, unique = TRUE)

  # Verify formatting layout structures of lib_all
  # Peak2 (250.5) has no match, it must read "unknown"
  peak2_row <- res$lib_all[res$lib_all$Alignment.ID == "Peak2", ]
  expect_equal(peak2_row$Compound, "unknown")

  # Peak1 must show its resolved compound assignment text string
  peak1_row <- res$lib_all[res$lib_all$Alignment.ID == "Peak1", ]
  expect_equal(peak1_row$Compound, "Glucose-Derivative")
})


test_that("libraryHits returns an empty structure gracefully if nothing matches", {
  # Force a tiny mass search criteria that catches zero elements
  res <- suppressMessages(libraryHits(mock_peaks, mock_library, RT.range = 0.001, mz.range = 0.0001))

  # Assertions
  expect_equal(nrow(res$lib_small), 0)
  expect_equal(length(res$unique_names), 0)
  expect_true(all(res$lib_all$Compound == "unknown"))
})
