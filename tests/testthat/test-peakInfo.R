library(testthat)

# ==============================================================================
# 1. SETUP SHARED MOCK DATA
# ==============================================================================
# 5 peaks across 3 samples (columns 2, 3, 4). n = 3.
mock_data <- data.frame(
  peak_id = c("Peak_01", "Peak_02", "Peak_03", "Peak_04", "Peak_05"),
  Sample1 = c(10.0,   0, 100.0, 15.11111, NA),
  Sample2 = c(20.0,  50, 105.0, 15.22222, NA),
  Sample3 = c(30.0, 100,  95.0, 15.33333, NA),
  `m.z`   = c(100.1, 200.2, 300.3, 400.4, 500.5),
  RT      = c(1.1, 2.2, 3.3, 4.4, 5.5),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ==============================================================================
# 2. TEST SUITE
# ==============================================================================

test_that("peakInfo calculates basic statistics correctly with accurate types", {
  res <- peakInfo(mock_data, n = 3)

  # Ensure numbers remained numbers (Not converted to text strings by cbind)
  expect_type(res$m.z, "double")
  expect_type(res$Mean, "double")
  expect_type(res$Max, "double")

  # Verify specific math calculations on Peak_01 (Values: 10, 20, 30)
  expect_equal(res["Peak_01", "Max"], 30)
  expect_equal(res["Peak_01", "Mean"], 20.0000)
  expect_equal(res["Peak_01", "SD"], 10.0000)
  expect_equal(res["Peak_01", "CV"], 0.500)
})


test_that("peakInfo handles zero values as missing data (NA)", {
  res <- peakInfo(mock_data, n = 3)

  # Peak_02 has values: 0 (becomes NA), 50, 100.
  # Missing percentage should be 1 out of 3 = 0.33
  expect_equal(res["Peak_02", "pctNAs"], 0.33)

  # Mean should be calculated using only the remaining two numbers (50 and 100)
  expect_equal(res["Peak_02", "Mean"], 75.0000)
})


test_that("peakInfo rounds decimals strictly according to documentation limits", {
  res <- peakInfo(mock_data, n = 3)

  # Peak_04 values have 5 decimal places.
  # Mean should round to exactly 4 places: mean(15.11111, 15.22222, 15.33333) = 15.2222
  expect_equal(res["Peak_04", "Mean"], 15.2222)
})


test_that("peakInfo isolates a single feature row when ID parameter is supplied", {
  # Query just Peak_03
  single_res <- peakInfo(mock_data, n = 3, ID = "Peak_03")

  expect_equal(nrow(single_res), 1)
  expect_equal(rownames(single_res), "Peak_03")
  expect_equal(single_res$Max, 105)
})


test_that("peakInfo throws an expected error if requested ID is missing", {
  expect_error(
    peakInfo(mock_data, n = 3, ID = "Fake_Peak_Name"),
    regexp = "could not be found in the dataset"
  )
})
