# ==============================================================================
# 1. SETUP SHARED MOCK DATA FOR THE TESTS
# ==============================================================================
# This setup tests varied background noise cases:
# - peak_id: Structured as purely numeric characters to test the Column 1 protection lock.
# - S1: Contains explicit NA values.
# - S2: Contains low-intensity baseline noise below threshold (500).
# - S3: Contains a high, valid target peak signal (50000).

mock_impute_peaks <- data.frame(
  peak_id = c(101,   102,   103),  # Numeric style names to test the safety lock
  S1      = c(NA,    20000, 30000),
  S2      = c(500,   15000, 25000),
  S3      = c(50000, 10000, 40000),
  `m.z`   = c(200.1, 300.2, 400.3),
  RT      = c(1.5,   2.5,   3.5),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2. THE FORMAL TEST SUITE
# ==============================================================================

test_that("peakImpute substitutes low values and NAs with a fraction of row minimums", {
  # For Row 1: Valid elements are strictly S3 (50000).
  # Minimum valid value = 50000. Imputation replacement value = 50000 * 0.1 = 5000.
  res <- peakImpute(mock_impute_peaks, threshold = 1000, fraction = 0.1, replace_low = TRUE, replace_na = TRUE)

  # Assertions
  expect_equal(res$S1[1], 5000) # NA replaced by row LOD fraction
  expect_equal(res$S2[1], 5000) # Low value (500 <= 1000) replaced by row LOD fraction
  expect_equal(res$S3[1], 50000) # Valid high element remains completely untouched
})

test_that("peakImpute safety lock protects structural IDs and metadata columns", {
  res <- peakImpute(mock_impute_peaks, threshold = 1000, fraction = 0.1)

  # Assertions
  expect_equal(res$peak_id[1], 101) # Column 1 must remain entirely untouched!
  expect_equal(res$RT[1], 1.5)      # Metadata tracks must remain untouched
})
