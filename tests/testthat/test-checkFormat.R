# ==============================================================================
# 1. SETUP SHARED MOCK DATA FOR THE TESTS
# ==============================================================================
# This setup creates 5 mock peak features across 4 columns (3 samples + 1 ID):
# - peak_id: Name of the peaks (Column 1)
# - S1-S3: Samples with their abundances (Columns 2 to 4)
# - m.z: mass to charge ratio (Column 5)
# - RT: retention time (Column 6)

mock_peaks <- data.frame(
  peak_id = c("Peak1", "Peak2", "Peak3_Isobaric", "Peak4", "Peak5"),
  S1      = c(10,      NA,     500,              15,    100),
  S2      = c(950,     50,      10,              15,    200),
  S3      = c(30,     100,      20,              15,    300),
  `m.z`   = c(593.128, 200.200, 593.130,         400.4, 500.5),
  RT      = c(4.2,     2.2,     8.5,             4.4,   5.5),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ==============================================================================
# 2. THE FORMAL TEST SUITE
# ==============================================================================

test_that("checkFormat returns TRUE for a perfectly structured dataset", {
  # mock_peaks contains 3 sample columns (S1, S2, S3), m.z, and RT.
  res <- checkFormat(mock_peaks, n = 3)

  # Assertions
  expect_true(res)
})


test_that("checkFormat catches missing metadata columns and structural failures", {
  # 1. Test failure if not a data frame
  expect_error(
    checkFormat(as.matrix(mock_peaks), n = 3),
    "Input must be a data.frame"
  )

  # 2. Test missing 'm.z' column
  no_mz <- mock_peaks
  no_mz$`m.z` <- NULL
  expect_error(
    checkFormat(no_mz, n = 3),
    "Column 'm.z' is missing"
  )

  # 3. Test missing 'RT' column
  no_rt <- mock_peaks
  no_rt$RT <- NULL
  expect_error(
    checkFormat(no_rt, n = 3),
    "Column 'RT' is missing"
  )
})


test_that("checkFormat validates sample column count and numeric data types", {
  # 1. Test when 'n' claims more sample columns than are actually present.
  # The function should complain about finding 4 non-metadata columns instead of 5.
  expect_error(
    checkFormat(mock_peaks, n = 5),
    "Expected 5 sample columns, but found 4"
  )

  # 2. Test when sample abundance columns contain non-numeric data (e.g. text characters)
  bad_types <- mock_peaks
  bad_types$S1 <- c("high", "low", "medium", "low", "high")
  expect_error(
    checkFormat(bad_types, n = 3),
    "Sample abundance column 'S1' must be numeric"
  )
})

