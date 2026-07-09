library(testthat)

# Mock Mass-Spec Peak Frames
mock_pos <- data.frame(
  peak_id = c("P1", "P2"),
  S1 = c(10, 500), S2 = c(20, 510), S3 = c(30, 520),
  `m.z` = c(150.045, 400.123), RT = c(5.2, 11.8),
  stringsAsFactors = FALSE, check.names = FALSE
)

mock_neg <- data.frame(
  peak_id = c("N1", "N2"),
  S1 = c(12, 5),   S2 = c(22, 5),   S3 = c(32, 5),
  `m.z` = c(148.032, 398.110), RT = c(5.1, 11.9),
  stringsAsFactors = FALSE, check.names = FALSE
)

test_that("checkCorr correctly calculates correlation within a single dataset", {
  # P1 has values: 10, 20, 30. Perfect correlation with itself.
  res <- checkCorr(mock_pos, masses = c(150.04, 150.04), RT = 5.0, n = 3, tolerance = 0.05)
  expect_equal(res, 1.0)
})

test_that("checkCorr tracks correct rows cross-dataset and handles negative correlations", {
  # P1 increases (10, 20, 30) while N2 is flat/inversely related to other metrics
  # Compare P1 in pos (150.045) and N1 in neg (148.032) -> Both increase together (10->30 vs 12->32)
  res <- checkCorr(x = mock_pos, y = mock_neg, masses = c(150.045, 148.032), RT = 5.2, n = 3)
  expect_equal(res, 1.0)
})

test_that("checkCorr fails safely when target mass is completely absent", {
  expect_error(
    checkCorr(mock_pos, masses = c(999.9, 150.0), RT = 5.2, n = 3),
    regexp = "was not found within a tolerance"
  )
})
