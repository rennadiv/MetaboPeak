# tests/testthat/test-barPeaks.R

library(testthat)
library(MetaboPeak)  # make sure your package is loaded

# ------------------------
# Toy dataset for testing
# ------------------------
toy_df <- data.frame(
  peak = c("p1", "p2","p3","p4"),
  sample1 = c(30, 20, 12, 18),
  sample2 = c(15, 22, 20, 17),
  sample3 = c(4, 21, 20, 17),
  sample4 = c(12, 20, 20, 17),
  sample5 = c(21, 22, 20, 16),
  m.z = c(355.12, 355.12, 170.42, 89.6),
  RT = c(8.12, 8.5, 3.4, 2.1)
)

# ------------------------
# Test 1: function returns structured output
# ------------------------
test_that("barPeaks returns structured output", {

  res <- barPeaks(toy_df, mass = "355.12", n = 5)

  expect_true(is.list(res))
  expect_true("mids" %in% names(res))
  expect_true(all(sapply(res$mids, is.numeric)))
})

# ------------------------
# Test 2: error on wrong input
# ------------------------
test_that("barPeaks errors with wrong input", {
  expect_error(barPeaks("not a dataframe", '355.12', 4))
  expect_error(barPeaks(toy_df, '355.12', 10))

})


# ------------------------
# Test 3: edge case with empty dataframe
# ------------------------
test_that("barPeaks handles empty data", {
  empty_df <- data.frame()
  expect_error(barPeaks(empty_df, '355.12', 4))
})

