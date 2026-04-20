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
  expect_error(barPeaks("not a dataframe", '355.12', 5))
  expect_error(barPeaks(toy_df, '355.12', 10))
  expect_error(barPeaks(toy_df, 100, 5))
})


# ------------------------
# Test 3: edge case with empty dataframe
# ------------------------
test_that("barPeaks handles empty data", {
  empty_df <- data.frame()
  expect_error(barPeaks(empty_df, '355.12', 5))
})

# ------------------------
# Test 4: function respects number of samples
# ------------------------
test_that("barPeaks respects number of samples", {
  res <- barPeaks(toy_df, mass = "355.12", n = 3)
  expect_equal(ncol(res$data), 3)
})

# ------------------------
# Test 5: selected peaks have correct mass
# ------------------------

test_that("Selected peaks have correct mass", {
  m <- "355.12"
  res <- barPeaks(toy_df, mass = m, n = 5)
  selected <- toy_df[toy_df$peak %in% res$peak_names, ]

  expect_true(all(selected$m.z == m))
})

# ------------------------
# Test 6: function selects lowest CV peaks
# ------------------------

test_that("barPeaks selects lowest CV peaks", {

  res <- barPeaks(toy_df, mass = "355.12", n = 5)

  # CV should be sorted (lowest first)
  expect_true(all(diff(res$CV) >= 0))
})

# ------------------------
# Test 7: function errors when mass is not present
# ------------------------

test_that("barPeaks errors when mass is not present", {
  expect_error(
    barPeaks(toy_df, mass = "999", n = 3),
    "This mass is not in the data frame"
  )
})

# ------------------------
# Test 8: function handles less than 3 peaks
# ------------------------

test_that("barPeaks handles less than 3 peaks", {

  res <- barPeaks(toy_df, mass = "170.42", n = 5)

  expect_true(length(res$mids) <= 2)
})
