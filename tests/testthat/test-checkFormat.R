library(testthat)
library(MetaboPeak)

# -------------------------
# Toy correct dataset
# -------------------------
toy_df <- data.frame(
  PeakID = c("p1", "p2"),
  sample1 = c(1, 2),
  sample2 = c(3, 4),
  sample3 = c(5, 6),
  m.z = c(100.1, 200.2),
  RT = c(5.0, 6.0)
)

# -------------------------
# 1. Correct input works
# -------------------------
test_that("checkFormat accepts valid dataset", {
  expect_true(checkFormat(toy_df, n = 3))
})

# -------------------------
# 2. Not a data frame
# -------------------------
test_that("checkFormat fails if input is not a data.frame", {
  expect_error(
    checkFormat(matrix(1:10, ncol = 5), n = 3),
    "data.frame"
  )
})

# -------------------------
# 3. Missing m.z column
# -------------------------
test_that("checkFormat fails if m.z is missing", {
  df <- toy_df
  df$m.z <- NULL

  expect_error(
    checkFormat(df, n = 3),
    "m.z"
  )
})

# -------------------------
# 4. Missing RT column
# -------------------------
test_that("checkFormat fails if RT is missing", {
  df <- toy_df
  df$RT <- NULL

  expect_error(
    checkFormat(df, n = 3),
    "RT"
  )
})

# -------------------------
# 5. Wrong number of samples
# -------------------------
test_that("checkFormat fails if n does not match columns", {
  expect_error(
    checkFormat(toy_df, n = 10),
    "sample"
  )
})

# -------------------------
# 6. Non-numeric sample columns
# -------------------------
test_that("checkFormat fails if sample columns are not numeric", {
  df <- toy_df
  df$sample2 <- c("a", "b")

  expect_error(
    checkFormat(df, n = 3),
    "numeric"
  )
})
