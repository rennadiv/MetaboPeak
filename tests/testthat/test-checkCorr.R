library(testthat)
library(MetaboPeak)

# ------------------------
# Toy dataset
# ------------------------
toy_df <- data.frame(
  ID = c("p1", "p2", "p3"),
  sample1 = c(1, 2, 3),
  sample2 = c(2, 4, 6),
  sample3 = c(3, 6, 9),
  m.z = c(100.1, 200.2, 300.3),
  RT = c(5.0, 5.0, 5.0)
)

# ------------------------
# 1. Basic functionality
# ------------------------
test_that("checkCorr returns numeric correlation", {

  res <- checkCorr(
    toy_df,
    masses = c(100.1, 200.2),
    RT = 5.0,
    n = 3
  )

  expect_true(is.numeric(res))
  expect_length(res, 1)
})

# ------------------------
# 2. Perfect correlation
# ------------------------
test_that("checkCorr detects perfect correlation", {

  res <- checkCorr(
    toy_df,
    masses = c(100.1, 200.2),
    RT = 5.0,
    n = 3
  )

  expect_equal(res, 1)  # perfectly correlated
})

# ------------------------
# 3. Works with two data frames
# ------------------------
test_that("checkCorr works with two data frames", {

  toy_df2 <- toy_df
  toy_df2$sample1 <- c(3, 2, 1)  # reversed

  res <- checkCorr(
    toy_df,
    toy_df2,
    masses = c(100.1, 100.1),
    RT = 5.0,
    n = 3
  )

  expect_true(is.numeric(res))
})

# ------------------------
# 4. Handles NA values
# ------------------------
test_that("checkCorr handles missing values", {

  df_na <- toy_df
  df_na$sample2[1] <- NA

  res <- checkCorr(
    df_na,
    masses = c(100.1, 200.2),
    RT = 5.0,
    n = 3
  )

  expect_true(is.numeric(res))
})

# ------------------------
# 5. Error: first mass not found
# ------------------------
test_that("checkCorr errors when first mass is missing", {

  expect_error(
    checkCorr(
      toy_df,
      masses = c(999, 200.2),
      RT = 5.0,
      n = 3
    ),
    "First mass is not in the data table"
  )
})

# ------------------------
# 6. Error: second mass not found
# ------------------------
test_that("checkCorr errors when second mass is missing", {

  expect_error(
    checkCorr(
      toy_df,
      masses = c(100.1, 999),
      RT = 5.0,
      n = 3
    ),
    "Second mass is not in the data table"
  )
})
