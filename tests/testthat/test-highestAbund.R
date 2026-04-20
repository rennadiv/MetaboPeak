library(testthat)
library(MetaboPeak)

# ------------------------
# Toy dataset
# ------------------------
toy_df <- data.frame(
  ID = c("p1", "p2"),
  sample1 = c(10, 5),
  sample2 = c(20, 15),
  sample3 = c(15, 25),
  m.z = c(100.1, 200.2),
  RT = c(5.0, 6.0)
)

# ------------------------
# 1. Basic functionality
# ------------------------
test_that("highestAbund returns correct sample name", {

  res <- highestAbund(toy_df, mass = "100.1", n = 3)

  expect_equal(res, "sample2")  # 20 is highest
})

# ------------------------
# 2. Returns character
# ------------------------
test_that("highestAbund returns character", {

  res <- highestAbund(toy_df, mass = "100.1", n = 3)

  expect_true(is.character(res))
  expect_length(res, 1)
})

# ------------------------
# 3. Handles NA values
# ------------------------
test_that("highestAbund handles NA values", {

  df_na <- toy_df
  df_na$sample2[1] <- NA

  res <- highestAbund(df_na, mass = "100.1", n = 3)

  expect_true(is.character(res))
})

# ------------------------
# 4. Error when multiple peaks match mass
# ------------------------
test_that("highestAbund errors when mass is not unique", {

  df_dup <- toy_df
  df_dup$m.z[2] <- 100.1  # duplicate mass

  expect_error(
    highestAbund(df_dup, mass = "100.1", n = 3),
    "not unique"
  )
})

# ------------------------
# 5. Error when mass not found (IMPORTANT)
# ------------------------
test_that("highestAbund errors when mass not found", {

  expect_error(
    highestAbund(toy_df, mass = "999", n = 3)
  )
})
