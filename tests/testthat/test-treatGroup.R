# ==============================================================================
# CREATE MOCK DATA FOR TESTING
# ==============================================================================
# Column 1: Patient ID
# Column 2: Drug Type (DrugA, DrugB, DrugC) -> 3 levels
# Column 3: Dosage (High, Low)      -> 2 levels
# Column 4: Measurement numeric data
mock_data <- data.frame(
  id = 1:6,
  drug = c("DrugA", "DrugA", "DrugB", "DrugB", "DrugC", "DrugC"),
  dosage = c("High", "Low", "High", "Low", "High","Low"),
  score = c(95, 82, 88, 74, 91, 85)
)

# ==============================================================================
# TEST CASE 1: Default Behavior (No Abbreviations)
# ==============================================================================
test_that("treatGroup combines columns with an underscore when abbreviations are NULL", {
  # Run the function
  result <- treatGroup(mock_data, number.of.treatments = 2)

  # Define exactly what we expect the output vector to be
  expected_treatment <- c("DrugA_High", "DrugB_High", "DrugC_High", "DrugA_Low","DrugB_Low", "DrugC_Low" )

  # Assertions
  expect_true("Treatment" %in% names(result))
  expect_equal(levels(result$Treatment),expected_treatment)
})

# ==============================================================================
# TEST CASE 2: Using Custom Abbreviations
# ==============================================================================

test_that("treatGroup applies custom abbreviations correctly when provided", {
  # Total 6 unique levels across columns 2 and 3 ("DrugA", "DrugB", "DrugC", "High", "Low")
  my_abbrevs <- c("A", "B", "C", "H", "L")

  # Run the function
  result <- treatGroup(mock_data, number.of.treatments = 2, my_abbreviations = my_abbrevs)

  # Define exactly what we expect the abbreviated codes to look like
  expected_treatment <- c("AH", "BH", "CH", "AL", "BL", "CL")

  # Assertions
  expect_equal(levels(result$Treatment), expected_treatment)
})

# ==============================================================================
# TEST CASE 3: Error Handling Safety Check
# ==============================================================================
test_that("treatGroup throws an explicit error if abbreviation length is incorrect", {
  # We need 4 abbreviations, but we only pass 3 here
  bad_abbrevs <- c("A", "B", "H")

  # expect_error checks if the code stops and looks for a specific error text fragment
  expect_error(
    treatGroup(mock_data, number.of.treatments = 2, my_abbreviations = bad_abbrevs),
    regexp = "does not match the total number of unique factor levels"
  )
})

# ==============================================================================
# TEST CASE 4: Handling NA cases
# ==============================================================================
test_that("treatGroup safely handles NA (missing data) values", {
  # Setup mock data where row 3 has an NA value in the dosage column
  na_mock_data <- data.frame(
    id = 1:3,
    drug = c("DrugA", "DrugA", "DrugB"),
    dosage = c("High", "Low", NA), # NA inserted here
    score = c(95, 82, 88)
  )

  # R's base interaction() function keeps NA values as NA by default
  result <- treatGroup(na_mock_data, number.of.treatments = 2)

  # Assertions
  # Row 3 should result in NA because one of its components is missing
  expect_true(is.na(result$Treatment[3]))
})

# ==============================================================================
# TEST CASE 5: Handling typo mistakes
# ==============================================================================
test_that("treatGroup accounts for typos as unique factor levels", {
  # Setup mock data where a human accidentally typed "Hgh" instead of "High"
  typo_mock_data <- data.frame(
    id = 1:4,
    drug = c("DrugA", "DrugA", "DrugB", "DrugB"),
    dosage = c("High", "Low", "Hgh", "Low"), # Typo level "Hgh" creates 3 levels total
    score = c(95, 82, 88, 91)
  )

  # Scenario A: Without abbreviations, it should just include the typo in the name
  result_no_abbrev <- treatGroup(typo_mock_data, number.of.treatments = 2)
  expect_equal(as.character(result_no_abbrev$Treatment[3]), c("DrugB_Hgh"))

  # Scenario B: If abbreviations are used, the typo increases the required vector size.
  # drug column has 1 level ("DrugA"), dosage column has 3 levels ("High", "Low", "Hgh").
  # Total required abbreviations = 4. Passing only 3 abbreviations should trigger our error.
  too_few_abbrevs <- c("A", "H", "L")

  expect_error(
    treatGroup(typo_mock_data, number.of.treatments = 2, my_abbreviations = too_few_abbrevs),
    regexp = "does not match the total number of unique factor levels"
  )
})

