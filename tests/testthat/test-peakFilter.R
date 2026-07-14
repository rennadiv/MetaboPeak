# ==============================================================================
# 1. SETUP SHARED MOCK DATA FOR THE TESTS
# ==============================================================================
# We simulate 5 mock peak features across 4 samples split into 2 tree species cohorts:
# - Beech (B1, B2) and Spruce (S1, S2)
# - Peak1: Complete in Spruce, but completely missing (0) in Beech (Biomarker).
# - Peak2: High missingness globally (75% zeros), no complete group.
# - Peak3: High variance (CV) in both groups, complete.
# - Peak4 & Peak5: Overlapping RTs (4.201 & 4.202) and strongly correlated.

mock_peaks <- data.frame(
  peak_id = c("Peak1_Biomarker", "Peak2_BadNoise", "Peak3_HighCV", "Peak4_Network1", "Peak5_Network2"),
  B1      = c(0,               0,               10,             500,              510),
  B2      = c(0,               0,               950,            510,              520),
  S1      = c(400,             0,               20,             520,              530),
  S2      = c(420,             15,              890,            530,              540),
  `m.z`   = c(150.01,          220.05,          310.12,         447.22,           120.02), # Peak4 has higher m.z than Peak5
  RT      = c(2.50,            14.20,           18.50,          4.201,            4.202),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mock_treatment <- data.frame(
  Sample = c("B1", "B2", "S1", "S2"),
  Treatment = c("Beech", "Beech", "Spruce", "Spruce"),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2. THE FORMAL TEST SUITE
# ==============================================================================

test_that("peakFilter protects cohort biomarkers when protect_complete_groups = TRUE", {
  # Peak1_Biomarker is 100% missing in Beech but 100% complete in Spruce.
  # Global missingness is 50%. Setting threshold to strict 30% (0.3) should keep it
  # because protect_complete_groups = TRUE acts as a biological shield.

  res <- peakFilter(
    x = mock_peaks,
    y = mock_treatment,
    fNA = c('T', 0.3),
    fCV = FALSE,
    fRT = FALSE,
    protect_complete_groups = TRUE
  )

  # Assertions
  expect_true("Peak1_Biomarker" %in% res$peak_id)
  expect_false("Peak2_BadNoise" %in% res$peak_id) # Dropped because it has no complete cohort
})


test_that("peakFilter strictly applies thresholds when protect_complete_groups = FALSE", {
  # Disabling group protection forces raw mathematical sweeping.
  # Peak1_Biomarker (50% missingness) will now fail the strict 30% (0.3) filter.

  res <- peakFilter(
    x = mock_peaks,
    y = mock_treatment,
    fNA = c('T', 0.3),
    fCV = FALSE,
    fRT = FALSE,
    protect_complete_groups = FALSE
  )

  # Assertions
  expect_false("Peak1_Biomarker" %in% res$peak_id)
})


test_that("peakFilter isolates high mass vertices in overlapping RT network windows", {
  # Peak4 and Peak5 share a tight elution pattern and correlate perfectly.
  # Peak4 features a higher mass vector (447.22 vs 120.02) and should be preserved.

  res <- peakFilter(
    x = mock_peaks,
    y = mock_treatment,
    fNA = FALSE,
    fCV = FALSE,
    fRT = c('T', 2, 0.05, 0.90),
    protect_complete_groups = FALSE
  )

  # Assertions
  expect_true("Peak4_Network1" %in% res$peak_id)
  expect_false("Peak5_Network2" %in% res$peak_id)
})


test_that("peakFilter handles parameter omissions dynamically in console logs", {
  # Wrapping the function call to read text logs outputted to console.
  # Disabling fCV should cleanly eliminate the CV tracking row from the console.

  console_output <- capture.output({
    res <- peakFilter(
      x = mock_peaks,
      y = mock_treatment,
      fNA = c('T', 0.5),
      fCV = FALSE,
      fRT = FALSE
    )
  })

  # Convert vector of text lines to a single string for pattern search
  full_log <- paste(console_output, collapse = "\n")

  # Assertions
  expect_match(full_log, "Removed by NA filter:")
  expect_no_match(full_log, "Removed by CV filter:")
  expect_match(full_log, "Final peak count:")
})
