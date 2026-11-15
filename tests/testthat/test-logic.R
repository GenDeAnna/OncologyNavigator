# test-logic.R - automatic checks for the matching rules
# Run from the Console with: testthat::test_dir("tests/testthat")

source("../../R/variants.R")
source("../../R/biomarkers.R")
source("../../R/therapy_match.R")
source("../../R/coverage.R")

rules <- readr::read_csv("../../data/therapy_rules.csv",
                         col_types = readr::cols(.default = "c", min_value = "d"),
                         show_col_types = FALSE)
catalog <- readr::read_csv("../../data/biomarker_catalog.csv", show_col_types = FALSE)

test_that("variants are classified into the right biomarkers", {
  expect_equal(classify_variant("EGFR", "Missense_Mutation", "p.L858R"), "EGFR_COMMON")
  expect_equal(classify_variant("EGFR", "In_Frame_Del", "p.E746_A750del"), "EGFR_COMMON")
  expect_equal(classify_variant("EGFR", "In_Frame_Ins", "p.D770_N771insSVD"), "EGFR_EX20INS")
  expect_equal(classify_variant("KRAS", "Missense_Mutation", "p.G12C"), "KRAS_G12C")
  expect_equal(classify_variant("ALK", "Fusion", "EML4-ALK fusion"), "ALK_FUSION")
  expect_true(is.na(classify_variant("KRAS", "Missense_Mutation", "p.G12D")))
  expect_true(is.na(classify_variant("TP53", "Missense_Mutation", "p.R273H")))
})

test_that("missing tests are flagged as gaps", {
  tested <- c("EGFR_COMMON", "ALK_FUSION")
  gaps <- find_testing_gaps(rules, catalog, "metastatic", "Nonsquamous", tested)
  expect_true("ROS1 rearrangement" %in% gaps$display_name)
  expect_false("ALK rearrangement" %in% gaps$display_name)
})

test_that("EGFR-positive tumors get osimertinib but not immunotherapy alone", {
  tested <- setdiff(catalog$biomarker_id, "PDL1")
  status <- build_biomarker_status(catalog, "EGFR_COMMON", c(tested, "PDL1"), 80)
  m <- match_therapies(rules, status, "metastatic", "Nonsquamous", "First-line")
  expect_true("osimertinib" %in% m$drug)
  expect_false("pembrolizumab" %in% m$drug)
})

test_that("PD-L1 cut-offs are respected", {
  tested <- c(setdiff(catalog$biomarker_id, "PDL1"), "PDL1")
  low  <- build_biomarker_status(catalog, character(0), tested, 10)
  high <- build_biomarker_status(catalog, character(0), tested, 60)
  m_low  <- match_therapies(rules, low,  "metastatic", "Squamous", "First-line")
  m_high <- match_therapies(rules, high, "metastatic", "Squamous", "First-line")
  expect_true("pembrolizumab" %in% m_low$drug)
  expect_false("atezolizumab" %in% m_low$drug)
  expect_true("atezolizumab" %in% m_high$drug)
})

test_that("later-line drugs are hidden for first-line patients", {
  tested <- setdiff(catalog$biomarker_id, "PDL1")
  status <- build_biomarker_status(catalog, "KRAS_G12C", tested, NA)
  first <- match_therapies(rules, status, "metastatic", "Nonsquamous", "First-line")
  later <- match_therapies(rules, status, "metastatic", "Nonsquamous", "Subsequent")
  expect_equal(nrow(first), 0)
  expect_true("sotorasib" %in% later$drug)
})
