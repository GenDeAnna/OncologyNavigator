# biomarkers.R
# Each biomarker's result (Positive / Negative / Not tested) and testing gaps.

setting_for_stage <- function(stage_group) {
  if (startsWith(stage_group, "Early")) "adjuvant" else "metastatic"
}

histology_matches <- function(rule_histology, patient_histology) {
  rule_histology == "any" | rule_histology == patient_histology
}

# "ALL" means every biomarker except PD-L1 (PD-L1 has its own column).
expand_tested <- function(tested_text, catalog) {
  if (is.na(tested_text) || tested_text == "") return(character(0))
  if (tested_text == "ALL") return(setdiff(catalog$biomarker_id, "PDL1"))
  trimws(strsplit(tested_text, ";")[[1]])
}

build_biomarker_status <- function(catalog, found_ids, tested_ids, pdl1_tps) {
  status <- data.frame(
    biomarker_id = catalog$biomarker_id,
    display_name = catalog$display_name,
    stringsAsFactors = FALSE
  )
  status$result <- ifelse(status$biomarker_id %in% found_ids, "Positive",
                   ifelse(status$biomarker_id %in% tested_ids, "Negative", "Not tested"))
  status$value <- NA_real_
  is_pdl1 <- status$biomarker_id == "PDL1"
  if (is.na(pdl1_tps)) {
    status$result[is_pdl1] <- "Not tested"
  } else {
    status$result[is_pdl1] <- if (pdl1_tps >= 1) "Positive" else "Negative"
    status$value[is_pdl1]  <- pdl1_tps
  }
  status
}

relevant_biomarkers <- function(rules, stage_setting, patient_histology) {
  keep <- rules$setting == stage_setting &
          histology_matches(rules$histology, patient_histology)
  unique(rules$biomarker_id[keep])
}

find_testing_gaps <- function(rules, catalog, stage_setting, patient_histology, tested_ids) {
  needed  <- relevant_biomarkers(rules, stage_setting, patient_histology)
  missing <- setdiff(needed, tested_ids)
  catalog[catalog$biomarker_id %in% missing, c("display_name", "test_method")]
}
