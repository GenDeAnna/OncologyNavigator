# therapy_match.R
# Finds the rules in therapy_rules.csv that fit this patient.

match_therapies <- function(rules, status, stage_setting, patient_histology, line) {
  positive <- status$biomarker_id[status$result == "Positive"]
  pdl1     <- status$value[status$biomarker_id == "PDL1"]
  pdl1_val <- if (length(pdl1) == 1 && !is.na(pdl1)) pdl1 else -Inf

  # 1. Right setting, right histology, biomarker is positive
  r <- rules[rules$setting == stage_setting &
             histology_matches(rules$histology, patient_histology) &
             rules$biomarker_id %in% positive, ]

  # 2. Line of therapy (only matters for advanced disease)
  if (stage_setting == "metastatic" && line == "First-line") {
    r <- r[r$line != "Subsequent", ]
  }

  # 3. PD-L1 cut-offs
  r <- r[is.na(r$min_value) | pdl1_val >= r$min_value, ]

  # 4. Exclusions (e.g. immunotherapy alone not for EGFR- or ALK-positive tumors)
  excluded <- vapply(r$exclude_if, function(x) {
    if (is.na(x) || x == "") return(FALSE)
    any(strsplit(x, ";")[[1]] %in% positive)
  }, logical(1))
  r <- r[!excluded, ]

  r[order(r$biomarker_id, r$drug), ]
}

format_matches <- function(matches, catalog) {
  if (nrow(matches) == 0) return(data.frame())
  data.frame(
    Biomarker = catalog$display_name[match(matches$biomarker_id, catalog$biomarker_id)],
    Therapy   = matches$drug,
    Brand     = matches$brand,
    Route     = matches$route,
    `Approved line` = matches$line,
    `Accelerated approval` = ifelse(is.na(matches$accelerated_approval), "", matches$accelerated_approval),
    `Rule verified` = matches$last_verified,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
