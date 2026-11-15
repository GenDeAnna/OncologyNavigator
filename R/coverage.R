# coverage.R
# Billing route, typical prior-auth pattern and a document checklist.
# Typical patterns for teaching purposes - always confirm with the payer.

coverage_table <- function(matches, payer_rules, payer) {
  if (nrow(matches) == 0) return(data.frame())
  pr <- payer_rules[payer_rules$payer_type == payer, ]
  idx <- match(matches$benefit, pr$benefit)
  data.frame(
    Therapy   = matches$drug,
    Benefit   = matches$benefit,
    `HCPCS code` = ifelse(is.na(matches$hcpcs), "Not applicable (oral)", matches$hcpcs),
    `Prior auth (typical)` = pr$pa_typical[idx],
    `Payer note` = pr$payer_note[idx],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

pa_checklist <- function(matches, line, stage_setting) {
  items <- c(
    "Pathology report confirming NSCLC and histology (squamous or nonsquamous)",
    "Stage documentation (imaging reports and stage)",
    "Biomarker report showing the alteration the drug label requires (lab, test name, date)",
    "ICD-10 diagnosis code (C34.x, malignant neoplasm of bronchus and lung)",
    "Performance status (ECOG) and relevant labs",
    "Treatment plan: drug, dose, schedule and planned start date"
  )
  if (nrow(matches) == 0) return(items)
  if (any(matches$benefit %in% c("Medical (Part B)", "Mixed"))) {
    items <- c(items, "HCPCS code and NDC for each infused drug, plus site of care (office or hospital outpatient)")
  }
  if (any(matches$benefit %in% c("Pharmacy (Part D)", "Mixed"))) {
    items <- c(items, "Specialty pharmacy referral, and a formulary exception request if the drug is not on formulary")
  }
  if (stage_setting == "metastatic" && line == "Subsequent") {
    items <- c(items, "Prior therapy history with dates and the reason for change (progression or side effects)")
  }
  if (any(matches$accelerated_approval == "Yes", na.rm = TRUE)) {
    items <- c(items, "At least one option has accelerated approval; some payers ask for extra clinical justification")
  }
  items
}

biomarker_test_note <- function(payer) {
  switch(payer,
    "Original Medicare" = "Medicare's national coverage determination for NGS (NCD 90.2) covers FDA-approved or FDA-cleared NGS companion diagnostic tests for patients with advanced cancer who have not had the same test for the same cancer. Prior authorization is generally not used.",
    "Medicare Advantage" = "The plan may require prior authorization for NGS testing. Under CMS-0057-F, plans must decide within 72 hours for expedited requests and 7 calendar days for standard requests for medical items and services such as this test.",
    "Medicaid" = "Coverage and prior authorization vary by state. CMS-0057-F decision timeframes (72 hours expedited, 7 calendar days standard) apply to Medicaid medical items and services such as this test.",
    "Commercial" = "Coverage varies by plan, and many states have biomarker testing coverage laws. CMS-0057-F timeframes apply only to qualified health plans on the federal exchange.",
    "Check the payer's policy."
  )
}
