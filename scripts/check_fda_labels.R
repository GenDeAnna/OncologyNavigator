# check_fda_labels.R
# Checks each rule in data/therapy_rules.csv against the drug's current FDA label
# (openFDA drug label API). Writes results/label_check.csv for review.

rules <- readr::read_csv("data/therapy_rules.csv",
                         col_types = readr::cols(.default = "c"),
                         show_col_types = FALSE)

# Words that should appear in the label's indication for each biomarker
keywords <- list(
  EGFR_COMMON  = c("exon 19", "L858R"),
  EGFR_EX20INS = c("exon 20 insertion"),
  ALK_FUSION   = c("\\bALK\\b"),
  ROS1_FUSION  = c("ROS1"),
  BRAF_V600E   = c("V600E"),
  KRAS_G12C    = c("G12C"),
  MET_EX14     = c("exon 14"),
  RET_FUSION   = c("\\bRET\\b"),
  NTRK_FUSION  = c("NTRK"),
  ERBB2_MUT    = c("ERBB2"),
  NRG1_FUSION  = c("NRG1"),
  PDL1         = c("PD-L1")
)

get_label <- function(brand) {
  query <- utils::URLencode(paste0('openfda.brand_name:"', brand, '"'), reserved = TRUE)
  url   <- paste0("https://api.fda.gov/drug/label.json?search=", query, "&limit=1")
  res   <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
  if (is.null(res)) return(list(found = FALSE, text = "", date = NA))
  r <- res$results
  text <- ""
  if (!is.null(r$indications_and_usage)) text <- paste(unlist(r$indications_and_usage), collapse = " ")
  list(found = TRUE, text = text, date = r$effective_time[1])
}

# Use the first brand named, e.g. "Tafinlar + Mekinist" -> "Tafinlar"
first_brand <- trimws(sub("[+].*", "", rules$brand))

cache  <- list()
report <- data.frame()

for (i in seq_len(nrow(rules))) {
  b <- first_brand[i]
  if (is.null(cache[[b]])) {
    cat("Checking", b, "...\n")
    cache[[b]] <- get_label(b)
    Sys.sleep(0.3)   # be polite to the free API
  }
  lab  <- cache[[b]]
  pats <- keywords[[rules$biomarker_id[i]]]
  
  hits <- sapply(pats, function(p) grepl(p, lab$text, ignore.case = TRUE))
  mentions_biomarker <- any(hits)
  mentions_nsclc     <- grepl("non-small cell lung|NSCLC", lab$text, ignore.case = TRUE)
  
  # Short excerpt around the first matching keyword
  excerpt <- ""
  if (mentions_biomarker) {
    pos <- regexpr(pats[which(hits)[1]], lab$text, ignore.case = TRUE)
    excerpt <- substr(lab$text, max(1, pos - 150), pos + 250)
  }
  
  report <- rbind(report, data.frame(
    rule_id            = rules$rule_id[i],
    drug               = rules$drug[i],
    brand_searched     = b,
    biomarker          = rules$biomarker_id[i],
    label_found        = lab$found,
    label_date         = lab$date,
    mentions_nsclc     = mentions_nsclc,
    mentions_biomarker = mentions_biomarker,
    status = ifelse(lab$found & mentions_nsclc & mentions_biomarker,
                    "PASS - skim excerpt", "REVIEW BY HAND"),
    excerpt            = excerpt
  ))
}

dir.create("results", showWarnings = FALSE)
readr::write_csv(report, "results/label_check.csv")

cat("\n===== LABEL CHECK SUMMARY =====\n")
print(table(report$status))
cat("\nRules needing review by hand:\n")
print(report[report$status == "REVIEW BY HAND", c("rule_id", "drug", "brand_searched", "label_found")])

View(report)
