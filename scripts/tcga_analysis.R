# tcga_analysis.R
# Applies the Navigator's variant classifier to real tumors from
# TCGA Lung Adenocarcinoma (PanCancer Atlas 2018, via cBioPortal).
# Data: run the download steps first; raw files live in rawdata/ (not on GitHub).

source("R/variants.R")

# ---- Part 1: mutations ------------------------------------------------------
# Read only the four columns the classifier needs (the full file is large)
mut <- readr::read_tsv(
  "rawdata/luad_tcga_pan_can_atlas_2018/data_mutations.txt",
  comment = "#",
  col_select = c(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification, HGVSp_Short),
  col_types = readr::cols(.default = "c"),
  show_col_types = FALSE
)

ann <- annotate_variants(mut)
actionable_mut <- ann[!is.na(ann$biomarker_id), ]

n_tumors <- length(unique(mut$Tumor_Sample_Barcode))

cat("\n===== PART 1: MUTATIONS =====\n")
cat("Mutations read:            ", nrow(mut), "\n")
cat("Tumors with mutation data: ", n_tumors, "\n")
cat("Actionable mutations found:", nrow(actionable_mut), "\n\n")
cat("Tumors per actionable biomarker (mutations only):\n")
print(sort(tapply(actionable_mut$sample, actionable_mut$biomarker_id,
                  function(x) length(unique(x))), decreasing = TRUE))

# ---- Part 2: fusions --------------------------------------------------------
sv <- readr::read_tsv(
  "rawdata/luad_tcga_pan_can_atlas_2018/data_sv.txt",
  col_types = readr::cols(.default = "c"),
  show_col_types = FALSE
)

needed_cols <- c("Sample_Id", "Site1_Hugo_Symbol", "Site2_Hugo_Symbol")
if (!all(needed_cols %in% names(sv))) {
  print(names(sv))
  stop("The fusion file has different column names - send DeAnna's helper the list above.")
}

fusion_genes <- c(ALK = "ALK_FUSION", ROS1 = "ROS1_FUSION", RET = "RET_FUSION",
                  NTRK1 = "NTRK_FUSION", NTRK2 = "NTRK_FUSION", NTRK3 = "NTRK_FUSION",
                  NRG1 = "NRG1_FUSION")

# A fusion counts if either partner gene is one of the targets
gene_hit <- ifelse(sv$Site1_Hugo_Symbol %in% names(fusion_genes), sv$Site1_Hugo_Symbol,
                   ifelse(sv$Site2_Hugo_Symbol %in% names(fusion_genes), sv$Site2_Hugo_Symbol, NA))

actionable_fus <- data.frame(
  sample       = sv$Sample_Id[!is.na(gene_hit)],
  partner_1    = sv$Site1_Hugo_Symbol[!is.na(gene_hit)],
  partner_2    = sv$Site2_Hugo_Symbol[!is.na(gene_hit)],
  biomarker_id = unname(fusion_genes[gene_hit[!is.na(gene_hit)]])
)
# Keep only tumors that are also in the mutation data, so the denominator matches
actionable_fus <- actionable_fus[actionable_fus$sample %in% mut$Tumor_Sample_Barcode, ]

cat("\n===== PART 2: FUSIONS =====\n")
print(actionable_fus)

# ---- Part 3: headline -------------------------------------------------------
all_hits <- rbind(actionable_mut[, c("sample", "biomarker_id")],
                  actionable_fus[, c("sample", "biomarker_id")])

tumors_any <- length(unique(all_hits$sample))

cat("\n===== PART 3: HEADLINE =====\n")
cat("Tumors analyzed:                         ", n_tumors, "\n")
cat("Tumors with >= 1 actionable biomarker:   ", tumors_any, "\n")
cat("Percent with an FDA-approved targeted option:",
    round(100 * tumors_any / n_tumors, 1), "%\n")
cat("Tumors with more than one actionable hit:",
    sum(table(unique(all_hits)$sample) > 1), "\n")

# ---- Part 4: save results ---------------------------------------------------
dir.create("results", showWarnings = FALSE)

summary_tbl <- aggregate(sample ~ biomarker_id, data = unique(all_hits), FUN = length)
names(summary_tbl)[2] <- "tumors"
summary_tbl$percent <- round(100 * summary_tbl$tumors / n_tumors, 1)
catalog <- readr::read_csv("data/biomarker_catalog.csv", show_col_types = FALSE)
summary_tbl$biomarker <- catalog$display_name[match(summary_tbl$biomarker_id, catalog$biomarker_id)]
summary_tbl <- summary_tbl[order(-summary_tbl$tumors), ]
readr::write_csv(summary_tbl, "results/tcga_luad_actionable_summary.csv")

png("results/tcga_luad_actionable.png", width = 1800, height = 1100, res = 200)
par(mar = c(5, 17, 4, 2))
barplot(rev(summary_tbl$percent), names.arg = rev(summary_tbl$biomarker),
        horiz = TRUE, las = 1, col = "#2C7FB8", border = NA,
        xlab = "Percent of tumors", cex.main = 0.9,
        main = paste0("Actionable biomarkers in TCGA lung adenocarcinoma (n = ", n_tumors,
                      "; ", round(100 * tumors_any / n_tumors, 1), "% with any)"))
dev.off()

cat("\nSaved results/tcga_luad_actionable_summary.csv and results/tcga_luad_actionable.png\n")
print(summary_tbl)
