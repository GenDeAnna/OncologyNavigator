# variants.R
# Reads an NGS results file in MAF format and turns each variant into a
# biomarker ID that the rest of the app understands (for example "EGFR_COMMON").

read_maf <- function(path) {
  maf <- readr::read_tsv(path, comment = "#",
                         col_types = readr::cols(.default = "c"),
                         show_col_types = FALSE)
  needed  <- c("Hugo_Symbol", "Variant_Classification", "HGVSp_Short")
  missing <- setdiff(needed, names(maf))
  if (length(missing) > 0) {
    stop("The file is missing these MAF columns: ", paste(missing, collapse = ", "))
  }
  if (!"Tumor_Sample_Barcode" %in% names(maf)) {
    maf$Tumor_Sample_Barcode <- "uploaded_sample"
  }
  maf
}

# Pull the first number out of a protein change, e.g. "p.E746_A750del" -> 746
first_codon <- function(hgvsp) {
  suppressWarnings(as.integer(stringr::str_extract(hgvsp, "[0-9]+")))
}

# Decide which biomarker (if any) each variant represents. NA = not actionable.
classify_variant <- function(gene, classification, hgvsp) {
  codon <- first_codon(hgvsp)
  hgvsp <- ifelse(is.na(hgvsp), "", hgvsp)
  erbb2_hotspots <- c("p.S310F", "p.S310Y", "p.L755S", "p.V777L", "p.V842I", "p.G776V")
  met_exon14 <- stringr::str_detect(tolower(hgvsp), "exon ?14")
  
  dplyr::case_when(
    gene == "EGFR" & hgvsp == "p.L858R" ~ "EGFR_COMMON",
    gene == "EGFR" & classification == "In_Frame_Del" & codon >= 729 & codon <= 761 ~ "EGFR_COMMON",
    gene == "EGFR" & classification == "In_Frame_Ins" & codon >= 762 & codon <= 823 ~ "EGFR_EX20INS",
    gene == "KRAS" & hgvsp == "p.G12C" ~ "KRAS_G12C",
    gene == "BRAF" & hgvsp == "p.V600E" ~ "BRAF_V600E",
    gene == "MET" & (classification == "Splice_Site" | met_exon14) ~ "MET_EX14",
    gene == "ERBB2" & (classification == "In_Frame_Ins" | hgvsp %in% erbb2_hotspots) ~ "ERBB2_MUT",
    classification == "Fusion" & gene == "ALK" ~ "ALK_FUSION",
    classification == "Fusion" & gene == "ROS1" ~ "ROS1_FUSION",
    classification == "Fusion" & gene == "RET" ~ "RET_FUSION",
    classification == "Fusion" & gene %in% c("NTRK1", "NTRK2", "NTRK3") ~ "NTRK_FUSION",
    classification == "Fusion" & gene == "NRG1" ~ "NRG1_FUSION",
    TRUE ~ NA_character_
  )
}

# Add a biomarker_id column and keep only the columns we show.
annotate_variants <- function(maf) {
  data.frame(
    sample         = maf$Tumor_Sample_Barcode,
    gene           = maf$Hugo_Symbol,
    classification = maf$Variant_Classification,
    protein_change = maf$HGVSp_Short,
    biomarker_id   = classify_variant(maf$Hugo_Symbol, maf$Variant_Classification, maf$HGVSp_Short),
    stringsAsFactors = FALSE
  )
}