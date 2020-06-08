report_title <- "BRC3"
report_author <- "Kevin Rue-Albrecht"

# location of the readqc_db
readqc_db="/gfs/work/kralbrecht/03_BRC3_Matthias_Sam_Zoe/readqc/csvdb"

# location of the scseq_db
scseq_db="/gfs/work/kralbrecht/03_BRC3_Matthias_Sam_Zoe/scrnaseq/csvdb"
name_field_titles=c("disease", "inflammation", "site", "gender", "id")
qc_name_field_titles=c(name_field_titles)
plot_groups = c("disease", "inflammation", "site", "gender")
experiment_groups = c("disease", "inflammation", "site", "gender")
plot_color="disease"
plot_shape="inflammation"
plot_label="id"
nreplicates=3

# exclude samples from analysis.
exclude <- unique(c(
  c("CD_I_LI_M_49", "CD_U_SI_F_7", "CD_U_SI_F_51", "CRC_N_LI_F_88",  "CRC_N_LI_M_62", "CRC_N_LI_M_70", "UC_I_LI_F_41"), # PCA
  c("CD_U_SI_F_7"), # GC > 50%
  c("OTHER_N_NA_NA_63"), # gender unknown
  c("OTHER_N_NA_NA_63", "UC_NA_LI_F_2"), # inflammed unknown
  c("CD_U_SI_F_7", "CD_U_SI_M_27", "OTHER_U_LI_M_39", "OTHER_I_LI_M_44", "OTHER_I_LI_F_53", "OTHER_U_SI_F_54", "OTHER_N_NA_NA_63", "OTHER_N_LI_M_79") # other samples to exclude: 7, 27, 39, 44, 53, 54, 63 and 79
))

# DESeq2 parameters
location_function <- "median"
fit_type <- "local"
# fit_type <- "parametric"

require(readxl)
xlsxFile <- "20181207_WholeTissue_RNAseq_metadata2.xlsx"
# excel_sheets(xlsxFile)
excelComparisons <- read.xlsx("20181207_WholeTissue_RNAseq_metadata2.xlsx", "comparisons_clean")
comparisonName <- "comparison_02" # This is the key!!!
experimentalDesignText <- unique(subset(excelComparisons, ContrastField == comparisonName, "ControlFor", drop=TRUE))
stopifnot(length(experimentalDesignText) == 1)
experimental_design <- paste0("~", gsub(",", "+", experimentalDesignText), "+", comparisonName)

contrasts <- list()
reference_levels <- list()
for (idx in which(excelComparisons$ContrastField == comparisonName)) {
    # Set up the constrasts
    contrastField <- excelComparisons$ContrastField[idx]
    contrastNumerator <- excelComparisons$Numerator[idx]
    contrastDenominator <- excelComparisons$Denominator[idx]
    contrastName <- paste(contrastField, contrastNumerator, "vs", contrastDenominator, sep = "_")
    contrasts[[contrastName]] <- c(contrastField, contrastNumerator, contrastDenominator)
    # Set up the reference levels
    reference_levels[[contrastField]] <- "A"
}

p_value_threshold <- 0.05
abs_fc_threshold <- 1.5
n_interesting_genes <- 15
hm_row_cex <- 0.8
hm_col_cex <- 0.8
