report_title <- "RNA-seq analysis of IFT88 knock-down timecourse"
report_author <- "Stephen Sansom"

# location of the readqc_db
readqc_db="/gfs/work/ssansom/angus_ift88/readqc/csvdb"

# location of the scseq_db
scseq_db="/gfs/work/kralbrecht/18_matthias/scrnaseq-full/csvdb"
name_field_titles=c("biosample", "celltype")
qc_name_field_titles=c(name_field_titles)
plot_groups = c("biosample","celltype")
experiment_groups = c("celltype")
plot_color="celltype"
plot_shape="celltype"
plot_label="biosample"
nreplicates=3

# exclude samples from analysis.
exclude <- c(
    "1_E", "12_M", "16_M", "23_E", "26_N", "3_E", "4_E", "5_E",  "6_E", "7_E", "8_E", "9_E", # QC; see "outliers_identify.R"
    "26N", "15N", "12M", "22M", "16M")

# DESeq2 parameters
location_function <- "median"
fit_type <- "local"
experimental_design <- "~experiment_group"
reference_levels <- list()
contrasts <- list(
    # comparison_01 = c("comparison_01", "B", "A"),
    comparison_02 = c("comparison_02", "B", "A"),
    comparison_03 = c("comparison_03", "B", "A"),
    comparison_04 = c("comparison_04", "B", "A"),
    comparison_05 = c("comparison_05", "B", "A"),
    comparison_06 = c("comparison_06", "B", "A"),
    comparison_07 = c("comparison_07", "B", "A"),
    comparison_08 = c("comparison_08", "B", "A"),
    comparison_09 = c("comparison_09", "B", "A"),
    comparison_10 = c("comparison_10", "B", "A"),
    comparison_11 = c("comparison_11", "B", "A"),
    comparison_12 = c("comparison_12", "B", "A"),
    comparison_13 = c("comparison_13", "B", "A"),
    comparison_14 = c("comparison_14", "B", "A"),
    comparison_15 = c("comparison_15", "B", "A"),
    comparison_16 = c("comparison_16", "B", "A")
)

p_value_threshold <- 0.05
abs_fc_threshold <- 1.5
n_interesting_genes <- 15
hm_row_cex <- 0.8
hm_col_cex <- 0.8