report_title <- "RNA-seq analysis of IFT88 knock-down timecourse"
report_author <- "Stephen Sansom"

# location of the readqc_db
readqc_db="/gfs/work/kralbrecht/16_matthias/readqc/csvdb"

# location of the scseq_db
scseq_db="/gfs/work/kralbrecht/16_matthias/scrnaseq/csvdb"
name_field_titles=c("treatment","replicate")
qc_name_field_titles=c(name_field_titles, "lane")
plot_groups = c("treatment")
experiment_groups = c("treatment")
plot_color="treatment"
plot_shape="treatment"
plot_label="replicate"
nreplicates=3

# exclude samples from analysis.
exclude <- c(
    "Unstimulated_3",
    "IL1beta_6",
    "TNF_6",
    "IL1betaTNF_6"
)

# DESeq2 parameters
location_function <- "median"
fit_type <- "local"
experimental_design <- "~experiment_group"
reference_levels <- list(treatment = "Unstimulated", replicate = 1)
contrasts <- list(
    "IL1betaTNF_vsIL1beta" = c("experiment_group","IL1betaTNF", "IL1beta"),
    "IL1betaTNF_vsTNF" = c("experiment_group","IL1betaTNF", "TNF"),
    "IL1betaTNF_vsUnstimulated" = c("experiment_group","IL1betaTNF", "Unstimulated"),
    "IL1beta_vsTNF" = c("experiment_group","IL1beta", "TNF"),
    "IL1beta_vsUnstimulated" = c("experiment_group","IL1beta", "Unstimulated"),
    "TNF_vsUnstimulated" = c("experiment_group","TNF", "Unstimulated")
)

p_value_threshold <- 0.05
abs_fc_threshold <- 2
n_interesting_genes <- 15
hm_row_cex <- 0.8
hm_col_cex <- 0.8