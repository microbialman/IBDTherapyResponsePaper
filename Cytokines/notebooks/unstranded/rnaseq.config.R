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
exclude <- c()

# DESeq2 parameters
location_function <- "median"
fit_type <- "local"
experimental_design <- "~experiment_group"
reference_levels <- list(treatment = "Unstimulated", replicate = 1)
contrasts <- list(
  "IL1betaOSM_vsIL1betaTNF" = c("experiment_group","IL1betaOSM", "IL1betaTNF"),
  "IL1betaOSM_vsOSMTNF" = c("experiment_group","IL1betaOSM", "OSMTNF"),
  "IL1betaOSM_vsIL1beta" = c("experiment_group","IL1betaOSM", "IL1beta"),
  "IL1betaOSM_vsOSM" = c("experiment_group","IL1betaOSM", "OSM"),
  "IL1betaOSM_vsTNF" = c("experiment_group","IL1betaOSM", "TNF"),
  "IL1betaOSM_vsUnstimulated" = c("experiment_group","IL1betaOSM", "Unstimulated"),
  "IL1betaTNF_vsOSMTNF" = c("experiment_group","IL1betaTNF", "OSMTNF"),
  "IL1betaTNF_vsIL1beta" = c("experiment_group","IL1betaTNF", "IL1beta"),
  "IL1betaTNF_vsOSM" = c("experiment_group","IL1betaTNF", "OSM"),
  "IL1betaTNF_vsTNF" = c("experiment_group","IL1betaTNF", "TNF"),
  "IL1betaTNF_vsUnstimulated" = c("experiment_group","IL1betaTNF", "Unstimulated"),
  "OSMTNF_vsIL1beta" = c("experiment_group","OSMTNF", "IL1beta"),
  "OSMTNF_vsOSM" = c("experiment_group","OSMTNF", "OSM"),
  "OSMTNF_vsTNF" = c("experiment_group","OSMTNF", "TNF"),
  "OSMTNF_vsUnstimulated" = c("experiment_group","OSMTNF", "Unstimulated"),
  "IL1beta_vsOSM" = c("experiment_group","IL1beta", "OSM"),
  "IL1beta_vsTNF" = c("experiment_group","IL1beta", "TNF"),
  "IL1beta_vsUnstimulated" = c("experiment_group","IL1beta", "Unstimulated"),
  "OSM_vsTNF" = c("experiment_group","OSM", "TNF"),
  "OSM_vsUnstimulated" = c("experiment_group","OSM", "Unstimulated"),
  "TNF_vsUnstimulated" = c("experiment_group","TNF", "Unstimulated")
)

p_value_threshold <- 0.05
abs_fc_threshold <- 2
n_interesting_genes <- 15
hm_row_cex <- 0.8
hm_col_cex <- 0.8