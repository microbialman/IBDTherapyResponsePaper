report_title <- "RNA-seq analysis of IFT88 knock-down timecourse"
report_author <- "Stephen Sansom"

# location of the readqc_db
readqc_db="/gfs/work/kralbrecht/16_matthias/readqc/csvdb"

# location of the scseq_db
name_field_titles=c("treatment","replicate","run")
qc_name_field_titles=c(name_field_titles)
plot_groups = c("treatment")
experiment_groups = c("treatment")
plot_color="treatment"
plot_shape="treatment"
plot_label="replicate"
nreplicates=3

# exclude samples from analysis.
exclude <- c()