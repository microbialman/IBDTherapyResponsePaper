
require(readxl)

excel_data <- read_excel(path = "/gfs/work/kralbrecht/16_matthias/data/WTCHG_info.xlsx")

data_dir <- "/gfs/work/kralbrecht/16_matthias/data"

data_fastq_dir <- file.path(data_dir, "fastqs")
dir.create(data_fastq_dir)

archive_fastq_dir <- c(
  "/gfs/archive/powrie/proj036/bsg-ftp.well.ox.ac.uk/190605_K00181_0145_AH5GMCBBXY",
  "/gfs/archive/powrie/proj036/bsg-ftp.well.ox.ac.uk/190611_K00181_0147_AH5J3VBBXY"
)

# Notes:
# 1. FASTQ files themselves should be '-' separated for pipeline_readqc
# 2. Folders representing samples should be '_' separated for pipeline_scrnaseq
# 3. Avoid '_' in FASTQ basenames to avoid issues with csvdb

for (i in seq_len(nrow(excel_data))) {
  i_sample_name <- excel_data[i, "Sample_name", drop=TRUE]
  i_wtchg_id <- excel_data[i, "WTCHG_ID", drop=TRUE]
  message(i_sample_name)
  data_sample_dir <- file.path(data_fastq_dir, sprintf("%s.fastq.1.gz", i_sample_name))
  dir.create(data_sample_dir)
  i_archive_fastqs <- list.files(path = archive_fastq_dir, pattern = sprintf("^%s.*fastq.gz$", i_wtchg_id), full.names = TRUE)
  i_new_basenames <- gsub("_", "", gsub("_([12]).fastq.gz", ".fastq.\\1.gz", basename(i_archive_fastqs)))
  file.symlink(from = i_archive_fastqs, to = file.path(data_sample_dir, paste(gsub("_", "-", i_sample_name), i_new_basenames, sep = "-")))
}
