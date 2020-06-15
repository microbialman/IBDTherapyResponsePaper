
src_dir <- "/gfs/work/kralbrecht/03_BRC3_Matthias_Sam_Zoe/data/fastq"
target_dir <- "/gfs/work/kralbrecht/03_BRC3_Matthias_Sam_Zoe/scrnaseq/data.dir"

# Only the forward read
src_fastq <- list.files(src_dir, "fastq.1.gz")

for (base_fastq in src_fastq){
  striped_fastq <- gsub(".1.gz$", ".gz", base_fastq)
  # We have multiple FASTQ files per sample, we make a folder per sample
  # that does not include the run index
  fastq_folder <- gsub("_[[:digit:]]+.fastq.gz", ".fastq.gz", striped_fastq)
  sample_folder <- file.path(target_dir, fastq_folder)
  if (!dir.exists(sample_folder)){
    dir.create(sample_folder)
  }
  message("sample folder: ", sample_folder)
  # Create the symlink
  message("symlink name: ", striped_fastq)
  file.symlink(
    file.path(src_dir, base_fastq),
    file.path(sample_folder, striped_fastq)
    )
}
