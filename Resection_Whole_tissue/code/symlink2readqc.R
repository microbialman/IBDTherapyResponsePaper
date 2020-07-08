
src_dir <- "/gfs/work/kralbrecht/03_BRC3_Matthias_Sam_Zoe/data/fastq"
target_dir <- "/gfs/work/kralbrecht/03_BRC3_Matthias_Sam_Zoe/readqc"

# Both forward and reverse (show that 2nd read is not great for 3' mRNA)
src_fastq <- list.files(src_dir, "fastq")

for (base_fastq in src_fastq){
  file.symlink(
    file.path(src_dir, base_fastq),
    file.path(target_dir, base_fastq)
  )
}
