
map_file <- "/gfs/work/kralbrecht/03_BRC3_Matthias_Sam_Zoe/data/symlinks.txt"
map_data <- read.table(map_file, as.is = TRUE)

archive_root <- "/gfs/archive/powrie/proj027"
data_dir <- "/gfs/work/kralbrecht/03_BRC3_Matthias_Sam_Zoe/data/fastq"

stopifnot(dir.create(data_dir))

for (row_index in seq_len(nrow(map_data))){
  # Fetch info
  sample_name <- map_data[row_index, "V3"]
  base_run1 <- map_data[row_index, "V1"]
  base_run2 <- map_data[row_index, "V2"]
  # Find source FASTQ files
  archive_run1 <- list.files(path = archive_root, pattern = base_run1, full.names = TRUE, recursive = TRUE)
  archive_run2 <- list.files(path = archive_root, pattern = base_run2, full.names = TRUE, recursive = TRUE)
  # Prepare symlinks filenames
  # Move the pairing information
  renamed_run1 <- gsub(pattern = "_([[:digit:]]).fastq.gz", ".fastq.\\1.gz", basename(archive_run1))
  renamed_run2 <- gsub(pattern = "_([[:digit:]]).fastq.gz", ".fastq.\\1.gz", basename(archive_run2))
  # Substitute sample names
  renamed_run1 <- gsub(pattern = base_run1, sample_name, basename(renamed_run1))
  renamed_run2 <- gsub(pattern = base_run2, sample_name, basename(renamed_run2))
  # Add run identifier
  renamed_run1 <- gsub(pattern = "(.fastq.*)", "_1\\1", basename(renamed_run1))
  renamed_run2 <- gsub(pattern = "(.fastq.*)", "_2\\1", basename(renamed_run2))
  # Generate the full symlink path
  symlink_run1 <- file.path(data_dir, renamed_run1)
  symlink_run2 <- file.path(data_dir, renamed_run2)
  # Make the symlinks
  file.symlink(archive_run1, symlink_run1)
  file.symlink(archive_run2, symlink_run2)
}

README <- cat(
"Filename fields are:
- disease
- inflammation_status
- site
- gender
- sample_number
- replicate
", file = file.path(data_dir, "README")
)
