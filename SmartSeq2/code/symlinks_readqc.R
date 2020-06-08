
# Batch 1 ----

map_file <- "/gfs/work/kralbrecht/18_matthias/data/bsg-ftp.well.ox.ac.uk/191119_K00181_0189_AHFTWYBBXY/sample_id_map.tsv"
map_table <- read.table(file = map_file, col.names = c("wtchg_id",  "sample_id"), as.is = TRUE)
head(map_table)

archive_folder <- "/gfs/work/kralbrecht/18_matthias/data/bsg-ftp.well.ox.ac.uk/191119_K00181_0189_AHFTWYBBXY"

map_table$fastq_1 <- sprintf("%s/%s%s", archive_folder, map_table$wtchg_id, "_1.fastq.gz")
head(map_table$fastq_1)
stopifnot(all(file.exists(map_table$fastq_1)))

map_table$fastq_2 <- sprintf("%s/%s%s", archive_folder, map_table$wtchg_id, "_2.fastq.gz")
head(map_table$fastq_2)
stopifnot(all(file.exists(map_table$fastq_2)))

readqc_folder <- "/gfs/work/kralbrecht/18_matthias/readqc-full"
stopifnot(dir.exists(readqc_folder))

map_table$sample_fields <- gsub("([[:digit:]]+)([[:alpha:]])", "\\1-\\2", map_table$sample_id)
map_table$sample_fields <- sprintf("%s-%s", map_table$sample_fields, gsub("WTCHG_([[:digit:]]+)_[[:digit:]]+", "\\1", map_table$wtchg_id))

map_table$symlink_1 <- sprintf("%s/%s%s", readqc_folder, map_table$sample_fields, ".fastq.1.gz")
map_table$symlink_2 <- sprintf("%s/%s%s", readqc_folder, map_table$sample_fields,  ".fastq.2.gz")

file.remove(list.files(readqc_folder, "fastq.1.gz", full.names = TRUE))
file.remove(list.files(readqc_folder, "fastq.2.gz", full.names = TRUE))

stopifnot(all(file.symlink(from = map_table$fastq_1, to = map_table$symlink_1)))
stopifnot(all(file.symlink(from = map_table$fastq_2, to = map_table$symlink_2)))

# Batch 2 ----

map_file <- "/gfs/work/kralbrecht/18_matthias/data/bsg-ftp.well.ox.ac.uk/191211_K00181_0191_AHFV2NBBXY/sample_id_map.tsv"
map_table <- read.table(file = map_file, col.names = c("wtchg_id",  "sample_id"), as.is = TRUE)
head(map_table)

archive_folder <- "/gfs/work/kralbrecht/18_matthias/data/bsg-ftp.well.ox.ac.uk/191211_K00181_0191_AHFV2NBBXY"

map_table$fastq_1 <- sprintf("%s/%s%s", archive_folder, map_table$wtchg_id, "_1.fastq.gz")
head(map_table$fastq_1)
stopifnot(all(file.exists(map_table$fastq_1)))

map_table$fastq_2 <- sprintf("%s/%s%s", archive_folder, map_table$wtchg_id, "_2.fastq.gz")
head(map_table$fastq_2)
stopifnot(all(file.exists(map_table$fastq_2)))

readqc_folder <- "/gfs/work/kralbrecht/18_matthias/readqc-full"
stopifnot(dir.exists(readqc_folder))

map_table$sample_fields <- gsub("([[:digit:]]+)([[:alpha:]])", "\\1-\\2", map_table$sample_id)
map_table$sample_fields <- sprintf("%s-%s", map_table$sample_fields, gsub("WTCHG_([[:digit:]]+)_[[:digit:]]+", "\\1", map_table$wtchg_id))

map_table$symlink_1 <- sprintf("%s/%s%s", readqc_folder, map_table$sample_fields, ".fastq.1.gz")
map_table$symlink_2 <- sprintf("%s/%s%s", readqc_folder, map_table$sample_fields,  ".fastq.2.gz")

# file.remove(list.files(readqc_folder, "fastq.1.gz", full.names = TRUE))
# file.remove(list.files(readqc_folder, "fastq.2.gz", full.names = TRUE))

stopifnot(all(file.symlink(from = map_table$fastq_1, to = map_table$symlink_1)))
stopifnot(all(file.symlink(from = map_table$fastq_2, to = map_table$symlink_2)))
