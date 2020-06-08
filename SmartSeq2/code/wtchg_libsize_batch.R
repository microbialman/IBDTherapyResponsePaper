library(forcats)
library(ggplot2)
library(cowplot)

x <- read.table("/gfs/work/kralbrecht/18_matthias/data/bsg-ftp.well.ox.ac.uk/wtchg_table.txt", as.is = TRUE)
colnames(x) <- c("wtchg", "sample", "count_text")
head(x)
dim(x)

x$count <- as.numeric(gsub("Count:", "", x$count_text))
x$lane <- factor(gsub("WTCHG_([[:digit:]]+)_[[:digit:]]+", "\\1", x$wtchg))

x$batch <- fct_collapse(x$lane,
    batch1 = "739087",
    batch2 = c("743310", "743311", "743312", "743313")
)

ggplot(x) +
    geom_jitter(aes(interaction(batch, lane), count)) +
    theme_cowplot()
