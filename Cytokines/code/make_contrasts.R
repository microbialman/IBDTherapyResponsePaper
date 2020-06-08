
# This script was used to prepare the list of contrasts for deseq2_featurecounts.Rmd

# x = ordered list of treatments

for (i in seq(1, length(x) - 1)) {
  for (j in seq(i+1, length(x))) {
    # cat(i, j, "\n")
    cat(sprintf(
      "\"%s_vs%s\" = c(\"experiment_group\",\"%s\", \"%s\"),\n",
      x[i], x[j], x[i], x[j]
    ))
  }
}
