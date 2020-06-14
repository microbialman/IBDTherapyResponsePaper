library(BiocParallel)
library(rmarkdown)

param <- BatchtoolsParam(
    workers=9, cluster="slurm", resources=list(ncpus=1, walltime=60*60*2, partition="compute", memory="16G"), # s*min*h
    registryargs=batchtoolsRegistryargs(packages=c("rmarkdown")),
    template="/gfs/home/kralbrecht/templates/slurm-3.5.tmpl", jobname="deseq2")

result <- bplapply(
    FUN=render,
    X=sprintf("deseq2_featurecounts.%02d.Rmd", 4:12),
    BPPARAM=param)
