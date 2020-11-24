# Module_Analyses

This contains the WGCNA analyses used to define the modules in the discovery resection dataset, any analyses done using the modules in the discovery whole tissue RNA-Seq data, and the module replications in the therapy response datasets.

Added by [Matthew Jackson](https://github.com/microbialman)

## Running code

All the code to run the analyses and generate figures etc. are found in the R markdown reports (.Rmd).

The input data tables used for these analyses can be found in Data/Input.
After running output tables and plots will be found in the Output and Plots folders respectively.
Markdown reports will also provide a brief discussion of the analyses used in the .html files generated.

Session info for the discovery_module_analysis run used for the paper is given below to show the software versions used.


```r
R version 3.5.0 (2018-04-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS release 6.7 (Final)

Matrix products: default
BLAS: /gfs/apps/apps/R-3.5.0/lib64/R/lib/libRblas.so
LAPACK: /gfs/apps/apps/R-3.5.0/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
   [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
    [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
     [9] LC_ADDRESS=C               LC_TELEPHONE=C
     [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] WGCNA_1.68                  fastcluster_1.1.25
  [3] dynamicTreeCut_1.63-1       forcats_0.4.0
   [5] stringr_1.4.0               dplyr_0.8.4
    [7] purrr_0.3.3                 readr_1.3.1
     [9] tidyr_1.0.2                 tibble_2.1.3
     [11] tidyverse_1.3.0             ggpubr_0.2.5
     [13] magrittr_1.5                knitr_1.28
     [15] biomaRt_2.36.1              pheatmap_1.0.12
     [17] ungeviz_0.1.0               snow_0.4-3
     [19] DESeq2_1.20.0               SummarizedExperiment_1.10.1
     [21] DelayedArray_0.6.6          BiocParallel_1.14.2
     [23] matrixStats_0.55.0          Biobase_2.40.0
     [25] GenomicRanges_1.32.7        GenomeInfoDb_1.16.0
     [27] IRanges_2.14.12             S4Vectors_0.18.3
     [29] BiocGenerics_0.26.0         RColorBrewer_1.1-2
     [31] reshape2_1.4.3              cowplot_1.0.0
     [33] ggdendro_0.1-20             ggplot2_3.2.1
     [35] here_0.1

loaded via a namespace (and not attached):
  [1] colorspace_1.4-1       ggsignif_0.6.0         rprojroot_1.3-2
    [4] htmlTable_1.13.3       XVector_0.20.0         base64enc_0.1-3
      [7] fs_1.3.1               rstudioapi_0.11        bit64_0.9-7
       [10] mvtnorm_1.0-12         AnnotationDbi_1.42.1   fansi_0.4.1
        [13] lubridate_1.7.4        xml2_1.2.2             codetools_0.2-16
	 [16] splines_3.5.0          doParallel_1.0.15      robustbase_0.93-5
	  [19] impute_1.54.0          geneplotter_1.58.0     Formula_1.2-3
	   [22] jsonlite_1.6.1         broom_0.5.4            annotate_1.58.0
	    [25] GO.db_3.6.0            cluster_2.1.0          dbplyr_1.4.2
	     [28] rrcov_1.5-2            compiler_3.5.0         httr_1.4.1
	      [31] backports_1.1.5        assertthat_0.2.1       Matrix_1.2-18
	       [34] lazyeval_0.2.2         cli_2.0.1              acepack_1.4.1
	        [37] htmltools_0.4.0        prettyunits_1.1.1      tools_3.5.0
		 [40] gtable_0.3.0           glue_1.3.1             GenomeInfoDbData_1.1.0
		  [43] Rcpp_1.0.3             cellranger_1.1.0       vctrs_0.2.2
		   [46] preprocessCore_1.42.0  nlme_3.1-144           iterators_1.0.12
		    [49] xfun_0.12              strapgod_0.0.4.9000    rvest_0.3.5
		     [52] lifecycle_0.1.0        XML_3.99-0.3           DEoptimR_1.0-8
		      [55] zlibbioc_1.26.0        MASS_7.3-51.5          scales_1.1.0
		       [58] hms_0.5.3              memoise_1.1.0          gridExtra_2.3
		        [61] rpart_4.1-15           latticeExtra_0.6-28    stringi_1.4.5
			 [64] RSQLite_2.2.0          highr_0.8              genefilter_1.62.0
			  [67] pcaPP_1.9-73           foreach_1.4.8          checkmate_2.0.0
			   [70] rlang_0.4.4            pkgconfig_2.0.3        bitops_1.0-6
			    [73] evaluate_0.14          lattice_0.20-38        htmlwidgets_1.5.1
			     [76] robust_0.4-18.2        bit_1.1-15.2           tidyselect_1.0.0
			      [79] plyr_1.8.5             R6_2.4.1               fit.models_0.5-14
			       [82] generics_0.0.2         Hmisc_4.3-1            DBI_1.1.0
			        [85] pillar_1.4.3           haven_2.2.0            foreign_0.8-72
				 [88] withr_2.1.2            survival_3.1-8         RCurl_1.98-1.1
				  [91] nnet_7.3-12            modelr_0.1.5           crayon_1.3.4
				   [94] progress_1.2.2         locfit_1.5-9.1         grid_3.5.0
				    [97] readxl_1.3.1           data.table_1.12.8      blob_1.2.1
				    [100] reprex_0.3.0           digest_0.6.24          xtable_1.8-4
				    [103] munsell_0.5.0
```