# Pakula_et_al_23
Repository assosiated with Pakula et al., 2023 'Involvment of ROS in adaptive reprogramming of neonatal glia after cerebellar neutron injury'.


# R SessionInfo for Ascl1_DESeq2.Rmd, GCP_DESeq2.Rmd, and Hopx_DESeq2.Rmd:
R version 4.2.0 (2022-04-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] rstatix_0.7.1                                   SeuratDisk_0.0.0.9020                           ggraph_2.1.0                                   
 [4] scran_1.24.1                                    scater_1.24.0                                   scuttle_1.6.3                                  
 [7] SingleCellExperiment_1.18.1                     miloR_1.4.0                                     edgeR_3.38.4                                   
[10] limma_3.52.4                                    TxDb.Hsapiens.UCSC.hg38.knownGene_3.15.0        org.Hs.eg.db_3.15.0                            
[13] sctransform_0.3.5                               forcats_0.5.2                                   purrr_0.3.5                                    
[16] readr_2.1.3                                     tidyverse_1.3.2                                 org.Mm.eg.db_3.15.0                            
[19] goseq_1.48.0                                    geneLenDataBase_1.32.0                          BiasedUrn_2.0.8                                
[22] EnhancedVolcano_1.14.0                          ggrepel_0.9.1                                   DESeq2_1.36.0                                  
[25] SummarizedExperiment_1.26.1                     MatrixGenerics_1.8.1                            matrixStats_0.62.0                             
[28] tibble_3.1.8                                    BSgenome.cellranger.arc.mm10.2020.A.2.0.0_2.0.0 chromVAR_1.18.0                                
[31] BSgenome.Mmusculus.UCSC.mm10_1.4.3              BSgenome_1.64.0                                 rtracklayer_1.56.1                             
[34] Biostrings_2.64.1                               XVector_0.36.0                                  TFBSTools_1.34.0                               
[37] JASPAR2020_0.99.10                              JASPAR2022_0.99.7                               BiocFileCache_2.4.0                            
[40] dbplyr_2.3.0                                    stringr_1.4.1                                   dplyr_1.0.10                                   
[43] tidyr_1.2.1                                     future_1.28.0                                   patchwork_1.1.2                                
[46] ggplot2_3.4.0                                   EnsDb.Mmusculus.v79_2.99.0                      ensembldb_2.20.2                               
[49] AnnotationFilter_1.20.0                         GenomicFeatures_1.48.4                          AnnotationDbi_1.58.0                           
[52] Biobase_2.56.0                                  sp_1.5-0                                        SeuratObject_4.1.2                             
[55] Seurat_4.1.0                                    Signac_1.8.0                                    GenomicRanges_1.48.0                           
[58] GenomeInfoDb_1.32.4                             IRanges_2.30.1                                  S4Vectors_0.34.0                               
[61] BiocGenerics_0.42.0                            

loaded via a namespace (and not attached):
  [1] rsvd_1.0.5                  ica_1.0-3                   RcppRoll_0.3.0              Rsamtools_2.12.0            lmtest_0.9-40               crayon_1.5.2               
  [7] spatstat.core_2.4-4         MASS_7.3-58.1               nlme_3.1-160                backports_1.4.1             reprex_2.0.2                rlang_1.0.6                
 [13] ROCR_1.0-11                 readxl_1.4.1                irlba_2.3.5.1               filelock_1.0.2              BiocParallel_1.30.3         rjson_0.2.21               
 [19] CNEr_1.32.0                 bit64_4.0.5                 glue_1.6.2                  poweRlaw_0.70.6             parallel_4.2.0              vipor_0.4.5                
 [25] spatstat.sparse_2.1-1       spatstat.geom_2.4-0         haven_2.5.1                 tidyselect_1.2.0            fitdistrplus_1.1-8          XML_3.99-0.11              
 [31] zoo_1.8-11                  GenomicAlignments_1.32.1    xtable_1.8-4                magrittr_2.0.3              evaluate_0.17               cli_3.4.1                  
 [37] zlibbioc_1.42.0             rstudioapi_0.14             miniUI_0.1.1.1              rpart_4.1.16                fastmatch_1.1-3             shiny_1.7.2                
 [43] BiocSingular_1.12.0         xfun_0.33                   cluster_2.1.4               caTools_1.18.2              tidygraph_1.2.2             KEGGREST_1.36.3            
 [49] listenv_0.8.0               TFMPvalue_0.0.8             png_0.1-7                   withr_2.5.0                 bitops_1.0-7                ggforce_0.4.1              
 [55] plyr_1.8.7                  cellranger_1.1.0            pracma_2.4.2                dqrng_0.3.0                 pillar_1.8.1                cachem_1.0.6               
 [61] fs_1.5.2                    hdf5r_1.3.7                 DelayedMatrixStats_1.18.1   vctrs_0.5.2                 ellipsis_0.3.2              generics_0.1.3             
 [67] tools_4.2.0                 beeswarm_0.4.0              munsell_0.5.0               tweenr_2.0.2                DelayedArray_0.22.0         fastmap_1.1.0              
 [73] compiler_4.2.0              abind_1.4-5                 httpuv_1.6.6                plotly_4.10.0               rgeos_0.5-9                 GenomeInfoDbData_1.2.8     
 [79] gridExtra_2.3               lattice_0.20-45             deldir_1.0-6                utf8_1.2.2                  later_1.3.0                 jsonlite_1.8.2             
 [85] scales_1.2.1                ScaledMatrix_1.4.1          pbapply_1.5-0               carData_3.0-5               sparseMatrixStats_1.8.0     genefilter_1.78.0          
 [91] lazyeval_0.2.2              promises_1.2.0.1            car_3.1-1                   R.utils_2.12.0              goftest_1.2-3               spatstat.utils_2.3-1       
 [97] reticulate_1.26             rmarkdown_2.17              cowplot_1.1.1               statmod_1.4.37              Rtsne_0.16                  uwot_0.1.14                
[103] igraph_1.3.5                survival_3.4-0              yaml_2.3.5                  htmltools_0.5.3             memoise_2.0.1               BiocIO_1.6.0               
[109] locfit_1.5-9.6              graphlayouts_0.8.4          viridisLite_0.4.1           digest_0.6.29               assertthat_0.2.1            mime_0.12                  
[115] rappdirs_0.3.3              RSQLite_2.2.18              future.apply_1.9.1          data.table_1.14.2           blob_1.2.3                  R.oo_1.25.0                
[121] splines_4.2.0               googledrive_2.0.0           ProtGenerics_1.28.0         RCurl_1.98-1.9              broom_1.0.2                 hms_1.1.2                  
[127] modelr_0.1.10               colorspace_2.0-3            ggbeeswarm_0.6.0            Rcpp_1.0.9                  RANN_2.6.1                  fansi_1.0.3                
[133] tzdb_0.3.0                  parallelly_1.32.1           R6_2.5.1                    grid_4.2.0                  ggridges_0.5.4              lifecycle_1.0.3            
[139] bluster_1.6.0               curl_4.3.3                  googlesheets4_1.0.1         leiden_0.4.3                Matrix_1.5-1                RcppAnnoy_0.0.19           
[145] RColorBrewer_1.1-3          htmlwidgets_1.5.4           beachmat_2.12.0             polyclip_1.10-0             biomaRt_2.52.0              timechange_0.2.0           
[151] seqLogo_1.62.0              rvest_1.0.3                 mgcv_1.8-40                 globals_0.16.1              spatstat.random_2.2-0       progressr_0.11.0           
[157] codetools_0.2-18            lubridate_1.9.0             GO.db_3.15.0                metapod_1.4.0               gtools_3.9.3                prettyunits_1.1.1          
[163] R.methodsS3_1.8.2           gtable_0.3.1                DBI_1.1.3                   tensor_1.5                  httr_1.4.4                  KernSmooth_2.23-20         
[169] stringi_1.7.8               progress_1.2.2              reshape2_1.4.4              farver_2.1.1                annotate_1.74.0             viridis_0.6.2              
[175] DT_0.27                     xml2_1.3.3                  BiocNeighbors_1.14.0        restfulr_0.0.15             geneplotter_1.74.0          scattermore_0.8            
[181] bit_4.0.4                   spatstat.data_2.2-0         pkgconfig_2.0.3             gargle_1.2.1                DirichletMultinomial_1.38.0 knitr_1.40  
