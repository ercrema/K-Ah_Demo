# R scripts for the paper 'Divergent Demographic Responses to the Kikai-Akahoya Volcanic Catastrophe in Mid-Holocene Japan'

This repository contains data and scripts used in the following paper:

Crema, E.R., Kuwahata, M., Uchiyama, J., Junno, A., Riede, F, and Jordan, P.D. (submitted). _Divergent Demographic Responses to the Kikai-Akahoya Volcanic Catastrophe in Mid-Holocene Japan_

### External Data ###
External data required for all analyses are either dynamically downloaded or contained in the _data_ folder. Radiocarbon dates for demographic inference were obtained from the [Database of Radiocarbon Dates Published in Japanese Archaeological Site Reports](https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd_en/db_param#:~:text=Thus%2C%20a%20systematic%20collation%20of,continually%20being%20maintained%20and%20updated.), curated by the National Museum of Japanese History. We used version 1.2.0 of the curated English translation of the database (see [Kudo et al 2023](https://doi.org/10.5334/joad.115) for details), dynamically downloaded within the R script file `01_dataprep.R`. Ashfall data for the Kikai-Akahoya eruption were obtained from the [Tephra Database](https://doi.org/10.5281/zenodo.5109160) described in [Uesawa et al 2022](https://doi.org/10.1186/s13617-022-00126-x). 

### Main Pipeline ###
The analyses are conducted in three steps: data preparation, core analysis, and post-processing. Data preparation is executed in the R script file `01_dataprep.R`, where external files are imported and pre-processed for the core analyses. Please note that executing `01_dataprep.R` requires an internet connection. The script generates three R image files (`01_dataprep_out_d500.RData`, `01_dataprep_out_d750.RData`, and `01_dataprep_out_d1000.RData`) as output, containing all the necessary objects for subsequent analysis and visualisation of the results across three different time windows. The core analyses consist of two series of Bayesian models, fitted using the NIMBLE package, and are included in the folder `runscripts`. The naming convention distinguished between the two models (`02a*` and `02b*`), the variant (`icar0` for ICAR only, `icar1` with distance as predictor, `icar2` as the changepoint model presented in the main manuscript, and `icar3` as the ash thickness as predictor), and the time window (`500`, `750`, and `1000`, with the three version only for the `icar2` model). These scripts require high computational costs and multicore processors (at least 4 cores) with an estimated runtime between 12 and 24 hours, depending on the CPU. The outputs of the two models are stored as R image files in the `Results` folder. Finally, the R script file `03_processresults.R` contains functions and commands to generate all figures and tables for both the main and supplementary materials (see section below for details). The script `0X_figure_chronology.R` can be executed independently, as it generates a chart detailing the chrono-typological sequence on the island of Kyushu, which is not required for any of the analyses.  


### Outputs in Relation to the Research Article ###

All output figures and tables generated using the R script `03_processresults.R` are contained in the directory `figures_and_tables`. Correspondence to figure and table numbers in the manuscript can be found in the table below:

| **This repository**           | **Paper** | **Type** | **Location** |
|-------------------------------|-----------|----------|--------------|
| figure_chronology.pdf         | Figure 2  | Figure   | Main         |
| samplemap_and_ashfall.pdf     | Figure 3  | Figure   | Main         |
| mean_posterior.pdf            | Figure 4  | Figure   | Main         |
| hex_focus_plot.pdf            | Figure 5  | Figure   | Main         |
| post_scatter_ashfall.pd       | Figure 6  | Figure   | Main         |
| beta_posterior.pdf            | Figure 7  | Figure   | ESM          |
| samplemap_500_and_1000.pdf    | Figure S1 | Figure   | ESM          |
| posterior_r1.pdf              | Figure S2 | Figure   | ESM          |
| posterior_r1.pdf              | Figure S3 | Figure   | ESM          |
| posterior_eta.pdf             | Figure S4 | Figure   | ESM          |
| decrease_probability.pdf      | Figure S5 | Table    | ESM          |
| mean_posterior_sensitivity.pdf| Figure S6 | Table    | ESM          |
| waic.csv                      | Table S1  | Table    | ESM          |
| waic.csv                      | Table S2  | Table    | ESM          |
| posteriors_model_a.csv        | Table S3  | Table    | ESM          |
| posteriors_model_b.csv        | Table S4  | Table    | ESM          |


Please note that Figure 1 in the main text was a photograph, while Figure 2 was edited in Inkscape.

### R Session Info ####
```
R version 4.5.1 (2025-06-13)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 20.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/London
tzcode source: system (glibc)

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] spdep_1.3-11        spData_2.3.4        nimbleCarbon_0.2.5 
 [4] nimble_1.3.0        rcarbon_1.5.1       ggridges_0.5.6     
 [7] tidyr_1.3.1         gridExtra_2.3       viridis_0.6.5      
[10] viridisLite_0.4.2   terra_1.8-50        latex2exp_0.9.6    
[13] tidybayes_3.0.7     rnaturalearth_1.0.1 dplyr_1.1.4        
[16] coda_0.19-4.1       sf_1.0-20           ggplot2_3.5.2      
[19] here_1.0.1         

loaded via a namespace (and not attached):
 [1] svUnit_1.0.6           tidyselect_1.2.1       farver_2.1.2          
 [4] pracma_2.4.4           spatstat.geom_3.3-6    spatstat.explore_3.4-2
 [7] tensorA_0.36.2.1       rpart_4.1.24           lifecycle_1.0.4       
[10] spatstat.data_3.1-6    magrittr_2.0.3         posterior_1.6.1       
[13] compiler_4.5.1         rlang_1.1.6            doSNOW_1.0.20         
[16] tools_4.5.1            igraph_2.1.4           knitr_1.50            
[19] sp_2.2-0               classInt_0.4-11        RColorBrewer_1.1-3    
[22] abind_1.4-8            KernSmooth_2.23-26     numDeriv_2016.8-1.1   
[25] withr_3.0.2            purrr_1.0.4            grid_4.5.1            
[28] polyclip_1.10-7        e1071_1.7-16           scales_1.4.0          
[31] iterators_1.0.14       spatstat.utils_3.1-3   spatstat_3.3-2        
[34] cli_3.6.5              generics_0.1.4         httr_1.4.7            
[37] DBI_1.2.3              proxy_0.4-27           stringr_1.5.1         
[40] splines_4.5.1          spatstat.model_3.3-5   s2_1.1.8              
[43] vctrs_0.6.5            boot_1.3-31            Matrix_1.7-3          
[46] jsonlite_2.0.0         arrayhelpers_1.1-0     tensor_1.5            
[49] ggdist_3.3.3           foreach_1.5.2          spatstat.univar_3.1-3 
[52] units_0.8-7            snow_0.4-4             goftest_1.2-3         
[55] glue_1.8.0             spatstat.random_3.3-3  codetools_0.2-19      
[58] distributional_0.5.0   stringi_1.8.7          gtable_0.3.6          
[61] deldir_2.0-4           tibble_3.2.1           pillar_1.10.2         
[64] R6_2.6.1               wk_0.9.4               rprojroot_2.0.4       
[67] evaluate_1.0.3         lattice_0.22-5         backports_1.5.0       
[70] class_7.3-23           Rcpp_1.0.14            spatstat.linnet_3.2-5 
[73] nlme_3.1-168           checkmate_2.3.2        spatstat.sparse_3.1-0 
[76] mgcv_1.9-1             xfun_0.52              pkgconfig_2.0.3  
```
### Funding ###
Funded by the Swedish Research Council (VR) project _Surviving the Apocalypse: multidimensional modeling of the impact of a prehistoric megadisaster on people’s lifeworlds, technologies and demography_ (Grant Code: 2024-00822), JSPS KAKENHI (Grant Numbers JP25K00523), and European Research Council Synergy project _FORAGER: Investigating alternative trajectories for human demographic growth in temperate northern Holocene societies (Grant Code: 101224035)


