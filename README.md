# R scripts for the paper 'Divergent Demographic Responses to the Kikai-Akahoya Mega-Eruption in Prehistoric Japan'

This repository contains data and scripts used in the following paper:

Crema, E.R. et al (In Prep). Divergent Demographic Responses to the Kikai-Akahoya Mega-Eruption in Prehistoric Japan

### External Data ###
External data required for all analyses are either dynamically downloaded or contained in the _data_ folder. Radiocarbon dates for demographic inference were obtained from the [Database of Radiocarbon Dates Published in Japanese Archaeological Site Reports](https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd_en/db_param#:~:text=Thus%2C%20a%20systematic%20collation%20of,continually%20being%20maintained%20and%20updated.), curated by the National Museum of Japanese History. We used version 1.2.0 of the curated English translation of the database (see [Kudo et al 2023](https://doi.org/10.5334/joad.115) for details), dynamically downloaded within the R script file `01_dataprep.R`. Ashfall data for the Kikai-Akahoya eruption were obtained from the [Tephra Database](https://doi.org/10.5281/zenodo.5109160) described in [Uesawa et al 2022](https://doi.org/10.1186/s13617-022-00126-x). 

### Main Pipeline ###
The analyses are conducted in three steps: data preparation, core analysis, and post-processing. Data preparation is executed in the R script file `01_dataprep.R`, where external files are imported and pre-processed for the core analyses. Please note that executing `01_dataprep.R` requires an internet connection. The script generates an R image file (`01_dataprep_out.RData`) as output, containing all the necessary objects for subsequent analysis and visualisation of the results. The core analyses consist of two Bayesian models fitted using the NIMBLE package and included in two separate R script files: `02a_fitmodel.R` and `02b_fitmodel.R`. Both scripts require high computational costs and multicore processors (at least 4 cores) with an estimated runtime between 12 and 24 hours, depending on the CPU. The output of the two models are stored in the R image files `/Results/02a_fitmodel_out.RData` and `/Results/02b_fitmodel_out.RData`. Finally, the R script file `03_plotresults.R` contains functions and commands to generate all figures and tables for both the main and supplementary materials (see section below for details).  


### Outputs in Relation to the Research Article ###


### File Structure ####

### 

