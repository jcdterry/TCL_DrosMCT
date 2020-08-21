# Data and Scripts for *The effect of natural enemies on the coexistence of competing species - an empirical test using Bayesian Modern Coexistence Theory*  (Terry, Chen & Lewis 2020)

## Contents:

- R markdown scripts for all analyses.

* `Fitting Core STAN models.rmd` Fits core STAN models and conducts LOO analysis. 
* `Main Text Analyses.rmd` Code for analyses and figures in main text 
* Scripts to generate appendices:
  + `SI1.rmd` Further Experimental Design Details
  + `SI2.rmd` Model Design - Allee Effect, Sex Ratio, Emergence Time Model
  + `SI3.rmd` Sensitivity to Priors
  + `SI4.rmd` Additional Results
  + `SI5.rmd` Fitting and Testing Count Model
  + `SI6.rmd` Intransitivity Analysis


Folders:

* `Scripts/` Source code for utility functions to write STAN models and analyse STAN outputs 
* `Data/` Raw data:
  + `d_both.csv` is the principal data file, arranged in long form with 1 row = total count of one species from 1 vial. It probably has too many columns, but represents the founder species in multiple ways for ease of use in different analyses
  + `Clean_IntraData.csv` and `CleanInterData.csv` include the same data but with 1 row = 1 vial. Where avaiable, counts are split by emergence time and sex. They also include data from the 'cold' (20 degree) treatment that yielded insufficent usable data for analysis.
  + `Data entry table.csv` includes the original layout of the vials within boxes in the incubator.
* `StanModels/`  Folder containing STAN models used to fit all models. Models are 'written' by the `Fitting Core STAN models.rmd` script, but are included here to avoid having to run the script to view the model. 
* `LOO/` Empty, Storage for outputs from `loo()` Too large for github, but will be archived at Zenodo. 
* `StanFits/` Empty, Storage for STAN fit objects used in main text analysis. Too large for github, but will be archived at Zenodo. 
* `StanFitsPriorTest/` Empty, Storage for STAN fit objects created while comparing different prior specifications 
* `Plots/` Empty, but would be populated by running markdown documents
* `Images/` Illustrative photos used in SI documents 


## Authorship and citation.


## Version information:

R version 3.5.3
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18362)

Other attached packages:

ggrepel_0.8.1
igraph_1.2.4
cowplot_0.9.4
loo_2.2.0
scales_1.0.0
knitr_1.23
shinystan_2.5.0
shiny_1.3.2
rstan_2.19.2
StanHeaders_2.19.0
forcats_0.4.0
stringr_1.4.0
dplyr_0.8.1
purrr_0.3.2
readr_1.3.1
tidyr_0.8.3
tibble_2.1.3
ggplot2_3.2.0
tidyverse_1.2.1
RevoUtils_11.0.3
RevoUtilsMath_11.0.0