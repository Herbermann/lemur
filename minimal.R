Rcpp::compileAttributes()  # scans your cpp files and regenerates R/RcppExports.R
devtools::load_all()       # compiles everything and loads your package

#remotes::install_version("harmony", version = "2.0.3")

library("tidyverse")
library("SingleCellExperiment")
library("lemur")
set.seed(42)

data("glioblastoma_example_data", package = "lemur")

fit <- lemur(glioblastoma_example_data, design = ~ patient_id + condition, 
             n_embedding = 15, test_fraction = 0.5)



fit <- align_harmony(fit, max_iter = 100)