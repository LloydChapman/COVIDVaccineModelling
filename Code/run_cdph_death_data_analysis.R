rm(list=ls())
source("/mnt/nlo_shared/code/analysis_functions.R")
# source("~/Dropbox/COVIDVaccineModelling/Code/analysis_functions.R")

dir <- "/mnt/nlo_shared/data/"
fnms <- list.files(dir,pattern = "processed_death_data_.*_1")
# dir <- "~/Dropbox/COVIDVaccineModelling/Data/"
# fnms <- list.files(dir,pattern = "processed_dummy_data_")

# Run Poisson regression
for (i in 1:length(fnms)){
  glm_fit <- run_death_regression(dir,fnms[i])
  saveRDS(glm_fit,file = paste0(dir,"death_regression_output",gsub("processed_death_data|\\.csv","",fnms[i]),".RDS"))
  # save(glm_fit,file = paste0(dir,"death_regression_dummy_output",gsub("processed_dummy_data|\\.csv","",fnms[i]),".RDS"))
}
