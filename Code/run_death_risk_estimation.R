rm(list=ls())

source("~/Dropbox/COVIDVaccineModelling/Code/analysis_functions.R")

dir <- "~/Dropbox/COVIDVaccineModelling/Data/"
fnms <- list.files(dir,pattern = "regression_output_.*_1.RDS")

for (i in 1:length(fnms)){
  glm_death_fit <- estimate_death_risk(dir,fnms[i])
  saveRDS(glm_death_fit,file = paste0(dir,"death_risk_regression_output",gsub("regression_output|\\.RDS","",fnms[i]),".RDS"))
  res <- calc_death_risk_probs(glm_death_fit)
  saveRDS(res,paste0("prob_death_given_clin_by_age_and_sex",gsub("regression_output|\\.RDS","",fnms[i]),".RDS"))
}


