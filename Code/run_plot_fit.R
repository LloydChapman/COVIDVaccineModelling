rm(list=ls())

source("prediction_functions.R")

dir <- "../Data/"
setwd(dir)  
fnms <- list.files(dir,pattern = "regression_output_.*.RData")
deaths_ratio <- readRDS("deaths_multiplier.RDS")
  
for (i in 1:length(fnms)){
  plot_fit(dir,fnms[i],deaths_ratio[i])
}  

for (i in 1:length(fnms)){
  p_death <- readRDS(paste0("prob_death_given_clin_by_age_and_sex",gsub("regression_output|\\.RData","",fnms[4]),".RDS"))
  View(p_death)
}