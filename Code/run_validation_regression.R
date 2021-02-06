rm(list=ls())

source("/mnt/nlo_shared/code/analysis_functions.R")

dir <- "/mnt/nlo_shared/data/"
deaths <- read.csv(paste0(dir,"processed_death_data_2020-02-05_3.csv"))

deaths$date_of_death <- as.Date(deaths$date_of_death)

# Subset data to fitting data
vldtn_end_date <- as.Date("2020-09-30")
deaths_fit <- deaths[deaths$date_of_death <= vldtn_end_date,]
deaths_test <- deaths[deaths$date_of_death > vldtn_end_date,]
fnm <- paste0("processed_death_data_2020-02-05_",vldtn_end_date,"_3.csv")
write.csv(deaths_fit,paste0(dir,fnm),row.names = F)
fnm1 <- paste0("processed_death_data_",vldtn_end_date,"_",max(deaths_test$date_of_death),"_3.csv")
write.csv(deaths_test,paste0(dir,fnm1),row.names = F)

# Run regression on fitting data
glm_fit <- run_death_regression(dir,fnm)
saveRDS(glm_fit,file = paste0(dir,"death_regression_output",gsub("processed_death_data|\\.csv","",fnm),".RDS"))

# Aggregate and save testing data
agg_deaths <- aggregate(n_deaths ~ county_res + age_cat + sex + race_ethnicity,deaths_test,sum)
write.csv(agg_deaths,paste0(dir,"agg_deaths",sub("processed_death_data|\\.csv","",fnm1)),row.names = F)