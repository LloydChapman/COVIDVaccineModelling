library(broom)

dir <- "~/Dropbox/COVIDVaccineModelling/Data/"
fnms <- list.files(dir,pattern = "death_regression_output_.*_1.RDS")

for (i in 1:length(fnms)){
  glm_fit <- readRDS(paste0(dir,fnms[i]))
  sink(file = paste0(dir,sub("\\.RDS","\\.txt",fnms[i])))
  print(summary(glm_fit))
  sink()
  x <- tidy(glm_fit)
  x$exp_estimate <- exp(x$estimate)
  write.csv(x,paste0(dir,"coeffs_",gsub("_regression_output|\\.RDS","",fnms[i]),".csv"),row.names = F)
}