rm(list=ls())

library(broom)

dir <- "../Data/"
fnms <- list.files(dir,pattern = "^regression_output_.*_1.RDS")

for (i in 1:length(fnms)){
  glm_fit <- readRDS(paste0(dir,fnms[i]))
  sink(file = paste0(dir,sub("\\.RDS","\\.txt",fnms[i])))
  print(summary(glm_fit))
  sink()
  x <- tidy(glm_fit)
  x$exp_estimate <- exp(x$estimate)
  write.csv(x,paste0(dir,"coeffs",gsub("regression_output|\\.RDS","",fnms[i]),".csv"),row.names = F)
}

