library(broom)

dir <- "../Data/"
fnms <- list.files(dir,pattern = "death_regression_output_.*_2.RDS")

y <- vector("list",length(fnms))
for (i in 1:length(fnms)){
  glm_fit <- readRDS(paste0(dir,fnms[i]))
  sink(file = paste0(dir,sub("\\.RDS","\\.txt",fnms[i])))
  print(summary(glm_fit))
  sink()
  x <- tidy(glm_fit)
  x <- cbind(x,exp_estimate=exp(x$estimate),exp(confint.default(glm_fit)))
  x$HR <- paste0(signif(x$exp_estimate,2)," (",signif(x$`2.5 %`,2),"-",signif(x$`97.5 %`,2),")")
  write.csv(x,paste0(dir,"coeffs_",gsub("_regression_output|\\.RDS","",fnms[i]),".csv"),row.names = F)
  y[[i]] <- x
}

res <- do.call(cbind,y)
write.csv(res,paste0(dir,"coeffs_death_all_2.csv"),row.names = F)
res1 <- cbind(res$term,res[,names(res)=="HR"])
write.csv(res1,paste0(dir,"HRs_2.csv"),row.names=F)
