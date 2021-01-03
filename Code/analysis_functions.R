calc_exposure_time <- function(x){
  j <- (x$first_report_date==max(x$first_report_date))
  susc <- x$population[j]-x$cum_cases[j]
  sus_exp_time <- susc*x$time[j] # exposure time of individuals still susceptible
  inf_exp_time <- sum(x$n * x$time) # exposure time of individuals that have been infected
  exp_time <- (sus_exp_time + inf_exp_time)/1e5 # total exposure time in units of 100,000 person-days
  return(list(susc=susc,exp_time=exp_time))
}

calc_survival_time <- function(x){
  j <- (x$date_of_death==max(x$date_of_death))
  alive <- x$population[j] - x$cum_deaths[j]
  alive_surv_time <- alive * x$time_death[j]
  dead_surv_time <- sum(x$n_deaths * x$time_death)
  surv_time <- (alive_surv_time + dead_surv_time)/1e6 # total survival time in units of 1,000,000 person-days
  return(list(alive=alive,surv_time=surv_time))
}

run_regression <- function(dir,fnm){
  # Load data
  y <- read.csv(paste0(dir,fnm),stringsAsFactors = F)
  y$first_report_date <- as.Date(y$first_report_date)
  
  # Calculate total exposure time for each risk group
  inc <- y[y$first_report_date==max(y$first_report_date),]
  inc$susc <- NA
  inc$exp_time <- NA
  for (i in 1:nrow(inc)){
    j <- (y$county_res==inc$county_res[i] & y$age_cat==inc$age_cat[i] & y$sex==inc$sex[i] & y$race_ethnicity==inc$race_ethnicity[i])
    res <- calc_exposure_time(y[j,])
    inc$susc[i] <- res$susc
    inc$exp_time[i] <- res$exp_time
  }
  # write.csv(inc,paste0(dir,gsub("processed_data","inc",fnm)),row.names = F)
  
  # Run Poisson regression
  glm_fit <- glm(cum_cases ~ offset(log(exp_time)) + county_res + age_cat + sex + race_ethnicity,family = poisson(link=log),data = inc)
  # Print summary
  print(summary(glm_fit))
  # Print hazard ratios
  print(exp(coef(glm_fit)))
  
  return(glm_fit)
}

run_death_regression <- function(dir,fnm){
  # Load data
  y <- read.csv(paste0(dir,fnm),stringsAsFactors = F)
  y$date_of_death <- as.Date(y$date_of_death)
  
  # Calculate total survival time for each risk group
  inc <- y[y$date_of_death==max(y$date_of_death),]
  inc$alive <- NA
  inc$surv_time <- NA
  for (i in 1:nrow(inc)){
    j <- (y$county_res==inc$county_res[i] & y$age_cat==inc$age_cat[i] & y$sex==inc$sex[i] & y$race_ethnicity==inc$race_ethnicity[i])
    res <- calc_survival_time(y[j,])
    inc$alive[i] <- res$alive
    inc$surv_time[i] <- res$surv_time
  }
  # write.csv(inc,paste0(dir,gsub("processed_data","inc",fnm)),row.names = F)
  
  # Run Poisson regression
  glm_fit <- glm(cum_deaths ~ offset(log(surv_time)) + county_res + age_cat + sex + race_ethnicity,family = poisson(link=log),data = inc)
  # Print summary
  print(summary(glm_fit))
  # Print hazard ratios
  print(exp(coef(glm_fit)))
  
  return(glm_fit)
}

estimate_death_risk <- function(dir,fnm){
  # Load regression output 
  load(paste0(dir,fnm))
  
  # Load aggregated incidence data
  x <- glm_fit$data
  # Extract fitted cases from regression output
  x$fitted_cases <- NA
  x$fitted_cases[!is.na(x$population)] <- fitted(glm_fit)
  
  # Fit Poisson model for death risk by age and sex
  glm_death_fit <- glm(cum_deaths ~ offset(log(fitted_cases)) + age_cat + sex,family = poisson(link=log),data = x)
  return(glm_death_fit)
}

calc_death_risk_probs <- function(glm_death_fit){
  x <- expand.grid(age_cat=unique(glm_death_fit$data$age_cat),sex=c(0,1))
  x$fitted_cases <- 100
  
  log_mu <- predict(glm_death_fit,newdata = x)
  x$lambda <- exp(log_mu - log(x$fitted_cases))
  
  res <- dcast(x,age_cat~sex,value.var = "lambda")
  row.names(res) <- res[,"age_cat"]
  res$age_cat <- NULL
  names(res) <- c("female","male")
  
  return(res)
}

process_regression_output <- function(dir,fnms){
  for (i in 1:length(fnms)){
    glm_fit <- readRDS(paste0(dir,fnms[i]))
    sink(file = paste0(dir,sub("\\.RDS","\\.txt",fnms[i])))
    print(summary(glm_fit))
    sink()
    x <- tidy(glm_fit)
    x$exp_estimate <- exp(x$estimate)
    write.csv(x,paste0(dir,"coeffs",gsub("regression_output|\\.RDS","",fnms[i]),".csv"),row.names = F)
  }
}