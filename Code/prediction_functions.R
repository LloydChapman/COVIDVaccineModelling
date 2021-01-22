calc_prop_v <- function(prop,v_e){
  prop_v <- prop * (1-v_e)
  return(prop_v)
}

calc_num_averted <- function(diff_prop,susc,idx){
  num_averted <- rep(0,length(diff_prop))
  num_averted[idx] <- diff_prop[idx] * susc[idx]
  return(num_averted)
}

plot_fit <- function(dir,fnm,deaths_ratio){
  load(paste0(dir,fnm))
  datestr <- gsub("regression_output|\\.RData","",fnm)
  
  x <- glm_fit$data
  x$fitted_cases <- NA
  x$fitted_cases[!is.na(x$population)] <- fitted(glm_fit)
  x <- x[!is.na(x$population),]
  x$sex <- ifelse(x$sex==0,"female","male")
  
  # Plot observed vs fitted cases
  dir1 <- "../Figures/FittedVsObservedCases/"
  fnm1 <- "fitted_vs_obsvd_cases_by_"
  plot_grouped_barchart(x,c("cum_cases","fitted_cases"),"county_res","County","cases",c("observed","fitted"),dir1,fnm1,datestr)
  plot_grouped_barchart(x,c("cum_cases","fitted_cases"),"age_cat","Age","cases",c("observed","fitted"),dir1,fnm1,datestr)
  plot_grouped_barchart(x,c("cum_cases","fitted_cases"),"sex","Sex","cases",c("observed","fitted"),dir1,fnm1,datestr)
  plot_grouped_barchart(x,c("cum_cases","fitted_cases"),"race_ethnicity","Race/ethnicity","cases",c("observed","fitted"),dir1,fnm1,datestr)
  
  # Predicted deaths from decision tree model with literature estimates
  p_death <- calc_hosp_ICU_death_risk()$p_death_given_clin_by_age_and_sex * deaths_ratio  
  x$pred_deaths <- NA
  for (i in 1:nrow(x)){
    x$pred_deaths[i] <- x$fitted_cases[i]*p_death[rownames(p_death)==x$age_cat[i],colnames(p_death)==x$sex[i]]
  }
  # Plot observed vs fitted deaths
  dir2 <- "../Figures/PredictedVsObservedDeathsDecisionTree/"
  fnm2 <- "predicted_vs_obsved_deaths_by_"
  plot_grouped_barchart(x,c("cum_deaths","pred_deaths"),"county_res","County","deaths",c("observed","predicted"),dir2,fnm2,datestr)
  plot_grouped_barchart(x,c("cum_deaths","pred_deaths"),"age_cat","Age","deaths",c("observed","predicted"),dir2,fnm2,datestr)
  plot_grouped_barchart(x,c("cum_deaths","pred_deaths"),"sex","Sex","deaths",c("observed","predicted"),dir2,fnm2,datestr)
  plot_grouped_barchart(x,c("cum_deaths","pred_deaths"),"race_ethnicity","Race/ethnicity","deaths",c("observed","predicted"),dir2,fnm2,datestr)
  
  # Predicted deaths from regression of CDPH deaths data
  p_death <- readRDS(paste0(dir,"prob_death_given_clin_by_age_and_sex",datestr,".RDS"))  
  x$pred_deaths <- NA
  for (i in 1:nrow(x)){
    x$pred_deaths[i] <- x$fitted_cases[i]*p_death[rownames(p_death)==x$age_cat[i],colnames(p_death)==x$sex[i]]
  }
  
  # Plot observed vs predicted deaths
  dir3 <- "../Figures/PredictedVsObservedDeaths/"
  fnm3 <- "predicted_vs_obsved_deaths_by_"
  plot_grouped_barchart(x,c("cum_deaths","pred_deaths"),"county_res","County","deaths",c("observed","predicted"),dir3,fnm3,datestr)
  plot_grouped_barchart(x,c("cum_deaths","pred_deaths"),"age_cat","Age","deaths",c("observed","predicted"),dir3,fnm3,datestr)
  plot_grouped_barchart(x,c("cum_deaths","pred_deaths"),"sex","Sex","deaths",c("observed","predicted"),dir3,fnm3,datestr)
  plot_grouped_barchart(x,c("cum_deaths","pred_deaths"),"race_ethnicity","Race/ethnicity","deaths",c("observed","predicted"),dir3,fnm3,datestr)
  
}

plot_grouped_barchart <- function(x,vrble,grp,xlbl,ylbl,lbls,dir,fnm,datestr=NULL){
  agg_x <- aggregate_by_grp(x,vrble,grp)
  agg_x_long <- melt(agg_x,id.vars = grp)
  p <- ggplot(agg_x_long,aes(x=agg_x_long[,grp],y=value,group=variable,fill=variable)) + geom_bar(position = "dodge",stat = "identity") + xlab(xlbl) + ylab(ylbl) + scale_fill_discrete(name = "",labels = lbls)
  if (grp == "county_res"){
    p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1)) 
  }
  dir.create(dir,recursive = T)
  if (is.null(datestr)){
    pdf(paste0(dir,fnm,grp,".pdf"),width = 7,height = 4)  
  } else {
    pdf(paste0(dir,fnm,grp,datestr,".pdf"),width = 7,height = 4)
  }
  print(p)
  dev.off()
}

predict_cases_hosps_deaths <- function(fnm,t_sim,r,cases_ratio,deaths_ratio,p_death){
  glm_fit <- readRDS(fnm)
  
  # Make data frame for predicting future cases
  x <- glm_fit$data
  x$exp_time <- x$susc * t_sim/1e5
  
  # # Predict cumulative cases for next t_sim days
  # x$cum_cases_pred <- predict(glm_fit,newdata = x,type = "response")
  # print(sum(x$cum_cases_pred,na.rm = T))
  # x$prop_case <- x$cum_cases_pred/x$susc 
  
  # # Sort 
  # x <- x[order(x$prop_case,decreasing = T),]
  # View(x)
  
  # Predict mean number of cases in each demographic subgroup
  log_mu <- predict(glm_fit,newdata = x)
  # View(log_mu)
  x$lambda <- exp(log_mu - log(x$exp_time))/1e5
  
  # Predict proportion of cases accounting for growth in force of infection
  x$prop_case <- 1-exp(-x$lambda*r*(1-r^t_sim)/(1-r) * cases_ratio) 
  x$cum_cases_pred <- x$prop_case * x$susc
  
  x$sex <- ifelse(x$sex==0,"female","male")
  res <- calc_hosp_ICU_death_risk()
  p_hosp <- res$p_hosp_given_clin
  # p_death <- res$p_death_given_clin_by_age_and_sex
  x$prop_hosp <- NA
  x$prop_death <- NA
  for (i in 1:nrow(x)){
    x$prop_hosp[i] <- x$prop_case[i]*p_hosp[rownames(p_hosp)==x$age_cat[i],colnames(p_hosp)==x$sex[i]] * deaths_ratio # N.B. Correction may be wrong here
    x$prop_death[i] <- x$prop_case[i]*p_death[rownames(p_death)==x$age_cat[i],colnames(p_death)==x$sex[i]] * deaths_ratio
  }
  
  x$cum_hosps_pred <- x$prop_hosp * x$susc
  x$cum_deaths_pred <- x$prop_death * x$susc
  
  return(x)
}

predict_deaths <- function(glm_fit,t_sim,r,deaths_ratio,IFR_long){
  # Make data frame for predicting future cases
  x <- glm_fit$data
  
  # Predict mean number of cases in each demographic subgroup
  log_mu <- predict(glm_fit,newdata = x)
  x$lambda <- exp(log_mu)/(x$surv_time*1e6)
  
  # Predict proportion of cases accounting for growth in force of infection
  if (r==1){
    x$prop_death <- 1-exp(-x$lambda*t_sim * deaths_ratio)
  } else {
    x$prop_death <- 1-exp(-x$lambda*r*(1-r^t_sim)/(1-r) * deaths_ratio)
  }
  # x1 <- x
  # x1$n_deaths <- NULL
  # names(x1)[names(x1)=="cum_deaths"] <- "n_deaths"
  # x1 <- backcalculate_infections(x1,IFR_long)
  # x <- merge(x,x1[,names(x1)!="n_deaths"],by=c("county_res","age_cat","sex","race_ethnicity"))
  # x$susc <- pmax(0,x$population - x$n)
  # x$cum_deaths_pred <- x$prop_death * x$susc
  x$cum_deaths_pred <- x$prop_death * x$alive
  
  return(x)
}

compare_predictions <- function(dir,fnm,x2,fdir,deaths_only=F){
  x <- read.csv(paste0(dir,fnm),stringsAsFactors = F)
  if (!deaths_only){
    agg_x <- aggregate(cbind(cum_cases_pred,cum_hosps_pred,cum_deaths_pred) ~ county_res,x,sum)
  } else {
    agg_x <- aggregate(cbind(cum_deaths_pred) ~ county_res,x,sum)
  }
  agg_x <- merge(agg_x,x2)
  dir.create(fdir,recursive = T)
  if (!deaths_only){
    cum_cases <- melt(agg_x[,c("county_res","cum_cases_pred","point")],id.vars = "county_res",value.name = "cases")
    pdf(paste0(fdir,"compare_cases_",sub(".csv","",fnm),"_county.pdf"),width = 7,height = 4)
    print(ggplot(cum_cases,aes(fill=variable,y=cases,x=county_res)) + geom_bar(position = "dodge", stat="identity") + xlab("County") + ylab("Predicted cases in next 4 weeks") + scale_fill_discrete(name = "",labels = c("Model prediction","Forecast Hub")) + theme(axis.text.x = element_text(angle = 45,hjust = 1)))
    dev.off()
    cum_cases_state <- data.frame(model=factor(c("Model prediction","Forecast Hub"),level=c("Model prediction","Forecast Hub")),cases=apply(agg_x[,c("cum_cases_pred","point")],2,sum))
    cases_ratio <- cum_cases_state$cases[2]/cum_cases_state$cases[1]
    pdf(paste0(fdir,"compare_cases_",sub(".csv","",fnm),"_state.pdf"),width = 4,height = 4)
    print(ggplot(cum_cases_state,aes(x=model,y=cases)) + geom_bar(stat = "identity") + theme(axis.title.x = element_blank()))
    dev.off()
  }
  cum_deaths <- melt(agg_x[,c("county_res","cum_deaths_pred","deaths")],id.vars = "county_res",value.name = "deaths")
  pdf(paste0(fdir,"compare_deaths_",sub(".csv","",fnm),"_county.pdf"),width = 7,height = 4)
  print(ggplot(cum_deaths,aes(fill=variable,y=deaths,x=county_res)) + geom_bar(position = "dodge", stat="identity") + xlab("County") + ylab("Predicted deaths in next 4 weeks") + scale_fill_discrete(name = "",labels = c("Model prediction","CalCAT")) + theme(axis.text.x = element_text(angle = 45,hjust = 1)))
  dev.off()
  cum_deaths_state <- data.frame(model=names(agg_x[,c("cum_deaths_pred","deaths")]),deaths=apply(agg_x[,c("cum_deaths_pred","deaths")],2,sum))
  cum_deaths_state <- data.frame(model=factor(c("Model prediction","CalCAT"),level=c("Model prediction","CalCAT")),deaths=apply(agg_x[,c("cum_deaths_pred","deaths")],2,sum))
  deaths_ratio <- cum_deaths_state$deaths[2]/cum_deaths_state$deaths[1]
  pdf(paste0(fdir,"compare_deaths_",sub(".csv","",fnm),"_state.pdf"),width = 4,height = 4)
  print(ggplot(cum_deaths_state,aes(x=model,y=deaths)) + geom_bar(stat = "identity") + theme(axis.title.x = element_blank()))
  dev.off()
  if (!deaths_only){
    return(list(cases_ratio=cases_ratio,deaths_ratio=deaths_ratio))  
  } else {
    return(list(deaths_ratio=deaths_ratio))
  }
}

add_risk_ests <- function(fnm,RR,df,p_death,lt_long,seroprev){
  glm_fit <- readRDS(fnm)
  
  # Make data frame for predicting future cases
  x <- glm_fit$data
  
  # Get infection rate based on demographic factors
  log_mu <- predict(glm_fit,newdata = x)
  x$lambda <- exp(log_mu)/(x$exp_time*1e5)
  
  # Calculate RR for non-special population
  # x$non_spec_pop <- x$population - x$hcw - x$prisoner - x$snf - x$teacher - x$homeless
  # x$RR0 <- x$population/(x$non_spec_pop + RR["HCW"]*x$hcw + RR["prisoner"]*x$prisoner + RR["SNF"]*x$snf + RR["teacher"]*x$teacher + RR["homeless"]*x$homeless)
  x$RR0 <- x$population/(x$non_spec_pop + RR["HCW"]*x$hcw + RR["prisoner"]*x$prisoner + RR["SNF"]*x$snf + RR["teacher"]*x$teacher + RR["homeless"]*x$homeless + RR["essential"]*x$essential)
  
  # Adjust relative risks for special populations based on RR for non-special population
  x$RR1 <- x$RR0 * RR["HCW"]
  x$RR2 <- x$RR0 * RR["prisoner"]
  x$RR3 <- x$RR0 * RR["SNF"]
  x$RR4 <- x$RR0 * RR["teacher"]
  x$RR5 <- x$RR0 * RR["homeless"]
  x$RR6 <- x$RR0 * RR["essential"]
  
  # Calculate past cases in each special population
  # x$cum_cases0 <- round(x$RR0*x$non_spec_pop/x$population * x$cum_cases)
  x$cases1 <- round(x$RR1 * x$hcw/x$population * x$cum_cases)
  x$cases2 <- round(x$RR2 * x$prisoner/x$population * x$cum_cases)
  x$cases3 <- round(x$RR3 * x$snf/x$population * x$cum_cases)
  x$cases4 <- round(x$RR4 * x$teacher/x$population * x$cum_cases)
  x$cases5 <- round(x$RR5 * x$homeless/x$population * x$cum_cases)
  x$cases6 <- round(x$RR6 * x$essential/x$population * x$cum_cases)
  # Assign remainder of cases to non-special population as it is generally the largest subgroup so rounding makes least difference
  # x$cases0 <- x$cum_cases - apply(x[,paste0("cases",1:5)],1,sum)
  x$cases0 <- x$cum_cases - apply(x[,paste0("cases",1:6)],1,sum)
  
  # Assign numbers to each risk group
  x$grp <- 1:nrow(x)
  
  # x_long <- melt(x,id.vars = c("grp","county_res","age_cat","sex","race_ethnicity","lambda"),measure.vars = c(paste0("cum_cases",0:5),paste0("RR",0:5)),variable.name = "special.population",value.name = "RR")
  # x_long$special.population <- as.numeric(sub("RR","",x_long$special.population))
  # x_long <- reshape(x[,c("county_res","age_cat","sex","race_ethnicity","time","lambda",paste0("RR",0:5),paste0("cases",0:5),"grp")],varying = list(paste0("RR",0:5),paste0("cases",0:5)),v.names = c("RR","cases"),timevar = "special.population",idvar = "grp",direction = "long")
  x_long <- reshape(x[,c("county_res","age_cat","sex","race_ethnicity","time","lambda",paste0("RR",0:6),paste0("cases",0:6),"grp")],varying = list(paste0("RR",0:6),paste0("cases",0:6)),v.names = c("RR","cases"),timevar = "special.population",idvar = "grp",direction = "long")
  x_long$special.population <- x_long$special.population - 1
  x_long$risk_grp <- 1:nrow(x_long)
  
  # Merge with synthetic population data frame keeping only rows in synthetic population data frame as San Luis Obispo County has no 80+ non-Hispanic black females (so introduces NA's)
  df <- merge(df,x_long,by=c("county_res","age_cat","sex","race_ethnicity","special.population"),all.x=T)
  
  # Merge with life table data frame
  df$sex <- ifelse(df$sex==1,"male","female")
  df <- merge(df,lt_long,by=c("age","sex","race_ethnicity"),all.x=T)
  
  # Merge with seroprevalence data frame
  age_cat_sero_lbs <- c(as.numeric(gsub("-.*|\\+|>","",seroprev$age_cat_sero)),max(df$age)+1)
  df$age_cat_sero <- cut(df$age,age_cat_sero_lbs,right = F,labels = seroprev$age_cat_sero)
  df <- merge(df,seroprev,by="age_cat_sero",all.x=T)
  
  # Get age-, sex- and comorbidity-dependent probabilities of hospitalisation and death
  # res <- calc_hosp_ICU_death_risk()
  # # p_hosp <- res$p_hosp_given_clin
  # # p_death_given_hosp <- res$p_death_given_hosp_by_age_and_sex
  # p_hosp <- as.data.frame(res$p_hosp_given_clin)
  # p_hosp$age_cat <- row.names(p_hosp)
  # p_hosp_long <- melt(p_hosp,id.vars = "age_cat",variable.name = "sex",value.name = "p_hosp")
  # p_death_given_hosp <- as.data.frame(res$p_death_given_hosp_by_age_and_sex)
  # p_death_given_hosp$age_cat <- row.names(p_death_given_hosp)
  # p_death_given_hosp_long <- melt(p_death_given_hosp,id.vars = "age_cat",variable.name = "sex",value.name = "p_death_given_hosp")
  
  # p_death <- as.data.frame(res$p_death_given_clin_by_age_and_sex)
  p_death$age_cat <- row.names(p_death)
  p_death_long <- melt(p_death,id.vars = "age_cat",variable.name = "sex",value.name = "p_death")
  
  # Merge with synthetic population data frame
  # df <- merge(df,p_hosp_long,by=c("age_cat","sex"),all.x=T)
  # df <- merge(df,p_death_given_hosp_long,by=c("age_cat","sex"),all.x=T)
  df <- merge(df,p_death_long,by=c("age_cat","sex"),all.x=T)
  
  return(df)
  # tstart <- Sys.time()
  # df$p_hosp <- NA
  # df$p_death <- NA
  # for (i in 1:nrow(df)){
  #   df$p_hosp[i] <- p_hosp[rownames(p_hosp)==df$age_cat[i],colnames(p_hosp)==df$sex[i]] * deaths_ratio # N.B. Correction may be wrong here
  #   df$p_death[i] <- p_death[rownames(p_death)==df$age_cat[i],colnames(p_death)==df$sex[i]] * deaths_ratio
  # }
  # tend <- Sys.time()
  # print(tend-tstart)
}

add_risk_ests_deaths <- function(glm_fit,alpha,df,RR,RR_LB,RR_UB,HR,lt_long,seroprev,IFR_long,IFR_ratio){
  # Add 1 for RR for non-special population to start of RR vector
  RR <- c(non_spec_pop=1,RR)
  RR_LB <- c(non_spec_pop=1,RR_LB)
  RR_UB <- c(non_spec_pop=1,RR_UB)
  
  # # Read in death regression output
  # glm_fit <- readRDS(fnm)
  
  # Make data frame for predicting future cases
  x <- as.data.table(glm_fit$data)
  
  # Get death rate and its 95% CI based on demographic factors
  z_star <- qnorm(1-alpha/2)
  preds <- predict(glm_fit,newdata=x,se.fit=T)
  log_mu <- preds$fit
  se_log_mu <- preds$se.fit
  log_mu_LB <- preds$fit - z_star * se_log_mu
  log_mu_UB <- preds$fit + z_star * se_log_mu
  x[,`:=`(log_lambda=log_mu-log(x[,surv_time]*1e6),se_log_lambda=se_log_mu,log_lambda_LB=log_mu_LB-log(x[,surv_time]*1e6),log_lambda_UB=log_mu_UB-log(x[,surv_time]*1e6))]
  
  # Merge into synthetic population data table
  demogrphcs <- c("county_res","age_cat","sex","race_ethnicity")
  df <- merge(df,x[,c(demogrphcs,"time_death","log_lambda","se_log_lambda","log_lambda_LB","log_lambda_UB"),with=F],by=demogrphcs)
  
  # Calculate RR for non-special population
  comorbs <- names(HR)
  spec_pop_comorbs <- c("special.population",comorbs)
  # Make data table of all unique combinations of special populations and comorbidities
  y <- CJ(special.population=unique(df[,special.population]),asthma=c(0,1),diabetes=c(0,1),smoker=c(0,1),heart.disease=c(0,1),heart.failure=c(0,1),hypertension=c(0,1),obesity=c(0,1))
  n_y <- nrow(y)
  
  # Make matrix of HRs for different comorbidities
  HR_mat <- matrix(rep(HR,each=n_y),nrow=n_y)
  # Multiply HRs together for each comorbidity combination
  y[,HR := apply(HR_mat * y[,comorbs,with=F],1,function(x){prod(x[x!=0])})] # Works as prod(numeric(0))=1, so gives right answer when all comorbidity statuses are 0
  # Multiply RR by overall HR for comorbidity combination
  y[,`:=`(RR = RR[special.population+1] * HR,RR_LB = RR_LB[special.population+1] * HR,RR_UB = RR_UB[special.population+1] * HR)]
  # Calculate standard error for relative risk
  y[,se_RR := pmin(RR-RR_LB,RR_UB-RR)/z_star]
  
  # Aggregate synthetic population data table by demographics, special populations, and comorbidities
  agg_df <- df[,.(population=.N),by=c(demogrphcs,spec_pop_comorbs)]
  
  # Make data table of all unique combinations of demographic factors, special populations, and comorbidities
  z <- CJ(county_res=unique(agg_df[,county_res]),age_cat=unique(agg_df[,age_cat]),sex=unique(agg_df[,sex]),race_ethnicity=unique(agg_df[,race_ethnicity]),special.population=unique(agg_df[,special.population]),asthma=c(0,1),diabetes=c(0,1),smoker=c(0,1),heart.disease=c(0,1),heart.failure=c(0,1),hypertension=c(0,1),obesity=c(0,1))
  
  # Merge with aggregated synthetic population data table
  agg_df <- merge(z,agg_df,all.x=T)
  agg_df[is.na(population),population:=0]
  # Add group index for special population and comorbidity combinations
  agg_df[,grp:=.GRP,by=c(demogrphcs,spec_pop_comorbs)]
  # Merge with RR data table
  agg_df <- merge(agg_df,y,by=spec_pop_comorbs,all.x=T)
  # Adjust RRs so that overall death rate remains the same
  agg_df <- merge(agg_df,agg_df[,.(RR0 = sum(population)/sum(population*RR)),by=demogrphcs],by=demogrphcs)
  agg_df <- agg_df[,`:=`(RR = RR*RR0,RR_LB = RR_LB*RR0,RR_UB = RR_UB*RR0,se_RR = se_RR*RR0)]
  
  # Merge with synthetic population data table by demographic, special population and comorbidity combination
  df <- merge(df,agg_df,by=c(demogrphcs,spec_pop_comorbs))
  
  # Add risk group variable for demographic and special population combination
  df[,risk_grp:=.GRP,by=c(demogrphcs,"special.population")]
  
  # Merge with life table data frame
  df[,sex:=ifelse(sex==1,"male","female")]
  df <- merge(df,lt_long,by=c("age","sex","race_ethnicity"),all.x=T)
  
  # Merge with seroprevalence data frame
  age_cat_sero_lbs <- c(as.numeric(gsub("-.*|\\+|>","",seroprev[,age_cat_sero])),max(df[,age])+1)
  df[,age_cat_sero:=cut(age,age_cat_sero_lbs,right = F,labels = seroprev[,age_cat_sero])]
  df <- merge(df,seroprev,by="age_cat_sero",all.x=T)
  
  # Merge with probability of past infection based on estimated infections from deaths and IFR
  x[,sex:=ifelse(sex==1,"male","female")]
  x <- merge(x,IFR_long[,.(age_cat,sex,median_perc,CI_95_LB,CI_95_UB)],by=c("age_cat","sex"))
  x[,IFR:=median_perc/100*IFR_ratio]
  x[,cum_inf:=cum_deaths/IFR]
  agg_x <- x[,.(IFR=mean(IFR),cum_inf=sum(cum_inf),population=sum(population,na.rm=T)),by=.(age_cat,sex)]
  agg_x[,p_past_inf:=cum_inf/population]
  df <- merge(df,agg_x[,.(age_cat,sex,IFR,p_past_inf)],by=c("age_cat","sex"))
  return(df)
}

backcalculate_infections <- function(x,IFR_long){
  agg_x <- aggregate(n_deaths ~ county_res + age_cat + sex + race_ethnicity,x,sum)
  # Merge IFRs by age and sex
  agg_x <- merge(agg_x,IFR_long[,c("age_cat","sex","median_perc","CI_95_LB","CI_95_UB")],by=c("age_cat","sex"),all.x = T)
  # Calculate estimated total number of infections using deaths and O'Driscoll IFR
  agg_x$n <- agg_x$n_deaths/(agg_x$median_perc/100)
  return(agg_x)
}

calc_cases_from_infections <- function(x,p_clin){
  x <- merge(x,p_clin,by="age_cat",all.x=T)
  x$n_cases <- x$n * x$p_clin
  return(x)
}

# Version for data.table input x
backcalculate_infections2 <- function(x,IFR_long,IFR_ratio){
  x <- merge(x,IFR_long[,.(age_cat,sex,median_perc,CI_95_LB,CI_95_UB)],by=c("age_cat","sex"),all.x = T)
  x[,`:=`(infection=death/(median_perc/100*IFR_ratio),infection_v=death_v/(median_perc/100*IFR_ratio))]
}

calc_cases_from_infections2 <- function(x,p_clin_dt){
  x <- merge(x,p_clin_dt,by="age_cat",all.x=T)
  x[,`:=`(case=infection*p_clin,case_v=infection_v*p_clin)]
}

calc_DALYs <- function(death,LE,case,p_mi,d,w){
  YLL <- death * LE
  YLD <- ((case & !death) * (p_mi * d[1] * w[1] + (1-p_mi) * d[2] * w[2]) + death * d[3] * w[3])/365
  # return(list(YLL=YLL,YLD=YLD))
  DALYs <- YLL + YLD
  return(DALYs)
}
  
calc_DALYs_from_deaths <- function(death,LE,d,w){
  YLL <- death * LE
  YLD <- death * d * w/365
  DALYs <- YLL + YLD
  return(DALYs)
}

calc_DALYs_from_infections_and_cases <- function(infection,case,d,w){
  DALYs <- ((infection - case) * d[1] * w[1] + case * d[2] * w[2])/365
  return(DALYs)
}

simulate_cases_hosps_deaths <- function(df,v_e,t_sim,r,cases_ratio,deaths_ratio,n_sim,p_mi,d,w){
  # # Simulate which individuals are already infected
  # risk_grps <- unique(df$risk_grp)
  # risk_grp <- df$risk_grp
  # 
  # tstart <- Sys.time()
  # past_inf <- rep(F,nrow(df))
  # for (i in risk_grps){
  #   j <- (risk_grp==i)
  #   n <- df$cases[j][1]
  #   # past_inf[j][sample.int(sum(j),n,useHash = T)] <- 1
  #   past_inf[j][dqsample.int(sum(j),n)] <- T
  # }
  # df$past_inf <- past_inf
  # tend <- Sys.time()
  # print(tend-tstart)
  # 
  # tstart1 <- Sys.time()
  # ldf <- split(df,grp)
  # for (i in 1:length(ldf)){
  #   n <- nrow(ldf[[i]])
  #   past_inf <- rep(F,n)
  #   past_inf[dqsample.int(n,ldf[[i]]$cases[1])] <- T
  #   ldf[[i]]$past_inf <- past_inf
  # }
  # df1 <- do.call(rbind,ldf)
  # tend1 <- Sys.time()
  # print(tend1-tstart1)
  
  # Total population
  N <- nrow(df)
  
  # ## No vaccination
  # # Simulate past cases
  # df$p_past_case <- 1 - exp(-df$lambda*df$time[1] * df$RR)
  # df$past_case <- (df$p_past_case > runif(N))
  # 
  # # Simulate new cases
  # df$p_case <- 1 - exp(-df$lambda*r*(1-r^t_sim)/(1-r) * df$RR * cases_ratio)
  # df$case <- (!df$past_case & df$p_case > runif(N))
  # 
  # # Simulate hospitalisations and deaths
  # df$hosp <- (df$case & df$p_hosp > runif(N))
  # df$death <- (df$hosp & df$p_death > runif(N))
  # 
  # ## All-or-nothing vaccination
  # # Simulate new cases
  # p_case_v <- (1 - exp(-(1-v_e)*df$lambda*r*(1-r^t_sim)/(1-r) * df$RR * cases_ratio))
  # df$case_v <- (!df$past_case & df$p_case_v > runif(N))
  # 
  # # Simulate hospitalisations and deaths
  # df$hosp_v <- (df$case_v & df$p_hosp > runif(N))
  # df$death_v <- (df$hosp_v & df$p_death > runif(N))
  # 
  # return(df)
  
  # p_past_case <- 1 - exp(-df$lambda*df$time[1] * df$RR)
  p_past_case <- df$seroprev * df$RR
  p_case <- 1 - exp(-df$lambda*r*(1-r^t_sim)/(1-r) * df$RR * cases_ratio)
  p_case_v <- (1 - exp(-(1-v_e)*df$lambda*r*(1-r^t_sim)/(1-r) * df$RR * cases_ratio))

  # out_nms <- c("past_case","case","hosp","death","case_v","hosp_v","death_v")
  # out_nms <- c("past_case","case","death","YLL","YLD","case_v","death_v","YLL_v","YLD_v")
  out_nms <- c("past_case","case","death","DALYs","case_v","death_v","DALYs_v")
  out <- array(F,dim = c(nrow(df),length(out_nms),n_sim))
  dimnames(out)[[1]] <- row.names(df)
  dimnames(out)[[2]] <- out_nms
  for (i in 1:n_sim){
    ## No vaccination
    # Simulate past cases
    past_case <- (p_past_case > runif(N))

    # Simulate new cases
    case <- (!past_case & p_case > runif(N))

    # # Simulate hospitalisations and deaths
    # hosp <- (case & (df$p_hosp * deaths_ratio) > runif(N))
    # death <- (hosp & df$p_death_given_hosp > runif(N))
    # Simulate deaths
    death <- (case & (df$p_death * deaths_ratio) > runif(N))
    # death <- (case & df$p_death > runif(N))
    # Calculate years of life lost (YLL) and years lived with disability (YLD)
    # res <- calc_DALYs(death,df$life_expectancy,case,p_mi,d,w)
    # YLL <- res$YLL
    # YLD <- res$YLD
    DALYs <- calc_DALYs(death,df$life_expectancy,case,p_mi,d,w)

    ## Vaccination
    # Simulate new cases
    case_v <- (!past_case & p_case_v > runif(N))

    # # Simulate hospitalisations and deaths
    # hosp_v <- (case_v & (df$p_hosp * deaths_ratio) > runif(N))
    # death_v <- (hosp_v & df$p_death_given_hosp > runif(N))
    # Simulate deaths
    death_v <- (case_v & (df$p_death * deaths_ratio) > runif(N))
    # death_v <- (case_v & df$p_death > runif(N))
    # Calculate years of life lost (YLL) and years lived with disability (YLD)
    # res_v <- calc_DALYs(death_v,df$life_expectancy,case_v,p_mi,d,w)
    # YLL_v <- res_v$YLL
    # YLD_v <- res_v$YLD
    DALYs_v <- calc_DALYs(death_v,df$life_expectancy,case_v,p_mi,d,w)
    
    # out[,,i] <- cbind(past_case,case,hosp,death,case_v,hosp_v,death_v)
    # out[,,i] <- cbind(past_case,case,death,YLL,YLD,case_v,death_v,YLL_v,YLD_v)
    out[,,i] <- cbind(past_case,case,death,DALYs,case_v,death_v,DALYs_v)
  }

  return(out)
}

# simulate_deaths <- function(df,lambda,RR,v_e,t_sim,r,deaths_ratio,n_sim,d,w){
simulate_deaths <- function(df,v_e,t_sim,r,deaths_ratio,n_sim,d,w){
  # Total population
  N <- nrow(df)
  
  p_past_inf <- df$p_past_inf
  log_lambda <- rtruncnorm(N,a = df$log_lambda_LB,b = df$log_lambda_UB,mean = df$log_lambda,sd = df$se_log_lambda)
  lambda <- exp(log_lambda)
  RR <- rtruncnorm(N,a = df$RR_LB - 1e-15,b = df$RR_UB + 1e-15,mean = df$RR,sd = df$se_RR)
  # lambda <- exp(df$log_lambda)
  # RR <- df$RR
  if (r!=1){
    p_death <- 1 - exp(-lambda*r*(1-r^t_sim)/(1-r) * RR * deaths_ratio)
    p_death_v <- (1 - exp(-(1-v_e)*lambda*r*(1-r^t_sim)/(1-r) * RR * deaths_ratio))    
  } else {
    p_death <- 1 - exp(-lambda*t_sim * RR * deaths_ratio)
    p_death_v <- (1 - exp(-(1-v_e)*lambda*t_sim * RR * deaths_ratio))    
  }
  
  # out_nms <- c("sim","past_inf","death","DALYs","death_v","DALYs_v")
  # out <- matrix(nrow=n_sim*nrow(df),ncol=length(out_nms))
  # # rownames(out) <- row.names(df)
  # colnames(out) <- out_nms
  # for (i in 1:n_sim){
    ## No vaccination
    # Simulate past infections
    past_inf <- (p_past_inf > runif(N))
    
    # Simulate new deaths
    death <- (!past_inf & (p_death > runif(N)))
    
    # Calculate DALYs from deaths
    DALYs <- calc_DALYs_from_deaths(death,df$life_expectancy,d,w)
    
    ## Vaccination
    # Simulate new deaths
    death_v <- (!past_inf & (p_death_v > runif(N)))

    # Calculate DALYs from deaths
    DALYs_v <- calc_DALYs_from_deaths(death_v,df$life_expectancy,d,w)
    
    # Bind into output array
    # out[,,i] <- cbind(past_inf,death,DALYs,death_v,DALYs_v)
    # out[(i-1)*N+(1:N),] <- cbind(rep(i,N),past_inf,death,DALYs,death_v,DALYs_v)
    out <- cbind(past_inf,death,DALYs,death_v,DALYs_v)
  # }
  
  return(out)
}

aggregate_predictions <- function(out,df,n_sim){
  # Bind column with risk group ID
  out <- abind(risk_grp=array(rep(df$risk_grp,n_sim),dim=c(nrow(out),1,n_sim)),out,along=2)
  # out_nms <- colnames(out)
  
  # Aggregate outcomes by risk groups
  agg_out <- array(dim=c(length(unique(df$risk_grp)),dim(out)[2],dim(out)[3]))
  for (i in 1:n_sim){
    # agg_out[,,i] <- as.matrix(aggregate(cbind(past_case,case,hosp,death,case_v,hosp_v,death_v) ~ risk_grp,out[,,i],sum))
    agg_out[,,i] <- as.matrix(aggregate(cbind(past_case,case,death,DALYs,case_v,death_v,DALYs_v) ~ risk_grp,out[,,i],sum))
    # agg_out[,,i] <- as.matrix(aggregate(out[,out_nms,i],by=list(out[,"risk_grp",i]),sum))
  }
  return(agg_out)
}

aggregate_predictions_deaths <- function(out,grp,n_sim){
  # Bind column with risk group ID
  # out <- cbind(grp=rep(df[,grp],n_sim),out)
  out <- cbind(grp,out)
  out_dt <- as.data.table(out)
  
  # Aggregate outcomes by risk groups
  # agg_out <- out_dt[,.(past_inf=sum(past_inf),death=sum(death),DALYs=sum(DALYs),death_v=sum(death_v),DALYs_v=sum(DALYs_v)),by=.(sim,grp)]
  agg_out <- out_dt[,.(past_inf=sum(past_inf),death=sum(death),DALYs=sum(DALYs),death_v=sum(death_v),DALYs_v=sum(DALYs_v)),by=.(grp)]
  
  return(agg_out)
}

acomb <- function(...){abind(...,along = 1)}

process_predictions <- function(out,df,n_sim){
  # # Collapse 3D (indvdls x outcomes x sims) array to 2D ((indvdlsxsims) x outcomes) matrix
  # out_long<-cbind(matrix(aperm(out,c(1,3,2)),nrow=prod(dim(out)[c(1,3)])),rep(seq_len(dim(out)[3]),each=dim(out)[1]))
  # # Name columns
  # colnames(out_long) <- c("past_case","case","hosp","death","case_v","hosp_v","death_v","sim")
  # # Convert to data frame
  # out_long <- data.frame(out_long)
  
  # Aggregate outcomes by risk groups
  agg_out <- aggregate_predictions(out,df,n_sim)
    
  # Calculate mean numbers for each outcome
  mean_agg_out <- apply(agg_out,c(1,2),mean)
  # colnames(mean_agg_out) <- out_nms
  colnames(mean_agg_out) <- c("risk_grp",colnames(out))
  mean_agg_out <- data.frame(mean_agg_out)

  # Merge county, age, sex and race/ethnicity information from synthetic population data frame and add risk group populations
  risk_grp_df <- df[match(unique(df$risk_grp),df$risk_grp),c("county_res","age_cat","sex","race_ethnicity","special.population","risk_grp")]
  risk_grp_pop <- as.data.frame(table(df$risk_grp),responseName = "population")
  risk_grp_df <- merge(risk_grp_df,risk_grp_pop,by.x="risk_grp",by.y="Var1")
  mean_agg_out <- merge(risk_grp_df,mean_agg_out,by="risk_grp")

  return(mean_agg_out)
}

# process_predictions <- function(agg_out,risk_grp_df,n_sim,grp){
#   agg_x <- array(dim = c(length(unique(risk_grp_df[,grp])),ncol(agg_out)-1,n_sim))
#   for (i in 1:n_sim){
#     x <- agg_out[,,i] 
#     colnames(x) <- c("risk_grp","past_case","case","hosp","death","case_v","hosp_v","death_v")
#     x <- merge(risk_grp_df,x,by="risk_grp")
#     tmp <- aggregate_by_grp(x,c("past_case","case","hosp","death","case_v","hosp_v","death_v"),grp)
#     agg_x[,,i] <- as.matrix(tmp[,2:ncol(tmp)])
#   }
#   rownames(agg_x) <- tmp[,1]
#   # colnames(agg_x) <- names(tmp[,2:ncol(tmp)])
#   q_agg_x <- apply(agg_x,c(1,2),function(x)(quantile(x,probs = c(0.5,0.25,0.75))))
#   q_agg_x_long <- melt(q_agg_x)
#   q_agg_x_wide <- acast(q_agg_x_long,Var2 ~ Var3 + Var1,value.var = "value")
#   mean_and_q_agg_x <- data.frame(rownames(q_agg_x_wide),cbind(apply(agg_x,c(1,2),mean),q_agg_x_wide))
#   names(mean_and_q_agg_x)[1] <- grp
#   
#   return(mean_and_q_agg_x)
# }


process_predictions_deaths <- function(out,grp,n_sim,risk_grp_dt,IFR_long,IFR_ratio,p_clin_dt,d,w){
  # Aggregate outcomes by risk groups
  agg_out <- aggregate_predictions_deaths(out,grp,n_sim)
  
  # Merge county, age, sex and race/ethnicity information from synthetic population data frame and add risk group populations
  agg_out <- merge(risk_grp_dt,agg_out,by="grp")
  
  # Backcalculate infections in each risk group using O'Driscoll IFR
  agg_out <- backcalculate_infections2(agg_out,IFR_long,IFR_ratio)
  
  # Calculate clinical cases from infections using Davies age-dependent clinical fraction
  agg_out <- calc_cases_from_infections2(agg_out,p_clin_dt)
  
  # Calculate and add DALYs from infections and cases
  agg_out[,`:=`(DALYs = DALYs + calc_DALYs_from_infections_and_cases(infection,case,d,w),DALYs_v = DALYs_v + calc_DALYs_from_infections_and_cases(infection_v,case_v,d,w))]
  
  # Calculate infections, cases, deaths and DALYs averted
  agg_out <- calc_infections_cases_deaths_DALYs_averted(agg_out)
  
  # # Calculate summary statistics across simulations
  # cols1 <- names(agg_out)[!(names(agg_out) %in% c(cols,"population","sim","past_inf","median_perc","CI_95_LB","CI_95_UB","p_clin"))]
  # res_mean <- agg_out[, lapply(.SD, mean), .SDcols = cols1, by = "grp"]
  # res_median <- agg_out[, lapply(.SD, median), .SDcols = cols1, by = "grp"][, setnames(.SD, cols1, paste(cols1, "median", sep = "_"))]
  # res_q_95_LB <- agg_out[, lapply(.SD, q_95_LB), .SDcols = cols1, by = "grp"][, setnames(.SD, cols1, paste(cols1, "q_95_LB", sep = "_"))]
  # res_q_95_UB <- agg_out[, lapply(.SD, q_95_UB), .SDcols = cols1, by = "grp"][, setnames(.SD, cols1, paste(cols1, "q_95_UB", sep = "_"))]
  # # Bind results together
  # res <- cbind(res_mean,res_median[,!"grp"],res_q_95_LB[,!"grp"],res_q_95_UB[,!"grp"])
  # 
  # # Merge with risk group data table
  # res <- merge(risk_grp_dt,res,by="grp")
  # 
  # return(res)
  return(agg_out)
}

calc_cases_hosps_deaths_averted <- function(x){
  x$prop_case <- x$case/x$population
  x$prop_hosp <- x$hosp/x$population
  x$prop_death <- x$death/x$population
  x$prop_case_v <- x$case_v/x$population
  x$prop_hosp_v <- x$hosp_v/x$population
  x$prop_death_v <- x$death_v/x$population
  
  x$diff_prop_case <- x$prop_case - x$prop_case_v
  x$diff_prop_hosp <- x$prop_hosp - x$prop_hosp_v
  x$diff_prop_death <- x$prop_death - x$prop_death_v
  
  x$cases_averted <- x$case - x$case_v
  x$hosps_averted <- x$hosp - x$hosp_v
  x$deaths_averted <- x$death - x$death_v
  
  return(x)
}

calc_cases_deaths_DALYs_averted <- function(x){
  x$prop_case <- x$case/x$population
  # x$prop_hosp <- x$hosp/x$population
  x$prop_death <- x$death/x$population
  x$DALYs_pp <- x$DALYs/x$population
  x$prop_case_v <- x$case_v/x$population
  # x$prop_hosp_v <- x$hosp_v/x$population
  x$prop_death_v <- x$death_v/x$population
  x$DALYs_pp_v <- x$DALYs_v/x$population
  
  x$diff_prop_case <- x$prop_case - x$prop_case_v
  # x$diff_prop_hosp <- x$prop_hosp - x$prop_hosp_v
  x$diff_prop_death <- x$prop_death - x$prop_death_v
  x$diff_DALYs_pp <- x$DALYs_pp - x$DALYs_pp_v
  
  x$cases_averted <- x$case - x$case_v
  # x$hosps_averted <- x$hosp - x$hosp_v
  x$deaths_averted <- x$death - x$death_v
  x$DALYs_averted <- x$DALYs - x$DALYs_v
  
  return(x)
}

calc_infections_cases_deaths_DALYs_averted <- function(x){
  # Calculate infections, cases, deaths and DALYs averted
  x[,`:=`(deaths_averted = death - death_v,DALYs_averted = DALYs - DALYs_v,infections_averted = infection - infection_v,cases_averted = case - case_v)]
  
  # Calculate proportion of infections, cases, deaths and DALYs averted
  x[,`:=`(diff_prop_death = deaths_averted/population,diff_DALYs_pp = DALYs_averted/population,diff_prop_infections = infections_averted/population,diff_prop_cases = cases_averted/population)]
}

optimize_allocation <- function(x,n_v,by=""){  
  # optimize_allocation <- function(df,n_v){
  #   df$population <- 1
  #   x <- aggregate(cbind(population,case,hosp,death,case_v,hosp_v,death_v) ~ county_res + age_cat + sex + race_ethnicity + special.population,df,sum)
  
  # x$prop_case <- x$case/x$population
  # x$prop_hosp <- x$hosp/x$population
  # x$prop_death <- x$death/x$population
  # x$prop_case_v <- x$case_v/x$population
  # x$prop_hosp_v <- x$hosp_v/x$population
  # x$prop_death_v <- x$death_v/x$population
  # 
  # x$diff_prop_case <- x$prop_case - x$prop_case_v
  # x$diff_prop_hosp <- x$prop_hosp - x$prop_hosp_v
  # x$diff_prop_death <- x$prop_death - x$prop_death_v
  
  # Sort in descending order of difference in proportion of deaths/DALYs
  if(by==""){
    x <- x[order(-x$diff_prop_death),]  
  } else if (by=="DALYs"){
    x <- x[order(-x$diff_DALYs_pp),]
  }
  
  # Overwrite any missing population totals with 0s
  x$population[is.na(x$population)] <- 0
  # Calculate cumulative population
  x$cum_pop <- cumsum(x$population)
  # idx <- (x$cum_pop < n_v)
  # x$cases_averted[idx] <- x$case[idx] - x$case_v[idx]
  # x$hosps_averted[idx] <- x$hosp[idx] - x$hosp_v[idx]
  # x$deaths_averted[idx] <- x$death[idx] - x$death_v[idx]
  
  return(x)
}

optimize_allocation_by_grp <- function(x,n_v,grp,by=""){
  if (by==""){
    agg_x <- aggregate_by_grp(x,"diff_prop_death",grp,FUN = "mean")
    rank <- order(order(-agg_x$x))
    x$rank <- rank[match(x[,grp],agg_x[,grp])]
    x <- x[order(x$rank,-x$diff_prop_death),]  
  } else if (by=="DALYs"){
    agg_x <- aggregate_by_grp(x,"diff_DALYs_pp",grp,FUN = "mean")
    rank <- order(order(-agg_x$x))
    x$rank <- rank[match(x[,grp],agg_x[,grp])]
    x <- x[order(x$rank,-x$diff_DALYs_pp),]  
  }
  
  # Overwrite any missing population totals with 0s
  x$population[is.na(x$population)] <- 0
  # Calculate cumulative population
  x$cum_pop <- cumsum(x$population)
  # idx <- (x$cum_pop <= n_v)
  # res <- apply(x[idx,c("cases_averted","hosps_averted","deaths_averted")],2,sum)
  # 
  # return(res)
  return(x)
}

order_by_risk <- function(x,by=""){
  if (by==""){
    ord <- x[order(-lambda_adj),grp]
  } else if (by=="DALYs"){
    ord <- x[order(-lambda_DALYs),grp]
  }
  return(ord)
}

order_by_risk_by_grp <- function(x,grp1,by=""){
  if (by==""){
    agg_x <- x[,.(lambda_adj=mean(lambda_adj)),by=grp1]
    agg_x[,rank := order(order(-lambda_adj))]
    x <- merge(x,agg_x[,!"lambda_adj"],by=grp1,all.x=T)
    # ord <- x[order(rank,-lambda_adj),grp]
    ord <- x[order(rank,sample.int(nrow(x))),grp]
  } else if (by=="DALYs"){
    agg_x <- x[,.(lambda_DALYs=mean(lambda_DALYs)),by=grp1]
    agg_x[,rank := order(order(-lambda_DALYs))]
    x <- merge(x,agg_x[,!"lambda_DALYs"],by=grp1,all.x=T)
    ord <- x[order(rank),grp]
    # ord <- x[order(rank,-lambda_adj),grp]
    ord <- x[order(rank,sample.int(nrow(x))),grp]
  }
  return(ord)
}

aggregate_by_grp <- function(x,vrble,grp,FUN = "sum"){
  if (length(grp)==1){
    agg_x <- aggregate(x[,vrble],by = list(x[,grp]),FUN = FUN)
  } else if (length(grp)==2){
    agg_x <- aggregate(as.vector(x[,vrble]),by = list(x[,grp[1]],x[,grp[2]]),FUN = FUN)
    names(agg_x)[names(agg_x)=="Group.2"] <- grp[2]  
  }
  names(agg_x)[names(agg_x)=="Group.1"] <- grp[1]
  return(agg_x)
}

reshape_and_plot <- function(x,grp,fdir,vrble=c("case","hosp","death")){
  vrble1 <- c(vrble,paste0(vrble,"_v"))
  if (is.null(grp)){
    agg_x_long <- data.frame(variable = vrble1,value = apply(x[,vrble1],2,sum))
  } else {
    agg_x <- aggregate_by_grp(x,vrble1,grp)
    agg_x_long <- melt(agg_x,id.vars = grp)
  }
  for (i in 1:length(vrble)){
    last_char <- substr(vrble[i],nchar(vrble[i]),nchar(vrble[i]))
    ylbl <- ifelse(last_char=="s",vrble[i],paste0(vrble[i],"s"))
    p <- ggplot(agg_x_long[grep(vrble[i],agg_x_long$variable),],aes(fill=variable,x=agg_x_long[grep(vrble[i],agg_x_long$variable),1],y=value)) + geom_bar(position="dodge",stat="identity") + xlab(grp) + ylab(ylbl) + scale_fill_discrete(name = "",labels = c("no vaccination","vaccination"))
    if (is.null(grp)){
      p <- p + theme(axis.text.x = element_blank())
    } else {
      if (grp == "county_res"){
        p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1)) 
      }
    }
    pdf(paste0(fdir,"pred_",ylbl,"_by_",grp,".pdf"),width = 7.5, height = 4)
    print(p)
    dev.off()
  }
}

plot_predictions <- function(x,n_v,fdir){
  # Overwrite cases, hospitalisations and deaths under vaccination for those not vaccinated due to limited doses by corresponding values under no vaccination
  idx <- (x$cum_pop <= n_v)
  x$case_v[!idx] <- x$case[!idx]
  x$hosp_v[!idx] <- x$hosp[!idx]
  x$death_v[!idx] <- x$death[!idx]
  dir.create(fdir,recursive = T)
  reshape_and_plot(x,NULL,fdir)
  reshape_and_plot(x,"county_res",fdir)
  reshape_and_plot(x,"age_cat",fdir)
  reshape_and_plot(x,"sex",fdir)
  reshape_and_plot(x,"race_ethnicity",fdir)
  spec_pops <- c("non_spec_pop","HCW","prisoner","SNF","teacher","homeless")
  x$spec_pop <- spec_pops[x$special.population+1]
  reshape_and_plot(x,"spec_pop",fdir)
}

plot_predictions1 <- function(x,n_v,fdir){
  # Overwrite cases, deaths and DALYs under vaccination for those not vaccinated due to limited doses by corresponding values under no vaccination
  idx <- (x$cum_pop <= n_v)
  x$case_v[!idx] <- x$case[!idx]
  x$death_v[!idx] <- x$death[!idx]
  x$DALYs_v[!idx] <- x$DALYs[!idx]
  dir.create(fdir,recursive = T)
  vrble <- c("case","death","DALYs")
  reshape_and_plot(x,NULL,fdir,vrble)
  reshape_and_plot(x,"county_res",fdir,vrble)
  reshape_and_plot(x,"age_cat",fdir,vrble)
  reshape_and_plot(x,"sex",fdir,vrble)
  reshape_and_plot(x,"race_ethnicity",fdir,vrble)
  spec_pops <- c("non_spec_pop","HCW","prisoner","SNF","education","homeless","EW")
  x$spec_pop <- spec_pops[x$special.population+1]
  reshape_and_plot(x,"spec_pop",fdir,vrble)
}

reshape_and_plot_four_week_ahead <- function(x,grp,fdir){
  dir.create(fdir,recursive = T)
  if (is.null(grp)){
    agg_x_long <- data.frame(variable = c("cum_cases_pred","cum_hosps_pred","cum_deaths_pred"),value = apply(x[,c("cum_cases_pred","cum_hosps_pred","cum_deaths_pred")],2,sum))
  } else {
    agg_x <- aggregate(x[,c("cum_cases_pred","cum_hosps_pred","cum_deaths_pred")],by = list(x[,grp]),FUN = sum)
    names(agg_x)[names(agg_x)=="Group.1"] <- grp
    agg_x_long <- melt(agg_x,id.vars = grp)
  }
  pdf(paste0(fdir,"pred_cases_by_",grp,".pdf"),width = 5.25, height = 3)
  p <- ggplot(agg_x_long[agg_x_long$variable=="cum_cases_pred",],aes(x=agg_x_long[agg_x_long$variable=="cum_cases_pred",1],y=value)) + geom_bar(stat="identity") + ylab("cases") + theme(axis.title.x = element_blank())
  if (is.null(grp)){
    p <- p + theme(axis.text.x = element_blank())
  } else {
    if (grp == "county_res"){
      p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1)) 
    }
  }
  print(p)
  dev.off()
  pdf(paste0(fdir,"pred_hosps_by_",grp,".pdf"),width = 5.25, height = 3)
  p <- ggplot(agg_x_long[agg_x_long$variable=="cum_hosps_pred",],aes(x=agg_x_long[agg_x_long$variable=="cum_hosps_pred",1],y=value)) + geom_bar(stat="identity") + ylab("hospitalisations") + theme(axis.title.x = element_blank())
  if (is.null(grp)){
    p <- p + theme(axis.text.x = element_blank())
  } else {
    if (grp == "county_res"){
      p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1)) 
    }
  }
  print(p)
  dev.off()
  pdf(paste0(fdir,"pred_deaths_by_",grp,".pdf"),width = 5.25, height = 3)
  p <- ggplot(agg_x_long[agg_x_long$variable=="cum_deaths_pred",],aes(x=agg_x_long[agg_x_long$variable=="cum_deaths_pred",1],y=value)) + geom_bar(stat="identity") + ylab("deaths") + theme(axis.title.x = element_blank())
  if (is.null(grp)){
    p <- p + theme(axis.text.x = element_blank())
  } else {
    if (grp == "county_res"){
      p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1)) 
    }
  }
  print(p)
  dev.off()
}

plot_four_week_ahead_predictions <- function(x,RR,fdir){
  # Overwrite any missing values for predicted cases, hospitalisations and deaths with 0s
  x$cum_cases_pred[is.na(x$cum_cases_pred)] <- 0
  x$cum_hosps_pred[is.na(x$cum_hosps_pred)] <- 0
  x$cum_deaths_pred[is.na(x$cum_deaths_pred)] <- 0
  
  # # Calculate RR for non-special population
  # x$non_spec_pop <- x$population - x$hcw - x$prisoner - x$snf - x$teacher-x$homeless
  # x$RR0 <- x$population/(x$non_spec_pop + RR["HCW"]*x$hcw + RR["prisoner"]*x$prisoner + RR["SNF"]*x$snf + RR["teacher"]*x$teacher + RR["homeless"]*x$homeless)
  # 
  # # Adjust relative risks for special populations based on RR for non-special population
  # x$RR1 <- x$RR0 * RR["HCW"]
  # x$RR2 <- x$RR0 * RR["prisoner"]
  # x$RR3 <- x$RR0 * RR["SNF"]
  # x$RR4 <- x$RR0 * RR["teacher"]
  # x$RR5 <- x$RR0 * RR["homeless"]
  # 
  # # Calculate predicted cases in each special population
  # # x$cases0 <- round(x$RR0*x$non_spec_pop/x$population * x$cum_cases_pred)
  # x$cases1 <- round(x$RR1 * x$hcw/x$population * x$cum_cases_pred)
  # x$cases2 <- round(x$RR2 * x$prisoner/x$population * x$cum_cases_pred)
  # x$cases3 <- round(x$RR3 * x$snf/x$population * x$cum_cases_pred)
  # x$cases4 <- round(x$RR4 * x$teacher/x$population * x$cum_cases_pred)
  # x$cases5 <- round(x$RR5 * x$homeless/x$population * x$cum_cases_pred)
  # # Assign remainder of cases to non-special population as it is generally the largest subgroup so rounding makes least difference
  # x$cases0 <- x$cum_cases_pred - apply(x[,paste0("cases",1:5)],1,sum)
  # 
  # # Calculate predicted hospitalisations in each special population
  # # x$hosps0 <- round(x$RR0*x$non_spec_pop/x$population * x$cum_hosps_pred)
  # x$hosps1 <- round(x$RR1 * x$hcw/x$population * x$cum_hosps_pred)
  # x$hosps2 <- round(x$RR2 * x$prisoner/x$population * x$cum_hosps_pred)
  # x$hosps3 <- round(x$RR3 * x$snf/x$population * x$cum_hosps_pred)
  # x$hosps4 <- round(x$RR4 * x$teacher/x$population * x$cum_hosps_pred)
  # x$hosps5 <- round(x$RR5 * x$homeless/x$population * x$cum_hosps_pred)
  # # Assign remainder of hospitalisations to non-special population as it is generally the largest subgroup so rounding makes least difference
  # x$hosps0 <- x$cum_hosps_pred - apply(x[,paste0("hosps",1:5)],1,sum)
  # 
  # # Calculate predicted deaths in each special population
  # # x$cum_deaths0 <- round(x$RR0*x$non_spec_pop/x$population * x$cum_deaths_pred)
  # x$deaths1 <- round(x$RR1 * x$hcw/x$population * x$cum_deaths_pred)
  # x$deaths2 <- round(x$RR2 * x$prisoner/x$population * x$cum_deaths_pred)
  # x$deaths3 <- round(x$RR3 * x$snf/x$population * x$cum_deaths_pred)
  # x$deaths4 <- round(x$RR4 * x$teacher/x$population * x$cum_deaths_pred)
  # x$deaths5 <- round(x$RR5 * x$homeless/x$population * x$cum_deaths_pred)
  # # Assign remainder of deaths to non-special population as it is generally the largest subgroup so rounding makes least difference
  # x$deaths0 <- x$cum_deaths_pred - apply(x[,paste0("deaths",1:5)],1,sum)
  # 
  # # Assign numbers to each risk group
  # x$grp <- 1:nrow(x)
  # 
  # # x_long <- melt(x,id.vars = c("grp","county_res","age_cat","sex","race_ethnicity","lambda"),measure.vars = c(paste0("cum_cases",0:5),paste0("RR",0:5)),variable.name = "special.population",value.name = "RR")
  # # x_long$special.population <- as.numeric(sub("RR","",x_long$special.population))
  # x_long <- reshape(x[,c("county_res","age_cat","sex","race_ethnicity","time","lambda",paste0("RR",0:5),paste0("cases",0:5),paste0("hosps",0:5),paste0("deaths",0:5),"grp")],varying = list(paste0("RR",0:5),paste0("cases",0:5),paste0("hosps",0:5),paste0("deaths",0:5)),v.names = c("RR","cum_cases_pred","cum_hosps_pred","cum_deaths_pred"),timevar = "special.population",idvar = "grp",direction = "long")
  # x_long$special.population <- x_long$special.population - 1
  # x_long$risk_grp <- 1:nrow(x_long)
  
  reshape_and_plot_four_week_ahead(x,NULL,fdir)
  reshape_and_plot_four_week_ahead(x,"county_res",fdir)
  reshape_and_plot_four_week_ahead(x,"age_cat",fdir)
  reshape_and_plot_four_week_ahead(x,"sex",fdir)
  reshape_and_plot_four_week_ahead(x,"race_ethnicity",fdir)
  # spec_pops <- c("non_spec_pop","HCW","incarcerated","SNF","teacher","homeless")
  # x_long$spec_pop <- spec_pops[x_long$special.population+1]
  # reshape_and_plot_four_week_ahead(x_long,"spec_pop",fdir)
}

reshape_and_plot_all_strtgs <- function(x,grp,fdir,vrble=c("case","hosp","death"),xlbl=NULL,lbls=c("no vaccination","random","special population","geographic","age"),plt="Spectral"){
  vrble1 <- c(vrble,paste0(vrble,"_v"))
  if (grp=="strategy"){
    agg_x <- aggregate_by_grp(x,vrble1,"strategy")
    agg_x_long <- melt(agg_x,id.vars = "strategy")
    # Remove duplicated "no vaccination rows"
    agg_x_long <- agg_x_long[!(duplicated(agg_x_long[,-1]) & (agg_x_long$variable %in% vrble)),]
    # Change strategy number for no vaccination to 0
    agg_x_long$strategy[agg_x_long$variable %in% vrble] <- 0
  } else {
    agg_x <- aggregate_by_grp(x,vrble1,c("strategy",grp))
    agg_x_long <- melt(agg_x,id.vars = c("strategy",grp))
    # Remove duplicated "no vaccination rows"
    agg_x_long <- agg_x_long[!(duplicated(agg_x_long[,-1]) & (agg_x_long$variable %in% vrble)),]
    # Change strategy number for no vaccination to 0
    agg_x_long$strategy[agg_x_long$variable %in% vrble] <- 0
  }
  for (i in 1:length(vrble)){
    last_char <- substr(vrble[i],nchar(vrble[i]),nchar(vrble[i]))
    ylbl <- ifelse(last_char=="s",vrble[i],paste0(vrble[i],"s"))
    p <- ggplot(agg_x_long[grep(vrble[i],agg_x_long$variable),],aes(x=as.factor(agg_x_long[grep(vrble[i],agg_x_long$variable),grp]),y=value,fill=as.factor(strategy))) + geom_bar(position="dodge",stat="identity") + ylab(ylbl)
    if (is.null(xlbl)){
      p <- p + xlab(grp)
    } else {
      p <- p + xlab(xlbl)
    }
    if (grp=="strategy"){
      p <- p + theme(legend.position = "none") + scale_x_discrete(breaks=0:max(agg_x_long$strategy),labels = lbls) + scale_fill_brewer(palette = plt)
      w <- 5.5
    } else {
      p <- p + scale_fill_brewer(name = "Strategy",labels = lbls,palette = plt)
      if (grp == "county_res"){
        p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1)) 
      }
      w <- 7.5
    }
    pdf(paste0(fdir,"pred_",ylbl,"_by_",grp,".pdf"),width = w, height = 4)
    print(p)
    dev.off()    
  }
}

plot_all_strtgs <- function(x,grp,fdir,vrble=c("case","death","DALYs"),xlbl=NULL,lbls=c("no vaccination","random","special populations","age","essential workers","comorbidities"),plt="Spectral"){
  vrble1 <- paste0(vrble,"_v")
  if (grp=="strategy"){
    x_long <- melt(x,id.vars = "strategy")
    vrble <- ifelse(sapply(1:length(vrble),function(i){substr(vrble[i],nchar(vrble[i]),nchar(vrble[i]))})=="s",vrble,paste0(vrble,"s"))
  } else {
    x_long <- melt(x,id.vars = c("strategy",grp))
    # Remove duplicated "no vaccination rows"
    x_long <- x_long[!(duplicated(x_long[,-c("strategy","value")]) & (variable %in% paste0(rep(vrble,each=3),c("","_q_90_LB","_q_90_UB")))),]
    # Change strategy number for no vaccination to 0
    x_long[variable %in% paste0(rep(vrble,each=3),c("","_q_90_LB","_q_90_UB")),strategy:=0]
  }
  for (i in 1:length(vrble)){
    last_char <- substr(vrble[i],nchar(vrble[i]),nchar(vrble[i]))
    ylbl <- ifelse(last_char=="s",vrble[i],paste0(vrble[i],"s"))
    q_long <- x_long[variable %in% paste0(rep(c(vrble[i],vrble1[i]),each=3),c("","_q_90_LB","_q_90_UB")),]
    q_long[variable==vrble1[i],variable:=vrble[i]]
    q_long[variable==paste0(vrble1[i],"_q_90_LB"),variable:=paste0(vrble[i],"_q_90_LB")]
    q_long[variable==paste0(vrble1[i],"_q_90_UB"),variable:=paste0(vrble[i],"_q_90_UB")]
    q <- dcast(q_long,... ~ variable,value.var = "value")
    p <- ggplot(q,aes(x=as.factor(q[[grp]]),y=q[[vrble[i]]],fill=as.factor(strategy))) + geom_bar(position="dodge",stat="identity") + ylab(ylbl) + geom_errorbar(aes(ymin=q[[paste0(vrble[i],"_q_90_LB")]],ymax=q[[paste0(vrble[i],"_q_90_UB")]]),width=0.3,position=position_dodge(0.9))
    if (is.null(xlbl)){
      p <- p + xlab(grp)
    } else {
      p <- p + xlab(xlbl)
    }
    if (grp=="strategy"){
      p <- p + theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1)) + scale_x_discrete(breaks=0:max(x_long[,strategy]),labels = lbls) + scale_fill_brewer(palette = plt)
      w <- 5.5
    } else {
      p <- p + scale_fill_brewer(name = "Strategy",labels = lbls,palette = plt)
      if (grp == "county_res"){
        p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1)) 
      }
      w <- 7.5
    }
    pdf(paste0(fdir,"pred_",ylbl,"_by_",grp,".pdf"),width = w, height = 4)
    print(p)
    dev.off()    
  }
}

plot_predictions_all_strtgs <- function(x,n_v,fdir){
  # Overwrite cases, hospitalisations and deaths under vaccination for those not vaccinated due to limited doses by corresponding values under no vaccination
  idx <- (x$cum_pop <= n_v)
  x$case_v[!idx] <- x$case[!idx]
  x$hosp_v[!idx] <- x$hosp[!idx]
  x$death_v[!idx] <- x$death[!idx]
  dir.create(fdir,recursive = T)
  reshape_and_plot_all_strtgs(x,"strategy",fdir)
  reshape_and_plot_all_strtgs(x,"county_res",fdir)
  reshape_and_plot_all_strtgs(x,"age_cat",fdir)
  reshape_and_plot_all_strtgs(x,"sex",fdir)
  reshape_and_plot_all_strtgs(x,"race_ethnicity",fdir)
  spec_pops <- c("non_spec_pop","HCW","prisoner","SNF","education","homeless")
  x$spec_pop <- spec_pops[x$special.population+1]
  reshape_and_plot_all_strtgs(x,"spec_pop",fdir)
}

plot_predictions_all_strtgs1 <- function(x,n_v,fdir,lbls,plt){
  # Overwrite cases, deaths and DALYs under vaccination for those not vaccinated due to limited doses by corresponding values under no vaccination
  idx <- (x$cum_pop <= n_v)
  x$case_v[!idx] <- x$case[!idx]
  x$death_v[!idx] <- x$death[!idx]
  x$DALYs_v[!idx] <- x$DALYs[!idx]
  dir.create(fdir,recursive = T)
  vrble <- c("case","death","DALYs")
  reshape_and_plot_all_strtgs(x,"strategy",fdir,vrble,"Strategy",lbls,plt)
  reshape_and_plot_all_strtgs(x,"county_res",fdir,vrble,"County",lbls,plt)
  reshape_and_plot_all_strtgs(x,"age_cat",fdir,vrble,"Age",lbls,plt)
  reshape_and_plot_all_strtgs(x,"sex",fdir,vrble,"Sex",lbls,plt)
  reshape_and_plot_all_strtgs(x,"race_ethnicity",fdir,vrble,"Race/ethnicity",lbls,plt)
  spec_pops <- c("non_spec_pop","HCW","prisoner","SNF","education","homeless","EW")
  x$spec_pop <- spec_pops[x$special.population+1]
  reshape_and_plot_all_strtgs(x,"spec_pop",fdir,vrble,"Special population",lbls,plt)
}

plot_predictions_all_strtgs2 <- function(res_ss,lx_ss,fdir,vrble,lbls,plt,strategies1,grp1,lbls1){
  dir.create(fdir,recursive = T)
  plot_all_strtgs(res_ss,"strategy",fdir,vrble,"Strategy",lbls,plt)
  for (i in 1:length(lx_ss)){
    plot_all_strtgs(lx_ss[[i]][strategy %in% strategies1],grp1[i],fdir,vrble,lbls1[i],lbls,plt)
  }
}

plot_predictions_vacc_eff_SA <- function(res,strategies1,v_e_nms,vrble,fdir){
  dir.create(fdir,recursive = T)
  setDT(res)
  res_long <- melt(res[strategy %in% c(0,strategies1),],id.vars = c("strategy","v_e"))
  res_long$v_e <- factor(res_long$v_e,levels = v_e_nms[c(2,1,3)])
  
  for (i in 1:length(vrble)){
    last_char <- substr(vrble[i],nchar(vrble[i]),nchar(vrble[i]))
    vrble[i] <- ifelse(last_char=="s",vrble[i],paste0(vrble[i],"s"))
    q_long <- res_long[variable %in% paste0(rep(vrble[i],each=3),c("","_q_90_LB","_q_90_UB")),]
    q <- dcast(q_long,... ~ variable,value.var = "value")
    p <- ggplot(q,aes(x=as.factor(strategy),y=q[[vrble[i]]],fill=as.factor(v_e))) + geom_bar(position="dodge",stat="identity") + xlab("Strategy") + ylab(vrble[i]) + geom_errorbar(aes(ymin=q[[paste0(vrble[i],"_q_90_LB")]],ymax=q[[paste0(vrble[i],"_q_90_UB")]]),width=0.3,position=position_dodge(0.9)) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_x_discrete(breaks=0:max(res_long[,strategy]),labels = lbls) + scale_fill_brewer(name = "Efficacy",labels = lbls2,palette = "YlOrRd")
    pdf(paste0(fdir,"pred_",vrble[i],"_by_strategy_vacc_eff_SA.pdf"),width = 5.5, height = 4)
    print(p)
    dev.off()
  }
}

plot_vacc_distn <- function(x,v,vrble,grp,fdir,xlbl=NULL,ylbl=NULL){
  agg_v <- aggregate_by_grp(v,vrble,grp)
  agg_x <- aggregate_by_grp(x,vrble,grp)
  agg_v <- merge(agg_v,agg_x,by=grp,all.x=T)
  agg_v$prop <- agg_v$x.x/agg_v$x.y
  # agg_v <- agg_v[agg_v$prop>0.001,]
  # if (length(grp)==1){
  #   p <- ggplot(agg_v,aes(x=agg_v[,grp],y=x)) + geom_bar(stat="identity") + xlab(grp) + ylab("vaccinations")
  #   if (grp == "county_res"){
  #     p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1))
  #   }
  # } else if (length(grp)==2) {
  #   p <- ggplot(agg_v,aes(x=agg_v[,grp[1]],y=agg_v[,grp[2]],fill=x)) + geom_tile() + xlab(grp[1]) + ylab(grp[2]) + labs(fill="vaccinations")
  # }
  # pdf(paste0(fdir,"vacc_distn_by_",paste(grp,collapse = "_"),".pdf"),width = 7, height = 4)
  # print(p)
  # dev.off()
  if (length(grp)==1){
    p <- ggplot(agg_v,aes(x=as.factor(agg_v[,grp]),y=prop)) + geom_bar(stat="identity") + ylab("proportion vaccinated")
    if (is.null(xlbl)){
      p <- p + xlab(grp[1])
    } else {
      p <- p + xlab(xlbl)
    }
    if (grp %in% c("county_res","spec_pop","race_ethnicity")){
      p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1))
    }
    w <- ifelse(grp=="county_res",7,5)
  } else if (length(grp)==2) {
    p <- ggplot(agg_v,aes(x=as.factor(agg_v[,grp[1]]),y=as.factor(agg_v[,grp[2]]),fill=prop)) + geom_tile() + labs(fill="proportion vaccinated") + scale_fill_viridis(discrete=FALSE,limits=c(0,1))
    if (is.null(xlbl)){
      p <- p + xlab(grp[1]) + ylab(grp[2])
    } else {
      p <- p + xlab(xlbl) + ylab(ylbl)
    }
    if (grp[1] == "county_res"){
      p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1))
    }
    w <- 7
  }
  pdf(paste0(fdir,"vacc_propn_by_",paste(grp,collapse = "_"),".pdf"),width = w, height = 4)
  print(p)
  dev.off()
}

plot_optimal_allocation <- function(x,n_v,fdir,spec_pops=c("non_spec_pop","HCW","prisoner","SNF","education","homeless")){
  idx <- (x$cum_pop < n_v)
  v <- x[idx,]
  dir.create(fdir,recursive = T)
  plot_vacc_distn(x,v,"population","county_res",fdir,"County")
  plot_vacc_distn(x,v,"population","age_cat",fdir,"Age")
  plot_vacc_distn(x,v,"population","sex",fdir,"Sex")
  plot_vacc_distn(x,v,"population","race_ethnicity",fdir,"Race/ethnicity")
  x$spec_pop <- spec_pops[x$special.population+1]
  v$spec_pop <- spec_pops[v$special.population+1]
  plot_vacc_distn(x,v,"population","spec_pop",fdir,"Special population")
  plot_vacc_distn(x,v,"population","comorb",fdir,"Comorbidity")
  
  # Two-way breakdowns
  # plot_vacc_distn(x,v,"population",c("county_res","age_cat"),fdir,"County","Age")
  plot_vacc_distn(x,v,"population",c("age_cat","county_res"),fdir,"Age","County")
  plot_vacc_distn(x,v,"population",c("age_cat","sex"),fdir,"Age","Sex")
  plot_vacc_distn(x,v,"population",c("age_cat","race_ethnicity"),fdir,"Age","Race/ethnicity")
  plot_vacc_distn(x,v,"population",c("age_cat","spec_pop"),fdir,"Age","Special population")
  plot_vacc_distn(x,v,"population",c("age_cat","comorb"),fdir,"Age","Comorbidity")
}

run_prioritisation_strategies <- function(x,n_v,fdir,by="",dir="",t_sim=180,plt="Spectral"){
  res_nms <- c("cases_averted","deaths_averted","DALYs_averted")
  fdir <- paste0(fdir,n_v,"doses/")
  
  # Strategy 0: no vaccination
  totals <- apply(x[,c("case","death","DALYs")],2,sum)
    
  # Strategy 1: random allocation
  set.seed(123)
  x1 <- x[sample.int(nrow(x)),]
  x1$population[is.na(x1$population)] <- 0
  x1$cum_pop <- cumsum(x1$population)
  idx <- (x1$cum_pop <= n_v)
  res1 <- apply(x1[idx,res_nms],2,sum)
  plot_predictions1(x1,n_v,fdir=paste0(fdir,"Strategy1_1",by,"/"))
  
  # Strategy 2: allocation by special population
  x2 <- optimize_allocation_by_grp(x,n_v,"special.population",by)
  idx <- (x2$cum_pop <= n_v)
  res2 <- apply(x2[idx,res_nms],2,sum)
  plot_predictions1(x2,n_v,fdir=paste0(fdir,"Strategy2_1",by,"/"))
  
  # Strategy 3: allocation by age group
  x3 <- optimize_allocation_by_grp(x,n_v,"age_cat",by)
  idx <- (x3$cum_pop <= n_v)
  res3 <- apply(x3[idx,res_nms],2,sum)
  plot_predictions1(x3,n_v,fdir=paste0(fdir,"Strategy3_1",by,"/"))
  
  # Strategy 4: optimal allocation
  x4 <- optimize_allocation(x,n_v,by)
  idx <- (x4$cum_pop <= n_v)
  res4 <- apply(x4[idx,res_nms],2,sum)
  
  # Combine results into single data frame
  res0 <- data.frame(cases_averted=0,deaths_averted=0,DALYs_averted=0)
  res <- data.frame(cbind(Strategy = 0:4,rbind(res0,res1,res2,res3,res4)))
  res$prop_cases_averted <- res$cases_averted/totals["case"]
  res$prop_deaths_averted <- res$deaths_averted/totals["death"]
  res$prop_DALYs_averted <- res$DALYs_averted/totals["DALYs"]
  res[,c("cases","deaths","DALYs")] <- t(totals-t(res[,res_nms]))
  write.csv(res,paste0("../Data/",dir,"pred_cases_deaths_DALYs_all_strategies_1_",n_v,"doses_",t_sim,"days",by,".csv"),row.names = F)
  
  # Combine allocations for all strategies into one data frame
  x1$rank <- NA
  x1$strategy <- 1
  x2$strategy <- 2
  x3$strategy <- 3
  x_all <- rbind(x1,x2,x3)
  
  # Plot all strategies in grouped bar plots
  lbls <- c("no vaccination","random","special population","age")
  plot_predictions_all_strtgs1(x_all,n_v,paste0(fdir,"Strategies1to3_1",by,"/"),lbls,plt)
  
  # Plot optimal vaccine allocation
  spec_pops <- c("non_spec_pop","HCW","prisoner","SNF","education","homeless","EW")
  plot_optimal_allocation(x4,n_v,paste0(fdir,"OptimalAllocation_1",by,"/"),spec_pops)
}

run_prioritisation_strategies1 <- function(x,n_v,fdir,by="",dir="",t_sim=180,plt="Spectral"){
  res_nms <- c("cases_averted","deaths_averted","DALYs_averted")
  fdir <- paste0(fdir,n_v,"doses/")
  
  x_HCW_SNF <- x[x$special.population %in% c(1,3),]
  x_other <- x[!(x$special.population %in% c(1,3)),]  
  
  # Strategy 0: no vaccination
  totals <- apply(x[,c("case","death","DALYs")],2,sum)
  
  # Strategy 1: random allocation
  set.seed(123)
  x1_HCW_SNF <- optimize_allocation(x_HCW_SNF,n_v,by)
  x1_other <- x_other[sample.int(nrow(x_other)),]
  x1_other$cum_pop <- 0
  x1 <- rbind(x1_HCW_SNF,x1_other)
  x1$population[is.na(x1$population)] <- 0
  x1$cum_pop <- cumsum(x1$population)
  idx <- (x1$cum_pop <= n_v)
  res1 <- apply(x1[idx,res_nms],2,sum)
  plot_predictions1(x1,n_v,fdir=paste0(fdir,"Strategy1_2",by,"/"))
  
  # Strategy 2: special populations
  x2_HCW_SNF <- x1_HCW_SNF
  x2_HCW_SNF$rank <- NA
  x2_other <- optimize_allocation_by_grp(x_other,n_v,"special.population",by)
  x2 <- rbind(x2_HCW_SNF,x2_other)
  x2$cum_pop <- cumsum(x2$population) # recalculate cumulative population to correct from splitting
  idx <- (x2$cum_pop <= n_v)
  res2 <- apply(x2[idx,res_nms],2,sum)
  plot_predictions1(x2,n_v,fdir=paste0(fdir,"Strategy2_2",by,"/"))
  
  # Strategy 3: age targeting
  x3_HCW_SNF <- x1_HCW_SNF
  x3_HCW_SNF$rank <- NA
  x3_other <- optimize_allocation_by_grp(x_other,n_v,"age_cat",by)
  x3 <- rbind(x3_HCW_SNF,x3_other)
  x3$cum_pop <- cumsum(x3$population) # recalculate cumulative population to correct from splitting
  idx <- (x3$cum_pop <= n_v)
  res3 <- apply(x3[idx,res_nms],2,sum)
  plot_predictions1(x3,n_v,fdir=paste0(fdir,"Strategy3_2",by,"/"))
  
  # Strategy 4: essential workers
  x4_HCW_SNF <- x1_HCW_SNF
  x4_EW <- optimize_allocation(x_other[x_other$special.population %in% c(4,6),],n_v,by)
  x4_other <- optimize_allocation(x_other[!(x_other$special.population %in% c(4,6)),],n_v,by)
  x4 <- rbind(x4_HCW_SNF,x4_EW,x4_other)
  x4$cum_pop <- cumsum(x4$population) # recalculate cumulative population to correct from splitting
  idx <- (x4$cum_pop <= n_v)
  res4 <- apply(x4[idx,res_nms],2,sum)
  plot_predictions1(x4,n_v,fdir=paste0(fdir,"Strategy4_2",by,"/"))
  
  # Strategy 5: optimal allocation
  x5 <- optimize_allocation(x,n_v,by)
  idx <- (x5$cum_pop <= n_v)
  res5 <- apply(x5[idx,res_nms],2,sum)
  
  # Combine results into single data frame
  res0 <- data.frame(cases_averted=0,deaths_averted=0,DALYs_averted=0)
  names(res0) <- res_nms
  res <- data.frame(cbind(Strategy = 0:5,rbind(res0,res1,res2,res3,res4,res5)))
  res$prop_cases_averted <- res$cases_averted/totals["case"]
  res$prop_deaths_averted <- res$deaths_averted/totals["death"]
  res$prop_DALYs_averted <- res$DALYs_averted/totals["DALYs"]
  res[,c("cases","deaths","DALYs")] <- t(totals-t(res[,res_nms]))
  write.csv(res,paste0(dir,"pred_cases_deaths_DALYs_all_strategies_2_",n_v,"doses_",t_sim,"days",by,".csv"),row.names = F)
  
  # Combine allocations for all strategies into one data frame
  x1$rank <- NA
  x1$strategy <- 1
  x2$strategy <- 2
  x3$strategy <- 3
  x4$rank <- NA
  x4$strategy <- 4
  x_all <- rbind(x1,x2,x3,x4)
  
  # Plot all strategies in grouped bar plots
  lbls <- c("no vaccination","random","special population","age","essential workers")
  plot_predictions_all_strtgs1(x_all,n_v,paste0(fdir,"Strategies1to4_2",by,"/"),lbls,plt)
  
  # Plot optimal vaccine allocation
  spec_pops <- c("non_spec_pop","HCW","prisoner","SNF","education","homeless","EW")
  plot_optimal_allocation(x5,n_v,paste0(fdir,"OptimalAllocation_2",by,"/"),spec_pops)
  
  return(list(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5))
}

comb <- function(x, ...) {
  lapply(seq_along(x),function(i){c(x[[i]],lapply(list(...),function(y){y[[i]]}))})
}

order_and_calc_impact <- function(x_HCW_SNF,x_other,strategy,agg_df_other,by,comorbs,n_v,res_nms){
  # ord <- order_by_risk(agg_df_HCW_SNF,by)
  # x_HCW_SNF <- x_HCW_SNF[match(ord,x_HCW_SNF$grp),]
  x_HCW_SNF <- x_HCW_SNF[sample.int(nrow(x_HCW_SNF)),]
  if (strategy==1){ # random
    x_other <- x_other[sample.int(nrow(x_other)),]
    x <- rbind(x_HCW_SNF,x_other)
  } else if (strategy==2){ # special population
    # ord <- order_by_risk_by_grp(agg_df_other,"special.population",by)
    # x_other <- x_other[match(ord,x_other$grp),]
    spec_pop_idx <- !(x_other$special.population %in% c(0,6,7))
    x_spec_pop <- x_other[spec_pop_idx,]
    x_other <- x_other[!spec_pop_idx,]
    agg_df_spec_pop <- agg_df_other[grp %in% x_spec_pop$grp,]
    ord <- order_by_risk_by_grp(agg_df_spec_pop,"special.population",by)
    x_spec_pop <- x_spec_pop[match(ord,x_spec_pop$grp),]
    x_other <- x_other[sample.int(nrow(x_other)),]
    x <- rbind(x_HCW_SNF,x_spec_pop,x_other)
  } else if (strategy==3){ # age
    ord <- order_by_risk_by_grp(agg_df_other,"age_cat",by)
    x_other <- x_other[match(ord,x_other$grp),]
    x <- rbind(x_HCW_SNF,x_other)
  } else if (strategy==4){ # essential worker
    EW_idx <- (x_other$special.population %in% c(6,7))
    x_EW <- x_other[EW_idx,]
    x_other <- x_other[!EW_idx,]
    x_EW <- x_EW[sample.int(nrow(x_EW)),]
    x_other <- x_other[sample.int(nrow(x_other)),]
    x <- rbind(x_HCW_SNF,x_EW,x_other) 
  } else if (strategy==5){ # comorbidity
    comorb_idx <- (rowSums2(as.matrix(x_other[,..comorbs]))!=0)
    x_comorb <- x_other[comorb_idx,]
    x_other <- x_other[!comorb_idx,]
    x_comorb <- x_comorb[sample.int(nrow(x_comorb)),]
    x_other <- x_other[sample.int(nrow(x_other)),]
    x <- rbind(x_HCW_SNF,x_comorb,x_other)
  } else if (strategy==6){ # optimal
    ord <- order_by_risk(agg_df_other,by)
    x_other <- x_other[match(ord,x_other$grp),]
    x <- rbind(x_HCW_SNF,x_other)
  }
  x$population[is.na(x$population)] <- 0
  x$cum_pop <- cumsum(x$population)
  idx <- (x$cum_pop <= n_v)
  x[!idx,c("infection_v","case_v","death_v","DALYs_v")] <- x[!idx,c("infection","case","death","DALYs")]
  x$strategy <- strategy
  res <- x[idx,lapply(.SD,sum),.SDcols = res_nms]
  return(list(x=x,res=res))
}

run_prioritisation_strategies2 <- function(x,n_v,fdir,agg_df,strategies,by="",comorbs=c("asthma","diabetes","smoker","heart.disease","heart.failure","hypertension","obesity"),dir="",t_sim=180,lbls=c("no vaccination","random","special populations","age","essential workers","comorbidities"),plt="Spectral",spec_pops=c("non-spec-pop","HCW","prisoner","SNF","education","homeless","frontline EW","non-frontline EW","ALF")){
  res_nms <- c("infections_averted","cases_averted","deaths_averted","DALYs_averted")
  # fdir <- paste0(fdir,n_v,"doses/")

  HCW_SNF_idx <- (x$special.population %in% c(1,3,8))
  x_HCW_SNF <- x[HCW_SNF_idx,]
  x_other <- x[!HCW_SNF_idx,]
  
  agg_df_HCW_SNF <- agg_df[grp %in% x_HCW_SNF$grp,]
  agg_df_other <- agg_df[grp %in% x_other$grp,]

  # Strategy 0: no vaccination
  totals <- apply(x[,c("infection","case","death","DALYs")],2,sum)
  
  # Prioritization strategies
  lx <- vector("list",length(strategies))
  lres <- vector("list",length(strategies))
  for (i in 1:length(strategies)){
    l <- order_and_calc_impact(x_HCW_SNF,x_other,strategies[i],agg_df_other,by,comorbs,n_v,res_nms)
    lx[[i]] <- l$x
    lres[[i]] <- l$res
  }
      
  # # Strategy 1: random allocation
  # set.seed(123)
  # # x1_HCW_SNF <- optimize_allocation(x_HCW_SNF,n_v,by)
  # ord1 <- order_by_risk(agg_df_HCW_SNF,by)
  # x1_HCW_SNF <- x_HCW_SNF[match(ord1,x_HCW_SNF$grp),]
  # x1_other <- x_other[sample.int(nrow(x_other)),]
  # # x1_other$cum_pop <- 0
  # x1 <- rbind(x1_HCW_SNF,x1_other)
  # x1$population[is.na(x1$population)] <- 0
  # x1$cum_pop <- cumsum(x1$population)
  # idx <- (x1$cum_pop <= n_v)
  # x1[!idx,c("infection_v","case_v","death_v","DALYs_v")] <- x1[!idx,c("infection","case","death","DALYs")]
  # # res1 <- apply(x1[idx,..res_nms],2,sum)
  # res1 <- x1[idx,lapply(.SD,sum),.SDcols = res_nms]
  # # plot_predictions1(x1,n_v,fdir=paste0(fdir,"Strategy1_2",by,"/"))
  # 
  # # Strategy 2: special populations
  # x2_HCW_SNF <- x1_HCW_SNF
  # # x2_HCW_SNF$rank <- NA
  # # x2_other <- optimize_allocation_by_grp(x_other,n_v,"special.population",by)
  # ord2 <- order_by_risk_by_grp(agg_df_other,"special.population",by)
  # x2_other <- x_other[match(ord2,x_other$grp),]
  # x2 <- rbind(x2_HCW_SNF,x2_other)
  # x2$population[is.na(x2$population)] <- 0
  # x2$cum_pop <- cumsum(x2$population) # recalculate cumulative population to correct from splitting
  # idx <- (x2$cum_pop <= n_v)
  # x2[!idx,c("infection_v","case_v","death_v","DALYs_v")] <- x2[!idx,c("infection","case","death","DALYs")]
  # # res2 <- apply(x2[idx,..res_nms],2,sum)
  # res2 <- x2[idx,lapply(.SD,sum),.SDcols = res_nms]
  # # plot_predictions1(x2,n_v,fdir=paste0(fdir,"Strategy2_2",by,"/"))
  # 
  # # Strategy 3: age targeting
  # x3_HCW_SNF <- x1_HCW_SNF
  # # x3_HCW_SNF$rank <- NA
  # # x3_other <- optimize_allocation_by_grp(x_other,n_v,"age_cat",by)
  # ord3 <- order_by_risk_by_grp(agg_df_other,"age_cat",by)
  # x3_other <- x_other[match(ord3,x_other$grp),]
  # x3 <- rbind(x3_HCW_SNF,x3_other)
  # x3$population[is.na(x3$population)] <- 0
  # x3$cum_pop <- cumsum(x3$population) # recalculate cumulative population to correct from splitting
  # idx <- (x3$cum_pop <= n_v)
  # x3[!idx,c("infection_v","case_v","death_v","DALYs_v")] <- x3[!idx,c("infection","case","death","DALYs")]
  # # res3 <- apply(x3[idx,..res_nms],2,sum)
  # res3 <- x3[idx,lapply(.SD,sum),.SDcols = res_nms]
  # # plot_predictions1(x3,n_v,fdir=paste0(fdir,"Strategy3_2",by,"/"))
  # 
  # # Strategy 4: essential workers
  # x4_HCW_SNF <- x1_HCW_SNF
  # # x4_EW <- optimize_allocation(x_other[x_other$special.population %in% c(4,6),],n_v,by)
  # # x4_other <- optimize_allocation(x_other[!(x_other$special.population %in% c(4,6)),],n_v,by)
  # EW_idx <- (x_other$special.population %in% c(4,6))
  # x4_EW <- x_other[EW_idx,]
  # x4_other <- x_other[!EW_idx,]
  # ord4_EW <- order_by_risk(agg_df_other[grp %in% x4_EW$grp,],by)
  # ord4_other <- order_by_risk(agg_df_other[grp %in% x4_other$grp,],by)
  # x4_EW <- x4_EW[match(ord4_EW,x4_EW$grp),]
  # x4_other <- x4_other[match(ord4_other,x4_other$grp),]
  # x4 <- rbind(x4_HCW_SNF,x4_EW,x4_other)
  # x4$population[is.na(x4$population)] <- 0
  # x4$cum_pop <- cumsum(x4$population) # recalculate cumulative population to correct from splitting
  # idx <- (x4$cum_pop <= n_v)
  # x4[!idx,c("infection_v","case_v","death_v","DALYs_v")] <- x4[!idx,c("infection","case","death","DALYs")]
  # # res4 <- apply(x4[idx,..res_nms],2,sum)
  # res4 <- x4[idx,lapply(.SD,sum),.SDcols = res_nms]
  # # plot_predictions1(x4,n_v,fdir=paste0(fdir,"Strategy4_2",by,"/"))
  # 
  # # Strategy 5: comorbidity targeting
  # x5_HCW_SNF <- x1_HCW_SNF
  # comorb_idx <- (rowSums2(as.matrix(x_other[,..comorbs]))!=0)
  # # x5_comorb <- optimize_allocation(x_other[comorb_idx,],n_v,by)
  # # x5_other <- optimize_allocation(x_other[!comorb_idx,],n_v,by)
  # x5_comorb <- x_other[comorb_idx,]
  # x5_other <- x_other[!comorb_idx,]
  # ord5_comorb <- order_by_risk(agg_df_other[grp %in% x5_comorb$grp],by)
  # ord5_other <- order_by_risk(agg_df_other[grp %in% x5_other$grp],by)
  # x5_comorb <- x5_comorb[match(ord5_comorb,x5_comorb$grp),]
  # x5_other <- x5_other[match(ord5_other,x5_other$grp),]
  # x5 <- rbind(x5_HCW_SNF,x5_comorb,x5_other)
  # x5$population[is.na(x5$population)] <- 0
  # x5$cum_pop <- cumsum(x5$population) # recalculate cumulative population to correct from splitting
  # idx <- (x5$cum_pop <= n_v)
  # x5[!idx,c("infection_v","case_v","death_v","DALYs_v")] <- x5[!idx,c("infection","case","death","DALYs")]
  # # res5 <- apply(x5[idx,..res_nms],2,sum)
  # res5 <- x5[idx,lapply(.SD,sum),.SDcols = res_nms]
  # # plot_predictions1(x5,n_v,fdir=paste0(fdir,"Strategy5_2",by,"/"))
  # 
  # # Strategy 6: optimal allocation
  # # x6 <- optimize_allocation(x,n_v,by)
  # ord6 <- order_by_risk(agg_df,by)
  # x6 <- x[match(ord6,x$grp),] 
  # x6$population[is.na(x6$population)] <- 0
  # x6$cum_pop <- cumsum(x6$population)
  # idx <- (x6$cum_pop <= n_v)
  # x6[!idx,c("infection_v","case_v","death_v","DALYs_v")] <- x6[!idx,c("infection","case","death","DALYs")]
  # # res6 <- apply(x6[idx,..res_nms],2,sum)
  # res6 <- x6[idx,lapply(.SD,sum),.SDcols = res_nms]
  
  # Combine results into single data frame
  res0 <- data.frame(matrix(0,nrow=1,ncol=length(res_nms)))
  names(res0) <- res_nms
  # res <- data.frame(cbind(strategy = 0:6,rbind(res0,res1,res2,res3,res4,res5,res6)))
  res <- data.frame(cbind(strategy = c(0,strategies),rbind(res0,do.call(rbind,lres))))
  res$prop_infections_averted <- res$infections_averted/totals["infection"]
  res$prop_cases_averted <- res$cases_averted/totals["case"]
  res$prop_deaths_averted <- res$deaths_averted/totals["death"]
  res$prop_DALYs_averted <- res$DALYs_averted/totals["DALYs"]
  res[,c("infections","cases","deaths","DALYs")] <- t(totals-t(res[,res_nms]))
  # write.csv(res,paste0("../Data/",dir,"pred_cases_deaths_DALYs_all_strategies_2_",n_v,"doses_",t_sim,"days",by,".csv"),row.names = F)
  
  # Combine allocations for all strategies into one data frame
  # # x1$rank <- NA
  # x1$strategy <- 1
  # x2$strategy <- 2
  # x3$strategy <- 3
  # # x4$rank <- NA
  # x4$strategy <- 4
  # # x5$rank <- NA
  # x5$strategy <- 5
  # x_all <- rbind(x1,x2,x3,x4,x5)
  x_all <- do.call(rbind,lx)
  # idx <- (x_all$cum_pop <= n_v)
  x_all$comorb <- as.integer(rowSums2(as.matrix(x_all[,..comorbs]))!=0)
  cols <- c("infection","case","death","DALYs","infection_v","case_v","death_v","DALYs_v")
  # x_all <- x_all[,lapply(.SD, sum), .SDcols = cols, by = .(strategy,county_res,age_cat,sex,race_ethnicity,special.population,comorb)]
  
  x_all_county <- x_all[,lapply(.SD, sum), .SDcols = cols, by = .(strategy,county_res)]
  x_all_age <- x_all[,lapply(.SD, sum), .SDcols = cols, by = .(strategy,age_cat)]
  x_all_sex <- x_all[,lapply(.SD, sum), .SDcols = cols, by = .(strategy,sex)]
  x_all_race_ethnicity <- x_all[,lapply(.SD, sum), .SDcols = cols, by = .(strategy,race_ethnicity)]
  x_all_spec_pop <- x_all[,lapply(.SD, sum), .SDcols = cols, by = .(strategy,special.population)]
  x_all_comorb <- x_all[,lapply(.SD, sum), .SDcols = cols, by = .(strategy,comorb)]
   
  # cols <- c("infection","case","death","DALYs","infection_v","case_v","death_v","DALYs_v")
  # x1$comorb <- as.numeric(rowSums2(as.matrix(x1[,..comorbs]))!=0)
  # x1 <- x1[,lapply(.SD, sum), .SDcols = cols, by = .(county_res,age_cat,sex,race_ethnicity,special.population,comorb)]
  # x2$comorb <- as.numeric(rowSums2(as.matrix(x2[,..comorbs]))!=0)
  # x2 <- x2[,lapply(.SD, sum), .SDcols = cols, by = .(county_res,age_cat,sex,race_ethnicity,special.population,comorb)]
  # x3$comorb <- as.numeric(rowSums2(as.matrix(x3[,..comorbs]))!=0)
  # x3 <- x3[,lapply(.SD, sum), .SDcols = cols, by = .(county_res,age_cat,sex,race_ethnicity,special.population,comorb)]
  # x4$comorb <- as.numeric(rowSums2(as.matrix(x4[,..comorbs]))!=0)
  # x4 <- x4[,lapply(.SD, sum), .SDcols = cols, by = .(county_res,age_cat,sex,race_ethnicity,special.population,comorb)]
  # x5$comorb <- as.numeric(rowSums2(as.matrix(x5[,..comorbs]))!=0)
  # x5 <- x5[,lapply(.SD, sum), .SDcols = cols, by = .(county_res,age_cat,sex,race_ethnicity,special.population,comorb)]
  # x6$comorb <- as.numeric(rowSums2(as.matrix(x6[,..comorbs]))!=0)
  # x6 <- x6[,lapply(.SD, sum), .SDcols = cols, by = .(county_res,age_cat,sex,race_ethnicity,special.population,comorb)]
  # 
  # # # Plot all strategies in grouped bar plots
  # # plot_predictions_all_strtgs1(x_all,n_v,paste0(fdir,"Strategies1to",length(lbls)-1,"_2",by,"/"),lbls,plt)
  # # 
  # # # Plot optimal vaccine allocation
  # # plot_optimal_allocation(x6,n_v,paste0(fdir,"OptimalAllocation_2",by,"/"),spec_pops)
  # # 
  # # return(list(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6))
  # return(list(res=res,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6))
  # return(res)
  return(list(res=res,x_all_county=x_all_county,x_all_age=x_all_age,x_all_sex=x_all_sex,x_all_race_ethnicity=x_all_race_ethnicity,x_all_spec_pop=x_all_spec_pop,x_all_comorb=x_all_comorb))
}

q_90_LB <- function(x){quantile(x,probs = 0.05)}

q_90_UB <- function(x){quantile(x,probs = 0.95)}

calc_summary_stats <- function(x,cols,grp1){
  x_mean <- x[, lapply(.SD, mean), .SDcols = cols, by = grp1]
  x_median <- x[, lapply(.SD, median), .SDcols = cols, by = grp1][, setnames(.SD, cols, paste(cols, "median", sep = "_"))]
  x_q_90_LB <- x[, lapply(.SD, q_90_LB), .SDcols = cols, by = grp1][, setnames(.SD, cols, paste(cols, "q_90_LB", sep = "_"))]
  x_q_90_UB <- x[, lapply(.SD, q_90_UB), .SDcols = cols, by = grp1][, setnames(.SD, cols, paste(cols, "q_90_UB", sep = "_"))]
  x_ss <- cbind(x_mean,x_median[,!grp1,with=F],x_q_90_LB[,!grp1,with=F],x_q_90_UB[,!grp1,with=F])
  return(x_ss)
}
