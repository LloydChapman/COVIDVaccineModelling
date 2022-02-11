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
  x$cases0 <- x$cum_cases - apply(x[,paste0("cases",1:6)],1,sum)
  
  # Assign numbers to each risk group
  x$grp <- 1:nrow(x)
  
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
}

add_risk_ests_deaths <- function(glm_fit,alpha,df,RR,RR_LB,RR_UB,HR,lt_long,seroprev,IFR_long,IFR_ratio,agg_deaths_age,frlty_idx,RR_death_essential){#
  # Add 1 for RR of infection for non-special population to start of RR vector
  RR <- c(non_spec_pop=1,RR)
  RR_LB <- c(non_spec_pop=1,RR_LB)
  RR_UB <- c(non_spec_pop=1,RR_UB)
  
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
  # x1 <- x[,c(demogrphcs,"time_death","log_lambda","se_log_lambda","log_lambda_LB","log_lambda_UB"),with=F]
  # setkeyv(x1,demogrphcs)
  # setkeyv(df,demogrphcs)
  # df[x1,(names(x1)):=mget(paste0("i.",names(x1)))]
  
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
  # y[,`:=`(RR = RR[special.population+1] * HR,RR_LB = RR_LB[special.population+1] * HR,RR_UB = RR_UB[special.population+1] * HR)]
  y[,`:=`(RR_inf = RR[special.population+1],RR_inf_LB = RR_LB[special.population+1],RR_inf_UB = RR_UB[special.population+1])]
  # # Calculate standard error for relative risk
  # y[,se_RR := pmin(RR-RR_LB,RR_UB-RR)/z_star]
  
  # Aggregate synthetic population data table by demographics, special populations, and comorbidities
  agg_df <- df[,.(population=.N),by=c(demogrphcs,spec_pop_comorbs)]
  
  # Make data table of all unique combinations of demographic factors, special populations, and comorbidities
  z <- CJ(county_res=unique(agg_df[,county_res]),age_cat=unique(agg_df[,age_cat]),sex=unique(agg_df[,sex]),race_ethnicity=unique(agg_df[,race_ethnicity]),special.population=unique(agg_df[,special.population]),asthma=c(0,1),diabetes=c(0,1),smoker=c(0,1),heart.disease=c(0,1),heart.failure=c(0,1),hypertension=c(0,1),obesity=c(0,1))
  
  # Merge with aggregated synthetic population data table
  agg_df <- merge(z,agg_df,all.x=T)
  agg_df[is.na(population),population:=0]
  
  # Merge with IFR
  agg_df[,sex:=ifelse(sex==1,"male","female")]
  agg_df <- merge(agg_df,IFR_long[,.(age_cat,sex,median_perc,CI_95_LB,CI_95_UB)],by=c("age_cat","sex"))

  # Calculate average IFR for each special population
  agg_df[,IFR_spec_pop_mean:=sum(population*median_perc)/(100*sum(population))*IFR_ratio,by=.(special.population)]
  # agg_df[special.population %in% c(3,8),IFR_spec_pop_mean:=frlty_idx*IFR_spec_pop_mean]
  agg_df[,RR_death_given_inf:=IFR_spec_pop_mean/agg_df[special.population==0,IFR_spec_pop_mean][1]]
  
  # Add group index for special population and comorbidity combinations
  agg_df[,grp:=.GRP,by=c(demogrphcs,spec_pop_comorbs)]
  # Merge with RR data table
  agg_df <- merge(agg_df,y,by=spec_pop_comorbs,all.x=T)
  
  # Merge with RR of death for LTCF residents
  agg_df <- merge(agg_df,agg_deaths_age[,.(age_cat,RR_LTCF,RR_LTCF_LB,RR_LTCF_UB)],by="age_cat",all.x=T)
  
  # Calculate RR of death
  agg_df[,`:=`(RR=RR_inf*RR_death_given_inf,RR_LB=RR_inf_LB*RR_death_given_inf,RR_UB=RR_inf_UB*RR_death_given_inf)]
  # Overwrite RR of death for LTCF residents
  agg_df[special.population %in% c(3,8),`:=`(RR=RR_LTCF,RR_LB=RR_LTCF_LB,RR_UB=RR_LTCF_UB)]
  # Overwrite RRs of death for essential workers
  agg_df[special.population==6,`:=`(RR=RR_death_essential[special.population==6,RR],RR_LB=RR_death_essential[special.population==6,RR_LB],RR_UB=RR_death_essential[special.population==6,RR_UB])]
  agg_df[special.population==7,`:=`(RR=RR_death_essential[special.population==7,RR],RR_LB=RR_death_essential[special.population==7,RR_LB],RR_UB=RR_death_essential[special.population==7,RR_UB])]
  # Multiply by HR for comorbidities
  agg_df[,`:=`(RR=RR*HR,RR_LB=RR_LB*HR,RR_UB=RR_UB*HR)]
  agg_df <- agg_df[,se_RR := pmin(RR-RR_LB,RR_UB-RR)/z_star]
  
  # Adjust RRs so that overall death rate remains the same
  # agg_df <- merge(agg_df,agg_df[,.(RR0 = sum(population)/sum(population*RR_inf)),by=demogrphcs],by=demogrphcs)
  agg_df[,RR0 := sum(population)/sum(population*RR),by=demogrphcs]
  agg_df <- agg_df[,`:=`(RR = RR*RR0,RR_LB = RR_LB*RR0,RR_UB = RR_UB*RR0,se_RR = se_RR*RR0)]
  x[,sex:=as.character(sex)][sex=="1",sex:="male"][sex=="0",sex:="female"]
  agg_df <- merge(agg_df,x[,c(demogrphcs,"cum_deaths"),with=F],by=demogrphcs,all.x=T)
  setnames(agg_df,"cum_deaths","cum_deaths_demogrphc")
  # Estimate cumulative deaths in each demographic-special-population-comorbidity group
  agg_df[,cum_deaths:=population*RR/sum(population*RR)*cum_deaths_demogrphc,by=demogrphcs]
  agg_df[is.na(cum_deaths),cum_deaths:=0]

  # Estimate past infections in each demographic-special-population-comorbidity group from deaths and IFR
  agg_df[,cum_inf:=cum_deaths/(median_perc/100*IFR_ratio)]
  # # Adjust estimate of cumulative infections in LTCFs by frailty index
  agg_df[special.population %in% c(3,8),cum_inf:=cum_inf/frlty_idx]
  # Calculate probability of past infection by age, sex and special population
  agg_df[,p_past_inf:=sum(cum_inf)/sum(population),by=.(age_cat,sex,special.population)]
  agg_df[p_past_inf>1,p_past_inf:=1]
  
  # Add risk group variable for demographic and special population combination
  agg_df[,risk_grp:=.GRP,by=c(demogrphcs,"special.population")]

  # Merge with synthetic population data table by demographic, special population and comorbidity combination
  df[,sex:=as.character(sex)][sex=="1",sex:="male"][sex=="0",sex:="female"]
  agg_df[,c("IFR_spec_pop_mean","RR_death_given_inf","RR_inf","RR_inf_LB","RR_inf_UB","RR_LTCF","RR_LTCF_LB","RR_LTCF_UB","cum_deaths","cum_deaths_demogrphc","cum_inf"):=NULL] # remove columns
  df <- merge(df,agg_df,by=c(demogrphcs,spec_pop_comorbs),all.x=T)
  # setkeyv(agg_df,c(demogrphcs,spec_pop_comorbs))
  # setkeyv(df,c(demogrphcs,spec_pop_comorbs))
  # df[agg_df,(names(agg_df)):=mget(paste0("i.",names(agg_df)))]
    
  # Merge with life table data frame
  df <- merge(df,lt_long,by=c("age","sex","race_ethnicity"),all.x=T)
  # setkeyv(lt_long,c("age","sex","race_ethnicity"))
  # setkeyv(df,c("age","sex","race_ethnicity"))
  # df[lt_long,(names(lt_long)):=mget(paste0("i.",names(lt_long)))]
  
  # Merge with seroprevalence data frame
  age_cat_sero_lbs <- c(as.numeric(gsub("-.*|\\+|>","",seroprev[,age_cat_sero])),max(df[,age])+1)
  df[,age_cat_sero:=cut(age,age_cat_sero_lbs,right = F,labels = seroprev[,age_cat_sero])]
  df <- merge(df,seroprev,by="age_cat_sero",all.x=T)
  # setkey(seroprev,age_cat_sero)
  # setkey(df,age_cat_sero)
  # df[seroprev,(names(seroprev)):=mget(paste0("i.",names(seroprev)))]
  
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

calc_QALYs_from_deaths <- function(death,LE,age_mid,Q_a){
  E <- mapply(function(x,y) sum(Q_a[(x+1):(x+round(y)+1)]),age_mid,LE) # add 1 to index as ages in Q_a start from 0
  QALYs <- death * E
  return(QALYs)
}

calc_DALYs_from_infections_and_cases <- function(infection,case,d,w){
  DALYs <- ((infection - case) * d[1] * w[1] + case * d[2] * w[2])/365
  return(DALYs)
}

calc_QALYs_from_infections_and_cases <- function(infection,case,u){
  QALYs <- (0.07*(infection - case) + case)*u
  return(QALYs)
}

simulate_cases_hosps_deaths <- function(df,v_e,t_sim,r,cases_ratio,deaths_ratio,n_sim,p_mi,d,w){
  # Total population
  N <- nrow(df)
  
  # p_past_case <- 1 - exp(-df$lambda*df$time[1] * df$RR)
  p_past_case <- df$seroprev * df$RR
  p_case <- 1 - exp(-df$lambda*r*(1-r^t_sim)/(1-r) * df$RR * cases_ratio)
  p_case_v <- (1 - exp(-(1-v_e)*df$lambda*r*(1-r^t_sim)/(1-r) * df$RR * cases_ratio))

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

    # Simulate deaths
    death <- (case & (df$p_death * deaths_ratio) > runif(N))
    DALYs <- calc_DALYs(death,df$life_expectancy,case,p_mi,d,w)

    ## Vaccination
    # Simulate new cases
    case_v <- (!past_case & p_case_v > runif(N))

    # Simulate deaths
    death_v <- (case_v & (df$p_death * deaths_ratio) > runif(N))
    DALYs_v <- calc_DALYs(death_v,df$life_expectancy,case_v,p_mi,d,w)
    
    out[,,i] <- cbind(past_case,case,death,DALYs,case_v,death_v,DALYs_v)
  }

  return(out)
}

simulate_deaths <- function(df,v_e,t_sim,r,deaths_ratio,n_sim,d,w){
  # Total population
  N <- nrow(df)
  
  p_past_inf <- df$p_past_inf
  # lambda <- exp(df$log_lambda)
  # RR <- df$RR
  log_lambda <- rtruncnorm(N,a = df$log_lambda_LB,b = df$log_lambda_UB,mean = df$log_lambda,sd = df$se_log_lambda)
  lambda <- exp(log_lambda)
  RR <- rtruncnorm(N,a = df$RR_LB - 1e-15,b = df$RR_UB + 1e-15,mean = df$RR,sd = df$se_RR)
  if (r!=1){
    p_death <- 1 - exp(-lambda*r*(1-r^t_sim)/(1-r) * RR * deaths_ratio)
    p_death_v <- (1 - exp(-(1-v_e)*lambda*r*(1-r^t_sim)/(1-r) * RR * deaths_ratio))    
  } else {
    p_death <- 1 - exp(-lambda*t_sim * RR * deaths_ratio)
    p_death_v <- (1 - exp(-(1-v_e)*lambda*t_sim * RR * deaths_ratio))    
  }
  
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
  out <- cbind(past_inf,death,DALYs,death_v,DALYs_v)
  
  return(out)
}

simulate_deaths2 <- function(df,v_e,t_sim,r,deaths_ratio,d,w,Q_a){
  # Total number of risk groups
  n_grp <- nrow(df)
  
  pop <- df$population
  p_past_inf <- df$p_past_inf
  log_lambda <- rtruncnorm(n_grp,a = df$log_lambda_LB,b = df$log_lambda_UB,mean = df$log_lambda,sd = df$se_log_lambda)
  lambda <- exp(log_lambda)
  RR <- rtruncnorm(n_grp,a = df$RR_LB - 1e-15,b = df$RR_UB + 1e-15,mean = df$RR,sd = df$se_RR)
  p_death <- 1 - exp(-lambda*t_sim * RR * deaths_ratio)
  p_death_v <- (1 - exp(-lambda*(t_sim-v_e/r*(1-exp(-r*t_sim))) * RR * deaths_ratio))
  p_death_i <- (1 - exp(-lambda*(t_sim-v_e/r*(1-exp(-r*t_sim))) * RR * deaths_ratio))
    
  ## No vaccination
  # Simulate past infections
  past_inf <- rbinom(n_grp,pop,p_past_inf)
  
  # Simulate new deaths
  death <- rbinom(n_grp,pop-past_inf,p_death)
  death_i <- rbinom(n_grp,past_inf,p_death_i)
  
  # Calculate DALYs from deaths
  DALYs <- calc_DALYs_from_deaths(death+death_i,df$life_expectancy,d,w)
  
  # Calculate QALYs lost from deaths
  QALYs <- calc_QALYs_from_deaths(death+death_i,df$life_expectancy,age_mid,Q_a)
  
  ## Vaccination
  # Simulate new deaths
  death_v <- rbinom(n_grp,pop-past_inf,p_death_v)
  
  # Calculate DALYs from deaths
  DALYs_v <- calc_DALYs_from_deaths(death_v+death_i,df$life_expectancy,d,w)
  
  # Calculate QALYs lost from deaths
  QALYs_v <- calc_QALYs_from_deaths(death_v+death_i,df$life_expectancy,age_mid,Q_a)
  
  # Bind into output array
  out <- cbind(past_inf,death,DALYs,QALYs,death_v,DALYs_v,QALYs_v,death_i)
  
  return(out)
}

aggregate_predictions <- function(out,df,n_sim){
  # Bind column with risk group ID
  out <- abind(risk_grp=array(rep(df$risk_grp,n_sim),dim=c(nrow(out),1,n_sim)),out,along=2)
  # out_nms <- colnames(out)
  
  # Aggregate outcomes by risk groups
  agg_out <- array(dim=c(length(unique(df$risk_grp)),dim(out)[2],dim(out)[3]))
  for (i in 1:n_sim){
    agg_out[,,i] <- as.matrix(aggregate(cbind(past_case,case,death,DALYs,case_v,death_v,DALYs_v) ~ risk_grp,out[,,i],sum))
  }
  return(agg_out)
}

aggregate_predictions_deaths <- function(out,grp,n_sim){
  # Bind column with risk group ID
  out <- cbind(grp,out)
  out_dt <- as.data.table(out)
  
  # Aggregate outcomes by risk groups
  agg_out <- out_dt[,.(past_inf=sum(past_inf),death=sum(death),DALYs=sum(DALYs),death_v=sum(death_v),DALYs_v=sum(DALYs_v)),by=.(grp)]
  
  return(agg_out)
}

aggregate_predictions_by_grp <- function(out,cols,grp1){
  out[,lapply(.SD, sum), .SDcols = cols, by = grp1]
}

aggregate_predictions_vldtn <- function(out,cols,grp){
  agg_out <- vector("list",length(grp))
  for (i in 1:length(grp)){
    agg_out[[i]] <- aggregate_predictions_by_grp(out,cols,grp[i])
  }
  return(agg_out)
}

acomb <- function(...){abind(...,along = 1)}

process_predictions <- function(out,df,n_sim){
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
  
  return(agg_out)
}

process_predictions_deaths2 <- function(out,grp,agg_df,v_e,v_ei,v_ec,r,r_i,t_sim,IFR_ratio,frlty_idx,d,w,u){
  out <- cbind(grp,out)
  
  # Merge county, age, sex and race/ethnicity information from synthetic population data frame and add risk group populations
  out <- merge(agg_df[,.(county_res,age_cat,sex,race_ethnicity,special.population,
                         asthma,diabetes,heart.disease,heart.failure,hypertension,obesity,smoker,
                         population,IFR=median_perc/100,p_clin,grp,risk_grp)],out,by="grp")
  
  # Calculate average vaccine efficacies against death and infection
  v_e_mean = v_e*(1-exp(-r*t_sim))/(r*t_sim)
  v_ei_mean = v_ei*(1-exp(-r_i*t_sim))/(r_i*t_sim)
  v_ec_mean = v_ec*(1-exp(-r_i*t_sim))/(r_i*t_sim)
  
  # Backcalculate infections in each risk group using O'Driscoll IFR
  out[,`:=`(infection=death/(IFR*IFR_ratio),
            infection_v=death_v/((1-v_e_mean)/(1-v_ei_mean)*IFR*IFR_ratio),
            infection_i=death_i/((1-v_e_mean)/(1-v_ei_mean)*IFR*IFR_ratio))]
  out[special.population %in% c(3,8),`:=`(infection=infection/frlty_idx,
                                          infection_v=infection_v/frlty_idx,
                                          infection_i=infection_i/frlty_idx)]
  
  # Calculate clinical cases from infections using Davies age-dependent clinical fraction
  out[,`:=`(case=infection*p_clin,
            case_v=(1-v_ec_mean)/(1-v_ei_mean)*infection_v*p_clin,
            case_i=(1-v_ec_mean)/(1-v_ei_mean)*infection_i*p_clin)]
  
  # Calculate and add DALYs from infections and cases
  out[,`:=`(DALYs = DALYs + calc_DALYs_from_infections_and_cases(infection+infection_i,case+case_i,d,w),
            DALYs_v = DALYs_v + calc_DALYs_from_infections_and_cases(infection_v+infection_i,case_v+case_i,d,w),
            QALYs = QALYs + calc_QALYs_from_infections_and_cases(infection+infection_i,case+case_i,u),
            QALYs_v = QALYs_v + calc_QALYs_from_infections_and_cases(infection_v+infection_i,case_v+case_i,u))]
  
  # Calculate infections, cases, deaths and DALYs averted
  out <- calc_infections_cases_deaths_DALYs_averted(out)
  
  return(out)
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
  # Calculate infections, cases, deaths, DALYs, and QALYs averted
  x[,`:=`(deaths_averted = death - death_v,
          DALYs_averted = DALYs - DALYs_v,
          QALYs_averted = QALYs - QALYs_v,
          infections_averted = infection - infection_v,
          cases_averted = case - case_v)]
  
  # Calculate proportion of infections, cases, deaths and DALYs averted
  x[,`:=`(diff_prop_death = deaths_averted/population,
          diff_DALYs_pp = DALYs_averted/population,
          diff_QALYs_pp = QALYs_averted/population,
          diff_prop_infections = infections_averted/population,
          diff_prop_cases = cases_averted/population)]
}

optimize_allocation <- function(x,n_v,by=""){  
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

  return(x)
}

order_by_risk <- function(x,by=""){
  if (by==""){
    ord <- x[order(-lambda_adj),grp]
  } else if (by=="DALYs"){
    ord <- x[order(-lambda_DALYs),grp]
  } else if (by=="QALYs"){
    ord <- x[order(-lambda_QALYs),grp]
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
    # ord <- x[order(rank),grp]
    # ord <- x[order(rank,-lambda_adj),grp]
    ord <- x[order(rank,sample.int(nrow(x))),grp]
  } else if (by=="QALYs"){
    agg_x <- x[,.(lambda_QALYs=mean(lambda_QALYs)),by=grp1]
    agg_x[,rank := order(order(-lambda_QALYs))]
    x <- merge(x,agg_x[,!"lambda_QALYs"],by=grp1,all.x=T)
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
  
  reshape_and_plot_four_week_ahead(x,NULL,fdir)
  reshape_and_plot_four_week_ahead(x,"county_res",fdir)
  reshape_and_plot_four_week_ahead(x,"age_cat",fdir)
  reshape_and_plot_four_week_ahead(x,"sex",fdir)
  reshape_and_plot_four_week_ahead(x,"race_ethnicity",fdir)
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

plot_all_strtgs_by_vrble <- function(x_long,vrble,grp,xlbl,lbls,plt){
  vrble1 <- paste0(vrble,"_v")
  last_char <- substr(vrble,nchar(vrble),nchar(vrble))
  ylbl <- ifelse(last_char=="s",vrble,paste0(vrble,"s"))
  q_long <- x_long[variable %in% paste0(rep(c(vrble,vrble1),each=3),c("","_q_95_LB","_q_95_UB")),]
  q_long[variable==vrble1,variable:=vrble]
  q_long[variable==paste0(vrble1,"_q_95_LB"),variable:=paste0(vrble,"_q_95_LB")]
  q_long[variable==paste0(vrble1,"_q_95_UB"),variable:=paste0(vrble,"_q_95_UB")]
  q <- dcast(q_long,... ~ variable,value.var = "value")
  p <- ggplot(q,aes(x=as.factor(q[[grp]]),y=q[[vrble]],fill=as.factor(strategy))) + 
    geom_bar(position="dodge",stat="identity") + ylab(ylbl) + 
    geom_errorbar(aes(ymin=q[[paste0(vrble,"_q_95_LB")]],ymax=q[[paste0(vrble,"_q_95_UB")]]),width=0.3,position=position_dodge(0.9))
  if (is.null(xlbl)){
    p <- p + xlab(grp)
  } else {
    p <- p + xlab(xlbl)
  }
  if (grp=="strategy"){
    p <- p + theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1)) + scale_x_discrete(breaks=0:max(x_long[,strategy]),labels = lbls) + scale_fill_brewer(palette = plt)
  } else {
    p <- p + scale_fill_brewer(name = "Strategy",labels = lbls,palette = plt)
    if (grp == "county_res"){
      p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1)) 
    }
  }
  return(p)
}

plot_all_strtgs <- function(x,grp,fdir,vrble=c("case","death","DALYs"),xlbl=NULL,lbls=c("no vaccination","random","special populations","age","essential workers","comorbidities"),plt="Spectral",by="DALYs"){
  if (grp=="strategy"){
    x_long <- melt(x,id.vars = "strategy")
    vrble <- ifelse(sapply(1:length(vrble),function(i){substr(vrble[i],nchar(vrble[i]),nchar(vrble[i]))})=="s",vrble,paste0(vrble,"s"))
  } else {
    x_long <- melt(x,id.vars = c("strategy",grp))
    # Remove duplicated "no vaccination rows"
    x_long <- x_long[!(duplicated(x_long[,-c("strategy","value")]) & (variable %in% paste0(rep(vrble,each=3),c("","_q_95_LB","_q_95_UB")))),]
    # Change strategy number for no vaccination to 0
    x_long[variable %in% paste0(rep(vrble,each=3),c("","_q_95_LB","_q_95_UB")),strategy:=0]
  }
  lp <- lapply(vrble,function(y) plot_all_strtgs_by_vrble(x_long,y,grp,xlbl,lbls,plt))
  for (i in 1:(length(lp)-1)){
    lp[[i]] <- lp[[i]] + theme(axis.title.x=element_blank(),axis.text.x=element_blank())
  }
  if (grp=="strategy"){
    w <- 4
    h <- 9
    rh <- 0.78
  } else {
    w <- 7.5
    h <- 12
    rh <- 0.9
  }
  p <- plot_grid(plotlist=lp,align="v",nrow=length(lp),rel_heights=c(rep(rh,length(lp)-1),1),labels="AUTO")
  ggsave(paste0(fdir,"pred_cases_deaths_",by,"_by_",grp,".pdf"),p,width = w,height = h)
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

plot_predictions_all_strtgs2 <- function(res_ss,lx_ss,fdir,vrble,lbls,plt,strategies1,grp1,lbls1,by="DALYs"){
  dir.create(fdir,recursive = T)
  plot_all_strtgs(res_ss,"strategy",fdir,vrble,"Strategy",lbls,plt,by)
  for (i in 1:length(lx_ss)){
    plot_all_strtgs(lx_ss[[i]][strategy %in% strategies1],grp1[i],fdir,vrble,lbls1[i],lbls,plt,by)
  }
}

plot_predictions_vacc_eff_SA <- function(res,strategies1,v_e_nms,vrble,lbls,lbls2){
  setDT(res)
  res_long <- melt(res[strategy %in% c(0,strategies1),],id.vars = c("strategy","v_e"))
  res_long$v_e <- factor(res_long$v_e,levels = v_e_nms[c(3,2,1)])
  
  p <- vector("list",length(vrble))
  for (i in 1:length(vrble)){
    last_char <- substr(vrble[i],nchar(vrble[i]),nchar(vrble[i]))
    vrble[i] <- ifelse(last_char=="s",vrble[i],paste0(vrble[i],"s"))
    q_long <- res_long[variable %in% paste0(rep(vrble[i],each=3),c("","_q_95_LB","_q_95_UB")),]
    q <- dcast(q_long,... ~ variable,value.var = "value")
    p[[i]] <- ggplot(q,aes(x=as.factor(strategy),y=q[[vrble[i]]],fill=as.factor(v_e))) + 
      geom_bar(position="dodge",stat="identity") + xlab("Strategy") + ylab(vrble[i]) + 
      geom_errorbar(aes(ymin=q[[paste0(vrble[i],"_q_95_LB")]],ymax=q[[paste0(vrble[i],"_q_95_UB")]]),width=0.3,position=position_dodge(0.9)) + 
      theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_x_discrete(breaks=0:max(res_long[,strategy]),labels = lbls) + 
      scale_fill_brewer(name = "Efficacy",labels = lbls2,palette = "YlOrRd")
  }
  return(p)
}

plot_vacc_distn <- function(x,v,vrble,grp,fdir,xlbl=NULL,ylbl=NULL){
  agg_v <- aggregate_by_grp(v,vrble,grp)
  agg_x <- aggregate_by_grp(x,vrble,grp)
  agg_v <- merge(agg_v,agg_x,by=grp,all.x=T)
  agg_v$prop <- agg_v$x.x/agg_v$x.y
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
    # w <- ifelse(grp=="county_res",7,5)
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
    # w <- 7
  }
  # pdf(paste0(fdir,"vacc_propn_by_",paste(grp,collapse = "_"),".pdf"),width = w, height = 4)
  # print(p)
  # dev.off()
  return(p)
}

plot_optimal_allocation <- function(x,n_v,fdir,spec_pops=c("non_spec_pop","HCW","prisoner","SNF","education","homeless")){
  idx <- (x$cum_pop < n_v)
  v <- x[idx,]
  dir.create(fdir,recursive = T)
  p1 <- plot_vacc_distn(x,v,"population","county_res",fdir,"County")
  p2 <- plot_vacc_distn(x,v,"population","age_cat",fdir,"Age")
  p3 <- plot_vacc_distn(x,v,"population","sex",fdir,"Sex")
  p4 <- plot_vacc_distn(x,v,"population","race_ethnicity",fdir,"Race/ethnicity")
  x$spec_pop <- spec_pops[x$special.population+1]
  v$spec_pop <- spec_pops[v$special.population+1]
  p5 <- plot_vacc_distn(x,v,"population","spec_pop",fdir,"Special population")
  p6 <- plot_vacc_distn(x,v,"population","comorb",fdir,"Comorbidity")
  p <- plot_grid(p2,p1,p3,p4,p6,p5,align="h",axis="b",nrow=3,ncol=2,labels="AUTO")
  ggsave(paste0(fdir,"vacc_propn_by_single_vrble.pdf"),p,width = 5*2, height = 4*3)
  
  # Two-way breakdowns
  # plot_vacc_distn(x,v,"population",c("county_res","age_cat"),fdir,"County","Age")
  p1 <- plot_vacc_distn(x,v,"population",c("age_cat","county_res"),fdir,"Age","County")
  p2 <- plot_vacc_distn(x,v,"population",c("age_cat","sex"),fdir,"Age","Sex")
  p3 <- plot_vacc_distn(x,v,"population",c("age_cat","race_ethnicity"),fdir,"Age","Race/ethnicity")
  p4 <- plot_vacc_distn(x,v,"population",c("age_cat","spec_pop"),fdir,"Age","Special population")
  p5 <- plot_vacc_distn(x,v,"population",c("age_cat","comorb"),fdir,"Age","Comorbidity")
  l = get_legend(p1)
  pl <- plot_grid(p1 + theme(legend.position = "none"),
                  p3 + theme(legend.position = "none"),
                  p4 + theme(legend.position = "none"),
                  align = "v",
                  axis = "l",
                  nrow = 3,
                  ncol = 1,
                  labels = c("A","C","E"))
  pr <- plot_grid(p2 + theme(legend.position = "none"),
                  p5 + theme(legend.position = "none"),
                  l,
                  align = "v",
                  axis = "l",
                  nrow = 3,
                  ncol = 1,
                  labels = c("B","D",""))
  p <- plot_grid(pl,pr,rel_widths = c(1,0.8))
  ggsave(paste0(fdir,"vacc_propn_by_two_vrbles.pdf"),p,width = 10, height = 4*3)
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

order_and_calc_impact <- function(x_HCW_LTCF,x_other,x_under20,strategy,agg_df_other,by,comorbs,n_v,res_nms){
  # ord <- order_by_risk(agg_df_HCW_LTCF,by)
  # x_HCW_LTCF <- x_HCW_LTCF[match(ord,x_HCW_LTCF$grp),]
  x_HCW_LTCF <- x_HCW_LTCF[sample.int(nrow(x_HCW_LTCF)),]
  if (strategy==1){ # random
    x_other <- x_other[sample.int(nrow(x_other)),]
    x <- rbind(x_HCW_LTCF,x_other,x_under20)
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
    x <- rbind(x_HCW_LTCF,x_spec_pop,x_other,x_under20)
  } else if (strategy==3){ # age
    ord <- order_by_risk_by_grp(agg_df_other,"age_cat",by)
    x_other <- x_other[match(ord,x_other$grp),]
    x <- rbind(x_HCW_LTCF,x_other,x_under20)
  } else if (strategy==4){ # essential worker
    EW_idx <- (x_other$special.population %in% c(6,7))
    x_EW <- x_other[EW_idx,]
    x_other <- x_other[!EW_idx,]
    x_EW <- x_EW[sample.int(nrow(x_EW)),]
    x_other <- x_other[sample.int(nrow(x_other)),]
    x <- rbind(x_HCW_LTCF,x_EW,x_other,x_under20) 
  } else if (strategy==5){ # comorbidity
    comorb_idx <- (rowSums2(as.matrix(x_other[,..comorbs]))!=0)
    x_comorb <- x_other[comorb_idx,]
    x_other <- x_other[!comorb_idx,]
    x_comorb <- x_comorb[sample.int(nrow(x_comorb)),]
    x_other <- x_other[sample.int(nrow(x_other)),]
    x <- rbind(x_HCW_LTCF,x_comorb,x_other,x_under20)
  } else if (strategy==6){ # age and county
    ord <- order_by_risk_by_grp(agg_df_other,c("age_cat","county_res"),by)
    x_other <- x_other[match(ord,x_other$grp),]
    x <- rbind(x_HCW_LTCF,x_other,x_under20)
  } else if (strategy==7){ # age and special population
    ord <- order_by_risk_by_grp(agg_df_other,c("age_cat","special.population"),by)
    x_other <- x_other[match(ord,x_other$grp),]
    x <- rbind(x_HCW_LTCF,x_other,x_under20)
  } else if (strategy==8){ # optimal
    ord <- order_by_risk(agg_df_other,by)
    x_other <- x_other[match(ord,x_other$grp),]
    x <- rbind(x_HCW_LTCF,x_other,x_under20)
  }
  x$population[is.na(x$population)] <- 0
  x$cum_pop <- cumsum(x$population)
  idx <- (x$cum_pop <= n_v)
  x[!idx,c("infection_v","case_v","death_v","DALYs_v","QALYs_v")] <- x[!idx,c("infection","case","death","DALYs","QALYs")]
  x$strategy <- strategy
  res <- x[idx,lapply(.SD,sum),.SDcols = res_nms]
  return(list(x=x,res=res))
}

run_prioritisation_strategies2 <- function(x,n_v,agg_df,strategies,by="",comorbs=c("asthma","diabetes","smoker","heart.disease","heart.failure","hypertension","obesity"),dir="",t_sim=180,lbls=c("no vaccination","random","special populations","age","essential workers","comorbidities"),plt="Spectral",spec_pops=c("non-spec-pop","HCW","prisoner","SNF","education","homeless","frontline EW","non-frontline EW","ALF")){
  cols <- c("infections","cases","deaths","DALYs","QALYs")
  res_nms <- paste0(cols,"_averted")
  
  HCW_LTCF_idx <- (x$special.population %in% c(1,3,8))
  under20_idx <- (x$age_cat %in% c("<10","10-19"))
  x_HCW_LTCF <- x[HCW_LTCF_idx & !under20_idx,]
  x_other <- x[!HCW_LTCF_idx & !under20_idx,]
  x_under20 <- x[under20_idx,] 
  
  agg_df_HCW_LTCF <- agg_df[grp %in% x_HCW_LTCF$grp,]
  agg_df_other <- agg_df[grp %in% x_other$grp,]

  # Strategy 0: no vaccination
  cols1 <- cols
  cols1[1:3] <- sub("(.*)s","\\1",cols[1:3])
  totals <- apply(x[,..cols1],2,sum)
  
  # Prioritization strategies
  lx <- vector("list",length(strategies))
  lres <- vector("list",length(strategies))
  for (i in 1:length(strategies)){
    l <- order_and_calc_impact(x_HCW_LTCF,x_other,x_under20,strategies[i],agg_df_other,by,comorbs,n_v,res_nms)
    lx[[i]] <- l$x
    lres[[i]] <- l$res
  }
  
  # Combine results into single data frame
  res0 <- data.frame(matrix(0,nrow=1,ncol=length(res_nms)))
  names(res0) <- res_nms
  res <- data.frame(cbind(strategy = c(0,strategies),rbind(res0,do.call(rbind,lres))))
  res$prop_infections_averted <- res$infections_averted/totals["infection"]
  res$prop_cases_averted <- res$cases_averted/totals["case"]
  res$prop_deaths_averted <- res$deaths_averted/totals["death"]
  res$prop_DALYs_averted <- res$DALYs_averted/totals["DALYs"]
  res$prop_QALYs_averted <- res$QALYs_averted/totals["QALYs"]
  res[,cols] <- t(totals-t(res[,res_nms]))

  # Combine allocations for all strategies into one data frame
  x_all <- do.call(rbind,lx)
  x_all$comorb <- as.integer(rowSums2(as.matrix(x_all[,..comorbs]))!=0)
  cols2 <- c(cols1,paste0(cols1,"_v"))
  # x_all <- x_all[,lapply(.SD, sum), .SDcols = cols2, by = .(strategy,county_res,age_cat,sex,race_ethnicity,special.population,comorb)]
  
  x_all_county <- x_all[,lapply(.SD, sum), .SDcols = cols2, by = .(strategy,county_res)]
  x_all_age <- x_all[,lapply(.SD, sum), .SDcols = cols2, by = .(strategy,age_cat)]
  x_all_sex <- x_all[,lapply(.SD, sum), .SDcols = cols2, by = .(strategy,sex)]
  x_all_race_ethnicity <- x_all[,lapply(.SD, sum), .SDcols = cols2, by = .(strategy,race_ethnicity)]
  x_all_spec_pop <- x_all[,lapply(.SD, sum), .SDcols = cols2, by = .(strategy,special.population)]
  x_all_comorb <- x_all[,lapply(.SD, sum), .SDcols = cols2, by = .(strategy,comorb)]
  
  return(list(res=res,x_all_county=x_all_county,x_all_age=x_all_age,x_all_sex=x_all_sex,x_all_race_ethnicity=x_all_race_ethnicity,x_all_spec_pop=x_all_spec_pop,x_all_comorb=x_all_comorb))
}

q_95_LB <- function(x){quantile(x,probs = 0.025)}

q_95_UB <- function(x){quantile(x,probs = 0.975)}

calc_summary_stats <- function(x,cols,grp1){
  x_mean <- x[, lapply(.SD, mean), .SDcols = cols, by = grp1]
  x_median <- x[, lapply(.SD, median), .SDcols = cols, by = grp1][, setnames(.SD, cols, paste(cols, "median", sep = "_"))]
  x_q_95_LB <- x[, lapply(.SD, q_95_LB), .SDcols = cols, by = grp1][, setnames(.SD, cols, paste(cols, "q_95_LB", sep = "_"))]
  x_q_95_UB <- x[, lapply(.SD, q_95_UB), .SDcols = cols, by = grp1][, setnames(.SD, cols, paste(cols, "q_95_UB", sep = "_"))]
  x_ss <- cbind(x_mean,x_median[,!grp1,with=F],x_q_95_LB[,!grp1,with=F],x_q_95_UB[,!grp1,with=F])
  return(x_ss)
}

calc_summary_stats2 <- function(x,cols,grp1){
  x_mean <- x[, lapply(.SD, mean), .SDcols = cols, by = grp1]
  x_median <- x[, lapply(.SD, median), .SDcols = cols, by = grp1][, setnames(.SD, cols, paste(cols, "median", sep = "_"))]
  x_q_95_LB <- x[, lapply(.SD, min), .SDcols = cols, by = grp1][, setnames(.SD, cols, paste(cols, "q_95_LB", sep = "_"))]
  x_q_95_UB <- x[, lapply(.SD, max), .SDcols = cols, by = grp1][, setnames(.SD, cols, paste(cols, "q_95_UB", sep = "_"))]
  x_ss <- cbind(x_mean,x_median[,!grp1,with=F],x_q_95_LB[,!grp1,with=F],x_q_95_UB[,!grp1,with=F])
  return(x_ss)
}

mean_and_CI <- function(x,l,u,f=1,d=1,method="round"){
  if (method=="signif"){
    paste0(signif(f*x,d)," (",signif(f*l,d),"-",signif(f*u,d),")")
  } else if (method=="round"){
    paste0(round(f*x,d)," (",round(f*l,d),"-",round(f*u,d),")")
  }
}

calc_outcome_averted <- function(p_past_inf,lambda,r,t_sim,population,v_e){
  out <- (1-p_past_inf)*(1-exp(-lambda*t_sim))*population  
  out_v <- (1-p_past_inf)*(1-exp(-lambda*(t_sim-v_e/r*(1-exp(-r*t_sim)))))*population
  out_averted <- out - out_v 
  return(out_averted)
}

num_and_perc <- function(x,p,f=1,d=1,method="round"){
  if (method=="signif"){
    paste0(signif(f*x,d)," (",round(100*p),")")
  } else if (method=="round"){
    paste0(round(f*x,d)," (",round(100*p),")")
  }
}

calc_RR_CI <- function(RR,a,b,c,d,alpha = 0.05){
  z <- qnorm(1-alpha/2)
  CI_low <- exp(log(RR) - z * sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d)))
  CI_upp <- exp(log(RR) + z * sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d)))
  list(CI_low=CI_low,CI_upp=CI_upp)
}