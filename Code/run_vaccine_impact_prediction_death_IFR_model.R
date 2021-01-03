rm(list=ls())

library(reshape2)
library(abind)
library(doParallel)
library(data.table)
library(matrixStats)
library(ggplot2)
library(viridis)
library(RColorBrewer)

source("processing_functions.R")
source("prediction_functions.R")
source("calc_hosp_ICU_death_risk.R")

# Register doParallel backend with maximum no. workers - 1
registerDoParallel(detectCores()-1)
# registerDoParallel(4)

## Inputs
dir0 <- "../Data/CalibrationDeathIFRModel/"
deaths_ratio <- readRDS(paste0(dir0,"deaths_multiplier.RDS"))

# Number of simulations
n_sim <- 1000 #500 #100

# Outcome to optimize allocation by
by <- "DALYs" # use "" for deaths

# Set vaccine parameters
v_e <- 0.9 # vaccine efficacy
n_v <- c(2e6,5e6,10e6) # number of vaccines
t_sim <- 180 #c(120,180,270) # prediction time horizon in days
r <- 1.1^(1/28) # per day multiplier for death risk = 10% increase every 4 weeks

# Median durations of illness and disability weights for mild, moderate and severe cases
d_mi <- 7 # days
w_mi <- 0.005 # Salomon Lancet 2012 Table 2
d_m <- 10 # days
w_m <- 0.053 # Salomon Lancet 2012 Table 2
d_s <- 18 # days Lewnard BMJ 2020 + Wang JAMA 2020 + Roy Soc report
w_s <- 0.210 # Salomon Lancet 2012 Table 2
d <- c(d_mi,d_m,d_s)
w <- c(w_mi,w_m,w_s)

# Vectors of labels for comorbidities, strategies, and special populations
demogrphcs <- c("county_res","age_cat","sex","race_ethnicity")
comorbs <- c("asthma","diabetes","smoker","heart.disease","heart.failure","hypertension","obesity")
lbls <- c("no vaccination","random","special populations","age","essential workers","comorbidities")
spec_pops <- c("non_spec_pop","HCW","prisoner","SNF","education","homeless","EW")

# Load synthetic population data frame with risk estimates added
df <- readRDS("../Data/CA_pop_with_risk_ests_deaths.RDS")
# # Subset to Alameda County for testing
# df1 <- df
# df <- df1[df1$county_res=="Alameda",]
# # Add ID variable
# df$id <- 1:nrow(df)

# Create data table of different risk groups
cols <- c(demogrphcs,"special.population",comorbs,"grp")
risk_grp_dt <- merge(unique(df[,..cols],by=cols),df[,.(population=.N),by=.(grp)],by="grp",all.x=T)

# Read in IFR data
IFR_long <- readRDS("../Data/IFR_by_age_ODriscoll.RDS")
IFR_long$sex <- NULL
names(IFR_long)[names(IFR_long)=="grp"] <- "sex"
setDT(IFR_long)
IFR_ratio <- readRDS("../Data/IFR_ratio.RDS") # estimated ratio of CA IFR to O'Driscoll IFR 

# Find unique risk groups
agg_df <- unique(df,by = "grp")
# Calculate RR-adjusted death risk
agg_df[,lambda_adj := lambda * RR]

# Calculate infection risk from death risk using IFR
agg_df[,lambda_inf := lambda_adj/IFR]

# Calculate clinical case risk from infection risk using age-dependent clinical fraction
p_clin <- calc_hosp_ICU_death_risk()$p_clin_given_infctn
p_clin_dt <- data.table(age_cat=names(p_clin),p_clin=p_clin)
agg_df <- merge(agg_df,p_clin_dt,by="age_cat",all.x=T)
agg_df[,lambda_case := lambda_inf*p_clin]

# Calculate mean life expectancy by age group, sex and race/ethnicity
agg_df_LE <- df[,.(mean_life_expectancy=mean(life_expectancy)),by=.(age_cat,sex,race_ethnicity)]
# Merge with risk group data table
agg_df <- merge(agg_df,agg_df_LE,by=c("age_cat","sex","race_ethnicity"),all.x=T)

# Calculate DALY risk
agg_df[,lambda_DALYs:=calc_DALYs_from_deaths(agg_df[,lambda_adj],agg_df[,mean_life_expectancy],d[3],w[3]) + calc_DALYs_from_infections_and_cases(agg_df[,lambda_inf],agg_df[,lambda_case],d,w)]

# # Split data frame into data frames for each county for next step
# ldf <- split(df,list(df$county_res,df$age_cat,df$sex,df$race_ethnicity),sep="")
# # Remove slashes from data frame names
# names(ldf) <- gsub("/","-",names(ldf))
# 

# Get risk group IDs from synthetic population data table
grp <- df[,grp]

# Delete synthetic population data table
rm(df)

## Simulate outcomes (cases, hospitalisations and deaths) under no vaccination and vaccination
dir1 <- "PredictionsDeathIFRModel/"
dir.create(paste0("../Data/",dir1))

# tstart <- Sys.time()
# 
# # for (j in 1:length(t_sim)){
# #   foreach(i=1:length(ldf)) %dopar% {
# #     out <- simulate_deaths(ldf[[i]],v_e,t_sim[j],r,deaths_ratio[1],n_sim,d[3],w[3])
# #     # Save output
# #     saveRDS(out,file = paste0("../Data/",dir1,"pred_outcomes_",names(ldf)[i],"_",100*v_e,"prcnt_eff_",t_sim[j],"days.RDS"))
# #   }
# # }
# 
# for (j in 1:length(t_sim)){
#   foreach(i=1:n_sim) %dopar% {
#     out <- simulate_deaths(df,v_e,t_sim[j],r,deaths_ratio[1],1,d[3],w[3])
#     # Save output
#     saveRDS(out,file = paste0("../Data/",dir1,"pred_outcomes_",i,"_",100*v_e,"prcnt_eff_",t_sim[j],"days.RDS"))
#   }
# }
# 
# tend <- Sys.time()
# print(tend-tstart)

## Aggregate output
tstart <- Sys.time()

# # lmean_agg_out <- vector("list",length(ldf))
# # for (i in 1:length(ldf)){
# # mean_agg_out <- vector("list",length(t_sim))
# for (j in 1:length(t_sim)){
#   x <- foreach(i=1:length(ldf),.combine = rbind) %dopar% {
#     # agg_out <- foreach(i=1:10,.combine = acomb,.multicombine = T) %dopar% {
#     out <- readRDS(paste0("../Data/",dir1,"pred_outcomes_",names(ldf)[i],"_",100*v_e,"prcnt_eff_",t_sim[j],"days.RDS"))
#     
#     # # Process simulated output
#     # # tstart2 <- Sys.time()
#     # # lmean_agg_out[[i]] <- process_predictions(out,ldf[[i]],n_sim)
#     # # aggregate_predictions(out,ldf[[i]],n_sim)
#     # # Aggregate outcomes over each risk group
#     # agg_out <- aggregate_predictions_deaths(out,ldf[[i]],n_sim)
#     process_predictions_deaths(out,ldf[[i]],n_sim,IFR_long,IFR_ratio,p_clin_dt,d,w)
#     # tend2 <- Sys.time()
#     # print(tend2-tstart2)
#   }  
#   # mean_agg_out <- do.call(rbind,lmean_agg_out)
#   
#   # Save output
#   saveRDS(x,paste0(dir1,"pred_vacc_impact_",100*v_e,"prcnt_eff_",t_sim[j],"days.RDS"))
#   # mean_agg_out[[j]] <- x
# }

# proc_res <- function(i){
#   out <- readRDS(paste0("../Data/",dir1,"pred_outcomes_",i,"_",100*v_e,"prcnt_eff_",t_sim[j],"days.RDS"))
#   
#   # Process simulated output
#   agg_out <- process_predictions_deaths(out,df,1,risk_grp_dt,IFR_long,IFR_ratio,p_clin_dt,d,w)
#   
#   # Run vaccination prioritization strategies
#   l <- run_prioritisation_strategies2(agg_out,n_v[k],fdir,agg_df,"",comorbs,dir1,t_sim[j],lbls,"YlGnBu",spec_pops)
#   
#   return(l)
# }

for (j in 1:length(t_sim)){
  for (k in 1:length(n_v)){
    # tstart <- Sys.time()
    # l <- vector("list",n_sim)
    # x <- foreach(i=1:n_sim,.combine = rbind) %dopar% {
    # l <- foreach(i=1:n_sim,.combine = 'comb',.multicombine = T,.init = list(list(),list())) %dopar% {
    # for (i in 1:n_sim){
    l <- foreach(i=1:n_sim) %dopar% {
      out <- readRDS(paste0("../Data/",dir1,"pred_outcomes_",i,"_",100*v_e,"prcnt_eff_",t_sim[j],"days.RDS"))

      # Process simulated output
      agg_out <- process_predictions_deaths(out,grp,1,risk_grp_dt,IFR_long,IFR_ratio,p_clin_dt,d,w)

      # Run vaccination prioritization strategies
      # l[[i]] <- run_prioritisation_strategies2(agg_out,n_v[k],fdir,agg_df,"",comorbs,dir1,t_sim[j],lbls,"YlGnBu",spec_pops)
      run_prioritisation_strategies2(agg_out,n_v[k],fdir,agg_df,by,comorbs,dir1,t_sim[j],lbls,"YlGnBu",spec_pops)
    }
    # l <- mclapply(1:n_sim,proc_res,mc.cores = 6)
    # tend <- Sys.time()
    # print(tend-tstart)
    
    # Row bind output lists
    res <- do.call(rbind,lapply(l,"[[",1))
    lx <- vector("list",length(l[[1]])-1)
    for (m in 1:(length(l[[1]])-1)){
      lx[[m]] <- do.call(rbind,lapply(l,"[[",m+1))
    }
    
    # Save output
    saveRDS(res,paste0("../Data/",dir1,"pred_vacc_impact_",n_v[k],"doses_",100*v_e,"prcnt_eff_",t_sim[j],"days",by,".RDS"))
    saveRDS(lx,paste0("../Data/",dir1,"pred_vacc_impact_by_covar",n_v[k],"doses_",100*v_e,"prcnt_eff_",t_sim[j],"days",by,".RDS"))
  }
}

tend <- Sys.time()
print(tend-tstart)

ord <- order_by_risk(agg_df,by)
x <- agg_df[match(ord,grp),]
x[,comorb:=as.numeric(rowSums2(as.matrix(x[,comorbs,with=F]))!=0)]
x[is.na(population),population:=0]
x[,cum_pop:=cumsum(population)]
setDF(x)

## Calculate summary statistics and plot results
for (j in 1:length(t_sim)){
  for (k in 1:length(n_v)){
    res <- readRDS(paste0("../Data/",dir1,"pred_vacc_impact_",n_v[k],"doses_",100*v_e,"prcnt_eff_",t_sim[j],"days",by,".RDS"))
    lx <- readRDS(paste0("../Data/",dir1,"pred_vacc_impact_by_covar",n_v[k],"doses_",100*v_e,"prcnt_eff_",t_sim[j],"days",by,".RDS"))

    # Calculate summary statistics across simulations
    setDT(res)
    names(res)[names(res)=="Strategy"] <- "strategy"
    cols1 <- names(res)[names(res)!="strategy"]
    res_ss <- calc_summary_stats(res,cols1,"strategy")
    write.csv(res_ss,paste0("../Data/",dir1,"pred_cases_deaths_DALYs_all_strategies_",n_v[k],"doses_",t_sim[j],"days",by,".csv"),row.names = F)

    grp1 <- c(demogrphcs,"special.population","comorb")
    lbls1 <- c("County","Age","Sex","Race/ethnicity","Special population","Comorbidity")
    lx_ss <- vector("list",length(lx))
    for (i in 1:length(lx)){
      cols2 <- names(lx[[i]])[!(names(lx[[i]]) %in% c("strategy",grp1[i]))]
      lx_ss[[i]] <- calc_summary_stats(lx[[i]],cols2,c("strategy",grp1[i]))
    }
    lx_ss[[5]][,special.population := spec_pops[special.population+1]]

    # Plot results
    fdir <- paste0("../Figures/",dir1,t_sim[j],"days/",n_v[k],"doses/")
    plot_predictions_all_strtgs2(res_ss[res_ss$strategy %in% 0:5,],lx_ss,paste0(fdir,"Strategies1to",length(lbls)-1,by,"/"),lbls,"YlGnBu",grp1,lbls1)
    
    # Plot optimal allocation
    plot_optimal_allocation(x,n_v[k],paste0(fdir,"OptimalAllocation",by,"/"),spec_pops) 
  }
}

# save(v_e,n_v,t_sim,r,dir1,file=paste0(dir1,"sim_pars.RData"))
# 
# for (j in 1:length(t_sim)){
#   # Calculate cases, hospitalisations and deaths averted
#   # x <- mean_agg_out[[j]][,lapply(.SD,sum),by=.(county_res,age_cat,sex,race_ethnicity,special.population)]
#   x <- readRDS(paste0("../Data/",dir1,"pred_vacc_impact_",100*v_e,"prcnt_eff_",t_sim[j],"days.RDS"))
#   setDF(x)
#   
#   # Set directory for saving figures into
#   fdir <- paste0("../Figures/",dir1,t_sim[j],"days/")
# 
#   # ## Calculate impact for different vaccination strategies with different numbers of initial vaccine doses
#   # res <- vector("list",length(n_v))
#   # resDALYs <- vector("list",length(n_v))
#   # for(i in 1:length(n_v)){
#   #   # Optimal allocation for minimising deaths
#   #   res[[i]] <- run_prioritisation_strategies(x,n_v[i],fdir,"",dir1,t_sim[j],"YlGnBu")
#   # 
#   #   # Optimal allocation for averting DALYs
#   #   resDALYs[[i]] <- run_prioritisation_strategies(x,n_v[i],fdir,"DALYs",dir1,t_sim[j],"YlGnBu")
#   # }
# 
#   ## Alternative "phase 1b" prioritization strategies with HCWs and SNF residents vaccinated first in "phase 1a"
#   res1 <- vector("list",length(n_v))
#   resDALYs1 <- vector("list",length(n_v))
#   for(i in 1:length(n_v)){
#     # Optimal allocation for minimising deaths
#     res1[[i]] <- run_prioritisation_strategies2(x,n_v[i],fdir,agg_df,"",comorbs,dir1,t_sim[j],lbls,"YlGnBu",spec_pops)
# 
#     # Optimal allocation for averting DALYs
#     resDALYs1[[i]] <- run_prioritisation_strategies2(x,n_v[i],fdir,agg_df,"DALYs",comorbs,dir1,t_sim[j],lbls,"YlGnBu",spec_pops)
#   }
# 
#   ## Table of % CA cases, deaths and DALYs prevented under age targeting, by age group
#   # Extract prioritization order data frame for age targeting (strategy 3)
#   x3 <- resDALYs1[[1]]$x3
#   # Remove HCWs and SNF residents
#   x3_other <- x3[!(x3$special.population %in% c(1,3)),]
#   # Aggregate over age groups
#   agg_x3_other <- aggregate(cbind(population,cases_averted,deaths_averted,DALYs_averted) ~ age_cat,x3_other,sum)
#   agg_x3_other$prop_CA_cases_averted <- agg_x3_other$cases_averted/sum(agg_x3_other$cases_averted)
#   agg_x3_other$prop_CA_deaths_averted <- agg_x3_other$deaths_averted/sum(agg_x3_other$deaths_averted)
#   agg_x3_other$prop_CA_DALYs_averted <- agg_x3_other$DALYs_averted/sum(agg_x3_other$DALYs_averted)
#   agg_x3_other <- agg_x3_other[nrow(agg_x3_other):1,]
#   write.csv(agg_x3_other,paste0("../Data/",dir1,"prop_CA_cases_deaths_DALYs_averted_by_age_",t_sim[j],"days.csv"),row.names = F)
# }

# # risk_grp_df <- df[match(unique(df$risk_grp),df$risk_grp),c("county_res","age_cat","sex","race_ethnicity","special.population","risk_grp")]
# # risk_grp_pop <- as.data.frame(table(df$risk_grp),responseName = "population")
# # risk_grp_df <- merge(risk_grp_df,risk_grp_pop,by.x="risk_grp",by.y="Var1")  
# # 
# # mean_and_q_agg_county <- process_predictions(agg_out,risk_grp_df,n_sim,"county_res")
# 
# ## Determine optimal vaccine allocation to maximise reduction in deaths
# x <- optimize_allocation(mean_agg_out,n_v)
# print(sum(x$cases_averted,na.rm = T))
# print(sum(x$hosps_averted,na.rm = T))
# print(sum(x$deaths_averted,na.rm = T))
# print(sum(x$case,na.rm = T))
# print(sum(x$hosp,na.rm = T))
# print(sum(x$death,na.rm = T))
# 
# # Save output
# save(mean_agg_out,x,file = paste0("pred_vacc_impact_",n_v,"doses_",100*v_e,"prcnt_eff",".RData"))
# 
# # Plot predictions
# fdir <- "~/Dropbox/COVIDVaccineModelling/Figures/"
# plot_predictions(x,n_v,fdir)
# 
# # Plot optimal vaccine allocation
# plot_optimal_allocation(x,n_v,fdir)
# 
