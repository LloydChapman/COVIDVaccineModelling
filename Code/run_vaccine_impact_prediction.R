rm(list=ls())

library(reshape2)
library(abind)
library(doParallel)
library(ggplot2)
library(viridis)
library(RColorBrewer)

source("~/Dropbox/COVIDVaccineModelling/Code/processing_functions.R")
source("~/Dropbox/COVIDVaccineModelling/Code/prediction_functions.R")
source("~/Dropbox/COVIDVaccineModelling/Code/calc_hosp_ICU_death_risk.R")

# Set working directory
dir <- "~/Dropbox/COVIDVaccineModelling/Data/"
setwd(dir)

# Register doParallel backend with maximum no. workers - 1
registerDoParallel(detectCores()-1)

## Inputs
dir0 <- "Calibration2/"
cases_ratio <- readRDS(paste0(dir0,"cases_multiplier1.RDS"))
deaths_ratio <- readRDS(paste0(dir0,"deaths_multiplier1.RDS"))

# Relative risks for HCWs, prisoners, SNF residents, teachers and homeless individuals
# RR <- c(HCW = 3,prisoner = 5.5,SNF = 9,teacher = 1.8,homeless = 1.65)
RR <- c(HCW = 3,prisoner = 5.5,SNF = 9,teacher = 1.8,homeless = 1.65,essential = 1.7)

# Number of simulations
n_sim <- 100

# Set vaccine parameters
v_e <- 0.9 # vaccine efficacy
n_v <- c(2e6,5e6,10e6) # number of vaccines
t_sim <- c(120,180,270) #180 # prediction time horizon in days
# r <- 1.1^(1/t_sim) # per day multiplier for force of infection = 10% increase over duration of simulation
r <- 1.1^(1/28) # per day multiplier for force of infection = 10% increase every 4 weeks

# Proportions, median durations of illness and disability weights for mild, moderate and severe cases
p_mi <- 0.4
d_mi <- 7 # days
w_mi <- 0.005 # Salomon Lancet 2012 Table 2
d_m <- 10 # days
w_m <- 0.053 # Salomon Lancet 2012 Table 2
d_s <- 18 # days Lewnard BMJ 2020 + Wang JAMA 2020 + Roy Soc report
w_s <- 0.210 # Salomon Lancet 2012 Table 2
d <- c(d_mi,d_m,d_s)
w <- c(w_mi,w_m,w_s)

# Load synthetic population data frame with risk estimates added
df <- readRDS("CA_pop_with_risk_ests1.RDS")
# # Subset to Alameda County for testing
# df1 <- df
# df <- df1[df1$county_res=="Alameda",]

# Split data frame into dataframes for each county for next step
ldf <- split(df,list(df$county_res,df$age_cat,df$sex),sep="")
rm(df)

## Simulate outcomes (cases, hospitalisations and deaths) under no vaccination and vaccination
# df <- simulate_cases_hosps_deaths(df,v_e,t_sim,r,cases_ratio[3])

# dir1 <- "Predictions/"
# dir1 <- "Predictions1/"
dir1 <- "Predictions2/"
# dir1 <- "Predictions2Corrected/"
dir.create(dir1)

tstart <- Sys.time()
for (j in 1:length(t_sim)){
  foreach(i=1:length(ldf)) %dopar% {
    out <- simulate_cases_hosps_deaths(ldf[[i]],v_e,t_sim[j],r,cases_ratio[3],deaths_ratio[3],n_sim,p_mi,d,w)
    # Save output
    saveRDS(out,file = paste0(dir1,"pred_outcomes_",names(ldf)[i],"_",100*v_e,"prcnt_eff_",t_sim[j],"days.RDS"))
  }
}
tend <- Sys.time()
print(tend-tstart)

## Aggregate output
tstart <- Sys.time()
# lmean_agg_out <- vector("list",length(ldf))
# for (i in 1:length(ldf)){
# mean_agg_out <- vector("list",length(t_sim))
for (j in 1:length(t_sim)){
  x <- foreach(i=1:length(ldf),.combine = rbind) %dopar% {
    # agg_out <- foreach(i=1:10,.combine = acomb,.multicombine = T) %dopar% {
    out <- readRDS(paste0(dir1,"pred_outcomes_",names(ldf)[i],"_",100*v_e,"prcnt_eff_",t_sim[j],"days.RDS"))
    
    # Process simulated output
    # tstart2 <- Sys.time()
    # lmean_agg_out[[i]] <- process_predictions(out,ldf[[i]],n_sim)
    # aggregate_predictions(out,ldf[[i]],n_sim)
    process_predictions(out,ldf[[i]],n_sim)
    # tend2 <- Sys.time()
    # print(tend2-tstart2)
  }  
# mean_agg_out <- do.call(rbind,lmean_agg_out)

# Save output
saveRDS(x,paste0(dir1,"pred_vacc_impact_",100*v_e,"prcnt_eff_",t_sim[j],"days.RDS"))
# mean_agg_out[[j]] <- x
}
tend <- Sys.time()
print(tend-tstart)

save(v_e,n_v,t_sim,r,dir1,file=paste0(dir1,"sim_pars.RData"))

for (j in 1:length(t_sim)){
  # Calculate cases, hospitalisations and deaths averted
  # x <- calc_cases_deaths_DALYs_averted(mean_agg_out[[j]])
  x <- readRDS(paste0(dir1,"pred_vacc_impact_",100*v_e,"prcnt_eff_",t_sim[j],"days.RDS"))
  # x <- readRDS(paste0("Predictions2/pred_vacc_impact_",100*v_e,"prcnt_eff_",t_sim[j],"days.RDS"))
  x <- calc_cases_deaths_DALYs_averted(x)
  
  # Set directory for saving figures into
  fdir <- paste0("~/Dropbox/COVIDVaccineModelling/Figures/",dir1,t_sim[j],"days/")
  
  ## Calculate impact for different vaccination strategies with different numbers of initial vaccine doses
  res <- vector("list",length(n_v))
  resDALYs <- vector("list",length(n_v))
  for(i in 1:length(n_v)){
    # Optimal allocation for minimising deaths
    res[[i]] <- run_prioritisation_strategies(x,n_v[i],fdir,"",dir1,t_sim[j],"YlGnBu")
    
    # Optimal allocation for averting DALYs
    resDALYs[[i]] <- run_prioritisation_strategies(x,n_v[i],fdir,"DALYs",dir1,t_sim[j],"YlGnBu")
  }
  
  ## Alternative "phase 1b" prioritization strategies with HCWs and SNF residents vaccinated first in "phase 1a"
  res1 <- vector("list",length(n_v))
  resDALYs1 <- vector("list",length(n_v))
  for(i in 1:length(n_v)){
    # Optimal allocation for minimising deaths
    res1[[i]] <- run_prioritisation_strategies1(x,n_v[i],fdir,"",dir1,t_sim[j],"YlGnBu")
    
    # Optimal allocation for averting DALYs
    resDALYs1[[i]] <- run_prioritisation_strategies1(x,n_v[i],fdir,"DALYs",dir1,t_sim[j],"YlGnBu")
  }  
  
  ## Table of % CA cases, deaths and DALYs prevented under age targeting, by age group
  # Extract prioritization order data frame for age targeting (strategy 3)
  x3 <- resDALYs1[[1]]$x3
  # Remove HCWs and SNF residents
  x3_other <- x3[!(x3$special.population %in% c(1,3)),]
  # Aggregate over age groups
  agg_x3_other <- aggregate(cbind(population,cases_averted,deaths_averted,DALYs_averted) ~ age_cat,x3_other,sum)
  agg_x3_other$prop_CA_cases_averted <- agg_x3_other$cases_averted/sum(agg_x3_other$cases_averted)
  agg_x3_other$prop_CA_deaths_averted <- agg_x3_other$deaths_averted/sum(agg_x3_other$deaths_averted)
  agg_x3_other$prop_CA_DALYs_averted <- agg_x3_other$DALYs_averted/sum(agg_x3_other$DALYs_averted)
  agg_x3_other <- agg_x3_other[nrow(agg_x3_other):1,]
  write.csv(agg_x3_other,paste0(dir1,"prop_CA_cases_deaths_DALYs_averted_by_age_",t_sim[j],"days.csv"),row.names = F)
}

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
