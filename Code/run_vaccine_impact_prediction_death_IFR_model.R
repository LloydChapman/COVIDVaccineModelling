rm(list=ls())

library(reshape2)
library(abind)
library(doParallel)
library(data.table)
library(truncnorm)
library(matrixStats)
library(ggplot2)
library(viridis)
library(RColorBrewer)

source("processing_functions.R")
source("prediction_functions.R")
source("calc_hosp_ICU_death_risk.R")

# Register doParallel backend with maximum no. workers - 1
# registerDoParallel(detectCores()-1)
registerDoParallel(4)

# Load synthetic population data frame with risk estimates added
df <- readRDS("../Data/CA_pop_with_risk_ests_deaths4.RDS")
# # Subset to Alameda County for testing
# df1 <- df
# df <- df1[df1$county_res=="Alameda",]

# Multiplier for death rate from calibration to CalCAT forecasts
deaths_ratio <- 1

# Number of simulations
n_sim <- 1000

# Outcome to optimize allocation by
by <- "DALYs" # use "" for deaths

# Vaccine parameters
# vaccine efficacy
v_e <- vector("list",3)
v_e[1:2] <- c(0.6,0.95)
n_v <- c(2e6,5e6,10e6) # number of vaccines
t_sim <- round(365/2) # prediction time horizon in days
r <- 1 # per day multiplier for death rate = no change

# Prioritization strategies
# 1: random allocation
# 2: special population targeting
# 3: age targeting
# 4: essential worker targeting
# 5: comorbidity targeting
# 6: age and county targeting
# 7: age and special population targeting
# 8: optimal allocation
strategies <- 1:8

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
spec_pops <- c("non-special-population","HCW","prisoner","SNF","education","homeless","frontline essential worker","non-frontline essential worker","ALF")

# Frailty index for LTCF residents
frlty_idx <- 3

# Read in IFR data
IFR_long <- readRDS("../Data/IFR_by_age_ODriscoll.RDS")
IFR_long$sex <- NULL
names(IFR_long)[names(IFR_long)=="grp"] <- "sex"
setDT(IFR_long)
IFR_ratio <- readRDS("../Data/IFR_ratio2.RDS") # estimated ratio of CA IFR to O'Driscoll IFR 

# Aggregate simulated population by risk group for use in simulations
cols <- c(demogrphcs,"special.population",comorbs,"grp")
agg_df <- df[,lapply(.SD,mean),.SDcols = names(df)[!(names(df) %in% c(cols,"age_cat_sero","id","alive.dead_stats","race","ethnicity","time_death"))],by=cols]

# Delete synthetic population data table
rm(df)
# Garbage collect
gc(T)

# Age-dependent vaccine efficacy
v_e_age <- rep(0.6,nrow(agg_df))
v_e_age[agg_df$age_cat=="60-69"] <- 0.5
v_e_age[agg_df$age_cat=="70-79"] <- 0.4
v_e_age[agg_df$age_cat=="80+"] <- 0.3
v_e[[3]] <- v_e_age
v_e_nms <- c(paste0(as.character(100*unlist(v_e[1:2])),"prcnt_eff"),paste0(100*v_e[[3]][agg_df$age_cat=="<10"][1],"prcnt_age_dep_eff"))

# Calculate RR-adjusted death risk
agg_df[,lambda_adj := exp(log_lambda) * RR * deaths_ratio[1]]

# Calculate infection risk from death risk using CA-corrected IFR
agg_df[,lambda_inf := lambda_adj/(median_perc/100*IFR_ratio)]

# Calculate clinical case risk from infection risk using age-dependent clinical fraction
p_clin <- calc_hosp_ICU_death_risk()$p_clin_given_infctn
p_clin_dt <- data.table(age_cat=names(p_clin),p_clin=p_clin)
agg_df <- merge(agg_df,p_clin_dt,by="age_cat",all.x=T)
agg_df[,lambda_case := lambda_inf*p_clin]

# Calculate DALY risk
agg_df[,lambda_DALYs:=calc_DALYs_from_deaths(agg_df[,lambda_adj],agg_df[,life_expectancy],d[3],w[3]) + calc_DALYs_from_infections_and_cases(agg_df[,lambda_inf],agg_df[,lambda_case],d,w)]

# Get risk group IDs from synthetic population data table
grp <- agg_df[,grp]

## Simulate outcomes (cases, hospitalisations and deaths) under no vaccination and vaccination
dir1 <- "PredictionsDeathIFRModel4/"
dir.create(paste0("../Data/",dir1))


tstart <- Sys.time()

for (i in 1:length(v_e)){
  for (j in 1:length(t_sim)){
    foreach(k=1:n_sim) %dopar% {
      out <- simulate_deaths2(agg_df,v_e[[i]],t_sim[j],r,deaths_ratio[1],d[3],w[3])
      # Save output
      saveRDS(out,file = paste0("../Data/",dir1,"pred_outcomes_",k,"_",v_e_nms[i],"_",t_sim[j],"days.RDS"))
    }
  }
}

tend <- Sys.time()
print(tend-tstart)

## Aggregate output
tstart <- Sys.time()

for (i in 1:length(v_e)){
  for (j in 1:length(t_sim)){
    for (l in 1:length(n_v)){
      lres <- foreach(k=1:n_sim) %dopar% {
        out <- readRDS(paste0("../Data/",dir1,"pred_outcomes_",k,"_",v_e_nms[i],"_",t_sim[j],"days.RDS"))

        # Process simulated output
        agg_out <- process_predictions_deaths2(out,grp,agg_df,IFR_ratio,d,w,frlty_idx)

        # Run vaccination prioritization strategies
        run_prioritisation_strategies2(agg_out,n_v[l],agg_df,strategies,by,comorbs,dir1,t_sim[j],lbls,"YlGnBu",spec_pops)
      }

      # Row bind output lists
      res <- do.call(rbind,lapply(lres,"[[",1))
      lx <- vector("list",length(lres[[1]])-1)
      for (m in 1:(length(lres[[1]])-1)){
        lx[[m]] <- do.call(rbind,lapply(lres,"[[",m+1))
      }

      # Save output
      saveRDS(res,paste0("../Data/",dir1,"pred_vacc_impact_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by,".RDS"))
      saveRDS(lx,paste0("../Data/",dir1,"pred_vacc_impact_by_covar_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by,".RDS"))
    }
  }
}

tend <- Sys.time()
print(tend-tstart)

## Calculate summary statistics and plot results
# Reorder population according to a priori optimal allocation order based on risk
x <- agg_df
HCW_LTCF_idx <- (x$special.population %in% c(1,3,8))
under20_idx <- (x$age_cat %in% c("<10","10-19"))
x_HCW_LTCF <- x[HCW_LTCF_idx & !under20_idx,]
x_other <- x[!HCW_LTCF_idx & !under20_idx,]
x_under20 <- x[under20_idx,]

ord <- order_by_risk(x_other,by)
x_other <- x_other[match(ord,x_other$grp),]
x <- rbind(x_HCW_LTCF,x_other,x_under20)

x[,comorb:=as.integer(rowSums2(as.matrix(x[,comorbs,with=F]))!=0)]
x[is.na(population),population:=0]
x[,cum_pop:=cumsum(population)]

setDF(x)

# Strategies and risk factors to plot performance for
strategies1 <- strategies[1:(length(strategies)-3)]
grp1 <- c(demogrphcs,"special.population","comorb")
lbls1 <- c("County","Age","Sex","Race/ethnicity","Special population","Comorbidity")
vrble <- c("case","death","DALYs")
for (i in 1:length(v_e)){
  for (j in 1:length(t_sim)){
    for (l in 1:length(n_v)){
      res <- readRDS(paste0("../Data/",dir1,"pred_vacc_impact_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by,".RDS"))
      lx <- readRDS(paste0("../Data/",dir1,"pred_vacc_impact_by_covar_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by,".RDS"))

      # Calculate summary statistics across simulations
      setDT(res)
      # names(res)[names(res)=="Strategy"] <- "strategy"
      cols1 <- names(res)[names(res)!="strategy"]
      res_ss <- calc_summary_stats(res,cols1,"strategy")
      write.csv(res_ss,paste0("../Data/",dir1,"pred_cases_deaths_DALYs_all_strategies_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by,".csv"),row.names = F)

      lx_ss <- vector("list",length(lx))
      for (m in 1:length(lx)){
        cols2 <- names(lx[[m]])[!(names(lx[[m]]) %in% c("strategy",grp1[m]))]
        lx_ss[[m]] <- calc_summary_stats(lx[[m]],cols2,c("strategy",grp1[m]))
      }
      lx_ss[[5]][,special.population := spec_pops[special.population+1]]

      # Plot results
      fdir <- paste0("../Figures/",dir1,t_sim[j],"days/",v_e_nms[i],"/",n_v[l],"doses/")
      plot_predictions_all_strtgs2(res_ss[res_ss$strategy %in% c(0,strategies1),],lx_ss,paste0(fdir,"Strategies1to",length(lbls)-1,by,"/"),vrble,lbls,"YlGnBu",strategies1,grp1,lbls1)

      # Plot optimal allocation
      # plot_optimal_allocation(x,n_v[l],paste0(fdir,"OptimalAllocation",by,"/"),spec_pops)
      plot_optimal_allocation(x[!(x$special.population %in% c(1,3,8)),],n_v[l],paste0(fdir,"OptimalAllocation",by,"/"),spec_pops)
    }
  }
}

# Plot results of vaccine sensitivity analysis
lbls2 <- c("90%","60%","60% age-dep")
for (j in 1:length(t_sim)){
  for (l in 1:length(n_v)){
    # Read in results for different vaccine efficacies
    lresSA <- vector("list",length(v_e_nms))
    for (i in 1:length(v_e_nms)){
      lresSA[[i]] <- read.csv(paste0("../Data/",dir1,"pred_cases_deaths_DALYs_all_strategies_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by,".csv"),stringsAsFactors = F)
      lresSA[[i]]$v_e <- v_e_nms[i]
    }
    # Bind results for different efficacies together
    res <- do.call(rbind,lresSA)

    # Plot results
    fdir <- paste0("../Figures/",dir1,t_sim[j],"days/vacc_eff_SA/",n_v[l],"doses/")
    plot_predictions_vacc_eff_SA(res,strategies1,v_e_nms,vrble,fdir,lbls,lbls2)
  }
}

for (j in 1:length(t_sim)){
  x1 <- rbind(x_other,x_under20)
  setDT(x1)
  x1[,`:=`(deaths_averted=calc_outcome_averted(p_past_inf,lambda_adj,r,t_sim[j],population,v_e[[2]]),infections_averted=calc_outcome_averted(p_past_inf,lambda_inf,r,t_sim[j],population,v_e[[2]]),cases_averted=calc_outcome_averted(p_past_inf,lambda_case,r,t_sim[j],population,v_e[[2]]),DALYs_averted=calc_outcome_averted(p_past_inf,lambda_DALYs,r,t_sim[j],population,v_e[[2]]))]
  agg_x1 <- x1[,lapply(.SD,sum),.SDcols=c("population","deaths_averted","infections_averted","cases_averted","DALYs_averted"),by=.(age_cat)]
  agg_x1_long <- melt(agg_x1,id.vars = "age_cat")
  agg_x1_long[,prop_CA:=value/sum(value),by=.(variable)]
  agg_x1 <- dcast(agg_x1_long, age_cat ~ variable,value.var = c("value","prop_CA"))
  names(agg_x1) <- gsub("value_","",names(agg_x1))
  agg_x1 <- agg_x1[c((nrow(agg_x1)-1):1,nrow(agg_x1)),]
  write.csv(agg_x1,paste0("../Data/",dir1,"prop_CA_cases_deaths_DALYs_averted_by_age_",t_sim[j],"days.csv"),row.names = F)
}
