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

## Inputs
agg_deaths_test <- read.csv("../Data/agg_deaths_2020-09-30_2020-12-30_3.csv")
agg_deaths_test$sex <- ifelse(agg_deaths_test$sex==1,"male","female")

# Load synthetic population data frame with risk estimates added
df <- readRDS("../Data/CA_pop_with_risk_ests_deaths_vldtn2.RDS")
# # Subset to Alameda County for testing
# df1 <- df
# df <- df1[df1$county_res=="Alameda",]

# Multiplier for death rate from calibration to CalCAT forecasts
deaths_ratio <- 1

# Number of simulations
n_sim <- 1000

# Simulation parameters
vldtn_end_date <- as.Date("2020-09-30")
v_e <- 0 # vaccine efficacy (dummy value)
t_sim <- as.numeric(as.Date("2020-12-30") - vldtn_end_date) # prediction time horizon in days
r <- 1 # per day multiplier for death rate = no change in death rate

# Median durations of illness and disability weights for mild, moderate and severe cases
d_mi <- 7 # days
w_mi <- 0.005 # Salomon Lancet 2012 Table 2
d_m <- 10 # days
w_m <- 0.053 # Salomon Lancet 2012 Table 2
d_s <- 18 # days Lewnard BMJ 2020 + Wang JAMA 2020 + Roy Soc report
w_s <- 0.210 # Salomon Lancet 2012 Table 2
d <- c(d_mi,d_m,d_s)
w <- c(w_mi,w_m,w_s)

# Vectors of labels for demographics, comorbidities, and special populations
demogrphcs <- c("county_res","age_cat","sex","race_ethnicity")
comorbs <- c("asthma","diabetes","smoker","heart.disease","heart.failure","hypertension","obesity")

# Frailty index for LTCF residents
frlty_idx <- 3

# Read in IFR ratio
IFR_ratio <- readRDS("../Data/IFR_ratio2.RDS") # estimated ratio of CA IFR to O'Driscoll IFR 

# Aggregate simulated population by risk group for use in simulations
cols <- c(demogrphcs,"special.population",comorbs,"grp")
agg_df <- df[,lapply(.SD,mean),.SDcols = names(df)[!(names(df) %in% c(cols,"age_cat_sero","id","alive.dead_stats","race","ethnicity","time_death"))],by=cols]

# Delete synthetic population data table
rm(df)
# Garbage collect
gc(T)

# Calculate clinical case risk from infection risk using age-dependent clinical fraction
p_clin <- calc_hosp_ICU_death_risk()$p_clin_given_infctn
p_clin_dt <- data.table(age_cat=names(p_clin),p_clin=p_clin)
agg_df <- merge(agg_df,p_clin_dt,by="age_cat",all.x=T)

# Get risk group IDs
grp <- agg_df[,grp]

## Simulate outcomes (cases, hospitalisations and deaths) under no vaccination and vaccination
dir1 <- "PredictionsDeathIFRModelVldtn2/"
dir.create(paste0("../Data/",dir1))

tstart <- Sys.time()

foreach(k=1:n_sim) %dopar% {
  out <- simulate_deaths2(agg_df,v_e,t_sim,r,deaths_ratio[1],d[3],w[3])
  # Save output
  saveRDS(out,file = paste0("../Data/",dir1,"pred_outcomes_",k,".RDS"))
}

tend <- Sys.time()
print(tend-tstart)

## Aggregate output
tstart <- Sys.time()

cols1 <- c("infection","case","death")
lres <- foreach(k=1:n_sim) %dopar% {
  out <- readRDS(paste0("../Data/",dir1,"pred_outcomes_",k,".RDS"))

  # Process simulated output
  agg_out <- process_predictions_deaths2(out,grp,agg_df,IFR_ratio,d,w,frlty_idx)

  # Aggregate outcomes over risk factors
  aggregate_predictions_vldtn(agg_out,cols1,demogrphcs)
}

# Row bind output lists
lx <- vector("list",length(lres[[1]]))
for (m in 1:(length(lres[[1]]))){
  lx[[m]] <- do.call(rbind,lapply(lres,"[[",m))
}

# Save output
saveRDS(lx,paste0("../Data/",dir1,"pred_outcomes_by_covar.RDS"))

tend <- Sys.time()
print(tend-tstart)

## Calculate summary statistics and plot results
lbls <- c("County","Age","Sex","Race/ethnicity")

lx <- readRDS(paste0("../Data/",dir1,"pred_outcomes_by_covar.RDS"))

# Calculate summary statistics across simulations
lx_ss <- vector("list",length(lx))
for (m in 1:length(lx)){
  cols2 <- names(lx[[m]])[!(names(lx[[m]]) %in% demogrphcs[m])]
  lx_ss[[m]] <- calc_summary_stats2(lx[[m]],cols2,demogrphcs[m])
  obsvd_deaths <- aggregate(agg_deaths_test[,"n_deaths"],by=list(agg_deaths_test[,demogrphcs[m]]),FUN=sum)
  names(obsvd_deaths) <- c(demogrphcs[m],"n_deaths")
  lx_ss[[m]] <- merge(lx_ss[[m]],obsvd_deaths,by=demogrphcs[m])
}

# Plot predicted vs observed deaths
fdir <- paste0("../Figures/",dir1)
dir.create(fdir)

for (m in 1:length(lx)){
  lx_ss_long <- melt(lx_ss[[m]])
  lx_ss_long <- lx_ss_long[variable %in% c("death","n_deaths")]
  lx_ss_long$variable <- factor(lx_ss_long$variable,levels = c("n_deaths","death"))
  w <- ifelse(demogrphcs[m]=="county_res",7,5)
  pdf(paste0(fdir,"pred_vs_obsvd_deaths_by_",demogrphcs[m],".pdf"),width = w,height = 4)
  print(ggplot() + geom_bar(aes(x=lx_ss_long[[demogrphcs[m]]],y=value,fill=variable),data = lx_ss_long,stat = "identity",position = "dodge") + geom_errorbar(aes(x=lx_ss[[m]][[demogrphcs[m]]],ymin=lx_ss[[m]][["death_q_95_LB"]],ymax=lx_ss[[m]][["death_q_95_UB"]]),data=lx_ss[[m]],width=0.3,position=position_nudge(0.22)) + xlab(lbls[m]) + ylab("deaths") + scale_fill_discrete(name="",labels = c("Observed","Predicted")) + theme(axis.text.x = element_text(angle = 45,hjust = 1)))
  dev.off()
}
