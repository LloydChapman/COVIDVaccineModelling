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
library(patchwork)

source("processing_functions.R")
source("prediction_functions.R")
source("calc_hosp_ICU_death_risk.R")

# Register doParallel backend with maximum no. workers - 1
# registerDoParallel(detectCores()-1)
registerDoParallel(2)

# Set plot theme
theme_set(theme_cowplot(font_size = 11) + theme(
  strip.background = element_blank(),
  plot.background = element_rect(fill="white"),
  legend.background = element_rect(fill="white"),
  panel.background = element_rect(fill="white")))

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
v_e <- 0 # vaccine efficacy against death (dummy value)
t_sim <- as.numeric(as.Date("2020-12-30") - vldtn_end_date) # prediction time horizon in days
r <- -log(1-(98.2-90.4)/98.2)/182  # rate of decline in protection against death from vaccination/previous infection per day

v_ei <- 0 # vaccine efficacy against infection (dummy value)
r_i <- -log(1-(92.4-69.7)/92.4)/182 # rate of decline in protection against infection from vaccination/previous infection per day
v_ec <- 0 # vaccine efficacy against clinical disease (dummy value)

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

# Age-specific quality of life weights
QOLW <- read.csv("../Data/Quality_of_life_weights.csv")
Q <- QOLW$Q

# QALY loss from clinical infection
u <- 0.007

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

age_cat <- agg_df$age_cat
age_cat[age_cat=="<10"] <- "0-9"
age_cat[age_cat=="80+"] <- "80-89"
age_lbs <- as.numeric(gsub("-.*|\\+|>","",age_cat))
age_ubs <- as.numeric(gsub(".*-|<","",age_cat)) + 1
age_mid <- (age_lbs + age_ubs)/2
Q_a <- c(rep(Q[1:(length(Q)-1)],each=10),rep(Q[length(Q)],max(age_mid)-max(age_lbs)+agg_df[age_cat=="80+",round(max(life_expectancy))]+1))

# Get risk group IDs
grp <- agg_df[,grp]

## Simulate outcomes (cases, hospitalisations and deaths) under no vaccination and vaccination
dir1 <- "PredictionsDeathIFRModelVldtn3/"
dir.create(paste0("../Data/",dir1))

tstart <- Sys.time()

foreach(k=1:n_sim) %dopar% {
  out <- simulate_deaths2(agg_df,v_e,t_sim,r,deaths_ratio[1],d[3],w[3],Q_a)
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
  agg_out <- process_predictions_deaths2(out,grp,agg_df,v_e,v_ei,v_ec,r,r_i,t_sim,IFR_ratio,frlty_idx,d,w,u)

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

plot_prev_vs_obs_deaths_by_demogrphc <- function(x_ss,demogrphc,lbl){
  x_ss_long <- melt(x_ss,id.vars=demogrphc)
  x_ss_long <- x_ss_long[variable %in% c("death","n_deaths")]
  x_ss_long$variable <- factor(x_ss_long$variable,levels = c("n_deaths","death"))
  p <- ggplot() + 
    geom_bar(aes(x=x_ss_long[[demogrphc]],y=value,fill=variable),data = x_ss_long,stat = "identity",position = "dodge") + 
    geom_errorbar(aes(x=x_ss[[demogrphc]],ymin=x_ss[["death_q_95_LB"]],ymax=x_ss[["death_q_95_UB"]]),data=x_ss,width=0.3,position=position_nudge(0.22)) + 
    xlab(lbl) + ylab("deaths") + scale_fill_discrete(name="",labels = c("Observed","Predicted")) + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  return(p)
}

lp <- mapply(plot_prev_vs_obs_deaths_by_demogrphc,lx_ss,demogrphcs,lbls,SIMPLIFY = F)

p <- (lp[[2]] + theme(legend.position="none") +
        lp[[3]] + theme(legend.position="none") + 
        lp[[4]] + theme(legend.position="none")) / 
  lp[[1]] + plot_layout(heights=c(0.8,1)) + plot_annotation(tag_levels="A")
ggsave(paste0(fdir,"pred_vs_obsvd_deaths_by_demogrphcs.pdf"),p,width=10,height=8)
