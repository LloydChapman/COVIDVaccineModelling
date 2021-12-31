rm(list=ls())

library(reshape2)
library(abind)
library(doParallel)
library(data.table)
library(truncnorm)
library(matrixStats)
library(ggplot2)
library(cowplot)
library(viridis)
library(RColorBrewer)

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
by <- c("DALYs","QALYs") # use "" for deaths

# Vaccine parameters
# base case initial two-dose vaccine efficacies against...
v_e0 = 0.95 # death
v_ei0 = 0.85 # infection
v_ec0 = 0.9 # clinical infection

# two-dose vaccine efficacy against death
v_e <- vector("list",4) # list for storing different vaccine efficacy profiles
v_e[1:3] <- c(0.2,0.6,v_e0)
n_v <- c(2e6,5e6,10e6) # number of vaccines
t_sim <- round(365/2) # prediction time horizon in days
r <- -log(1-(98.2-90.4)/98.2)/182  # rate of decline in protection against death from vaccination/previous infection per day

# two-dose vaccine efficacy against infection
v_ei <- vector("list",length(v_e))
v_ei[1:3] <- lapply(v_e[1:3],function(x) v_ei0/v_e0*x) # assume constant scaling of efficacy against infection vs efficacy against death as efficacy against death varies
r_i <- -log(1-(92.4-69.7)/92.4)/182 # rate of decline in protection against infection from vaccination/previous infection per day

# two-dose vaccine efficacy against clinical infection
v_ec <- vector("list",length(v_e))
v_ec[1:3] <- lapply(v_e[1:3],function(x) v_ec0/v_e0*x)

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

# Age-dependent vaccine efficacy against death
v_e_age <- rep(0.6,nrow(agg_df))
v_e_age[agg_df$age_cat=="60-69"] <- 0.5
v_e_age[agg_df$age_cat=="70-79"] <- 0.4
v_e_age[agg_df$age_cat=="80+"] <- 0.3
v_e[[4]] <- v_e_age
v_e_nms <- c(paste0(as.character(100*unlist(v_e[1:3])),"prcnt_eff"),paste0(100*v_e[[4]][agg_df$age_cat=="<10"][1],"prcnt_age_dep_eff"))

# Age-dependent vaccine efficacy against infection
v_ei[[4]] <- v_ei0/v_e0 * v_e[[4]]

# Age-dependent vaccine efficacy against clinical infection
v_ec[[4]] <- v_ec0/v_e0 * v_e[[4]]

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
agg_df[,lambda_DALYs:=calc_DALYs_from_deaths(agg_df[,lambda_adj],agg_df[,life_expectancy],d[3],w[3]) + 
         calc_DALYs_from_infections_and_cases(agg_df[,lambda_inf],agg_df[,lambda_case],d,w)]

# Calculate QALY risk
age_cat <- agg_df$age_cat
age_cat[age_cat=="<10"] <- "0-9"
age_cat[age_cat=="80+"] <- "80-89"
age_lbs <- as.numeric(gsub("-.*|\\+|>","",age_cat))
age_ubs <- as.numeric(gsub(".*-|<","",age_cat)) + 1
age_mid <- (age_lbs + age_ubs)/2
Q_a <- c(rep(Q[1:(length(Q)-1)],each=10),rep(Q[length(Q)],max(age_mid)-max(age_lbs)+agg_df[age_cat=="80+",round(max(life_expectancy))]+1))
agg_df[,lambda_QALYs:=calc_QALYs_from_deaths(agg_df[,lambda_adj],agg_df[,life_expectancy],age_mid,Q_a) +
         calc_QALYs_from_infections_and_cases(agg_df[,lambda_inf],agg_df[,lambda_case],u)]

# Get risk group IDs from synthetic population data table
grp <- agg_df[,grp]

## Simulate outcomes (cases, hospitalisations and deaths) under no vaccination and vaccination
dir1 <- "PredictionsDeathIFRModel5/"
dir.create(paste0("../Data/",dir1))


tstart <- Sys.time()

for (i in 1:length(v_e)){
  for (j in 1:length(t_sim)){
    foreach(k=1:n_sim) %dopar% {
      out <- simulate_deaths2(agg_df,v_e[[i]],t_sim[j],r,deaths_ratio[1],d[3],w[3],Q_a)
      # Save output
      saveRDS(out,file = paste0("../Data/",dir1,"pred_outcomes_",k,"_",v_e_nms[i],"_",t_sim[j],"days.RDS"))
    }
  }
}

tend <- Sys.time()
print(tend-tstart)

## Aggregate output
tstart <- Sys.time()

for (b in 1:length(by)){
  for (i in 1:length(v_e)){
    for (j in 1:length(t_sim)){
      for (l in 1:length(n_v)){
        tstart1 = Sys.time()
        lres <- foreach(k=1:n_sim) %dopar% {
          out <- readRDS(paste0("../Data/",dir1,"pred_outcomes_",k,"_",v_e_nms[i],"_",t_sim[j],"days.RDS"))
          
          # Process simulated output
          agg_out <- process_predictions_deaths2(out,grp,agg_df,v_e[[i]],v_ei[[i]],v_ec[[i]],r,r_i,t_sim,IFR_ratio,frlty_idx,d,w,u)
          
          # Run vaccination prioritization strategies
          run_prioritisation_strategies2(agg_out,n_v[l],agg_df,strategies,by[b],comorbs,dir1,t_sim[j],lbls,"YlGnBu",spec_pops)
        }
        
        # Row bind output lists
        res <- do.call(rbind,lapply(lres,"[[",1))
        lx <- vector("list",length(lres[[1]])-1)
        for (m in 1:(length(lres[[1]])-1)){
          lx[[m]] <- do.call(rbind,lapply(lres,"[[",m+1))
        }
        
        # Save output
        saveRDS(res,paste0("../Data/",dir1,"pred_vacc_impact_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by[b],".RDS"))
        saveRDS(lx,paste0("../Data/",dir1,"pred_vacc_impact_by_covar_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by[b],".RDS"))
        
        tend1 = Sys.time()
        print(tend1-tstart1)
      }
    }
  }
}

tend <- Sys.time()
print(tend-tstart)

## Calculate summary statistics and plot results
for (b in 1:length(by)){
  # Reorder population according to a priori optimal allocation order based on risk
  x <- copy(agg_df)
  HCW_LTCF_idx <- (x$special.population %in% c(1,3,8))
  under20_idx <- (x$age_cat %in% c("<10","10-19"))
  x_HCW_LTCF <- x[HCW_LTCF_idx & !under20_idx,]
  x_other <- x[!HCW_LTCF_idx & !under20_idx,]
  x_under20 <- x[under20_idx,]
  
  ord <- order_by_risk(x_other,by[b])
  x_other <- x_other[match(ord,x_other$grp),]
  x <- rbind(x_HCW_LTCF,x_other,x_under20)
  
  x[,comorb:=as.integer(rowSums2(as.matrix(x[,comorbs,with=F]))!=0)]
  x[is.na(population),population:=0]
  x[,cum_pop:=cumsum(population)]
  
  setDF(x)
  
  
  ## Plot outcomes for different prioritization strategies
  # Strategies and risk factors to plot performance for
  strategies1 <- strategies[1:(length(strategies)-3)]
  grp1 <- c(demogrphcs,"special.population","comorb")
  lbls1 <- c("County","Age","Sex","Race/ethnicity","Special population","Comorbidity")
  vrble <- c("case","death",by[b])
  
  for (i in 1:length(v_e)){
    for (j in 1:length(t_sim)){
      for (l in 1:length(n_v)){
        res <- readRDS(paste0("../Data/",dir1,"pred_vacc_impact_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by[b],".RDS"))
        lx <- readRDS(paste0("../Data/",dir1,"pred_vacc_impact_by_covar_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by[b],".RDS"))

        # Calculate summary statistics across simulations
        setDT(res)
        # names(res)[names(res)=="Strategy"] <- "strategy"
        cols1 <- names(res)[names(res)!="strategy"]
        res_ss <- calc_summary_stats(res,cols1,"strategy")
        write.csv(res_ss,paste0("../Data/",dir1,"pred_cases_deaths_",by[b],"_all_strategies_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by[b],".csv"),row.names = F)
        tbl1 <- res_ss[strategy!=0,.(`Strategy`=strategy,
                                     `Cases averted, mean (95% UI)`=mean_and_CI(cases_averted,cases_averted_q_95_LB,cases_averted_q_95_UB,d=2,method="signif"),
                                     `Percentage of cases averted, mean (95% UI)`=mean_and_CI(prop_cases_averted,prop_cases_averted_q_95_LB,prop_cases_averted_q_95_UB,f=100,d=0),
                                     `Deaths averted, mean (95% UI)`=mean_and_CI(deaths_averted,deaths_averted_q_95_LB,deaths_averted_q_95_UB,d=2,method="signif"),
                                     `Percentage of deaths averted, mean (95% UI)`=mean_and_CI(prop_deaths_averted,prop_deaths_averted_q_95_LB,prop_deaths_averted_q_95_UB,f=100,d=0),
                                     `DALYs averted, mean (95% UI)`=mean_and_CI(DALYs_averted,DALYs_averted_q_95_LB,DALYs_averted_q_95_UB,d=2,method="signif"),
                                     `Percentage of DALYs averted, mean (95% UI)`=mean_and_CI(prop_DALYs_averted,prop_DALYs_averted_q_95_LB,prop_DALYs_averted_q_95_UB,f=100,d=0))]
        write.csv(tbl1,paste0("../Data/",dir1,"tbl1_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by[b],".csv"),row.names = F)

        lx_ss <- vector("list",length(lx))
        for (m in 1:length(lx)){
          cols2 <- names(lx[[m]])[!(names(lx[[m]]) %in% c("strategy",grp1[m]))]
          lx_ss[[m]] <- calc_summary_stats(lx[[m]],cols2,c("strategy",grp1[m]))
        }
        lx_ss[[5]][,special.population := spec_pops[special.population+1]]

        # Plot results
        fdir <- paste0("../Figures/",dir1,t_sim[j],"days/",v_e_nms[i],"/",n_v[l],"doses/")
        plot_predictions_all_strtgs2(res_ss[res_ss$strategy %in% c(0,strategies1),],lx_ss,paste0(fdir,"Strategies1to",length(lbls)-1,by[b],"/"),vrble,lbls,"YlGnBu",strategies1,grp1,lbls1,by[b])

        # Plot optimal allocation
        # plot_optimal_allocation(x,n_v[l],paste0(fdir,"OptimalAllocation",by[b],"/"),spec_pops)
        plot_optimal_allocation(x[!(x$special.population %in% c(1,3,8)),],n_v[l],paste0(fdir,"OptimalAllocation",by[b],"/"),spec_pops)
      }
    }
  }
  
  
  ## Plot results of vaccine sensitivity analysis
  lbls2 <- c("95%","60%","20%") #,"60% age-dep")
  for (j in 1:length(t_sim)){
    lp <- vector("list",length(n_v))
    for (l in 1:length(n_v)){
      # Read in results for different vaccine efficacies
      lresSA <- vector("list",length(v_e_nms[1:3]))
      for (i in 1:length(v_e_nms[1:3])){
        lresSA[[i]] <- read.csv(paste0("../Data/",dir1,"pred_cases_deaths_",by[b],"_all_strategies_",n_v[l],"doses_",v_e_nms[i],"_",t_sim[j],"days",by[b],".csv"),stringsAsFactors = F)
        lresSA[[i]]$v_e <- v_e_nms[i]
      }
      # Bind results for different efficacies together
      res <- do.call(rbind,lresSA)

      # Plot results
      # fdir <- paste0("../Figures/",dir1,t_sim[j],"days/vacc_eff_SA/",n_v[l],"doses/")
      lp[[l]] <- plot_predictions_vacc_eff_SA(res,strategies1,v_e_nms[1:3],vrble,lbls,lbls2)
    }
    p <- lapply(lp,"[[",3)
    for (i in 1:(length(p)-1)){
      p[[i]] <- p[[i]] + theme(legend.position="none")
      p[[i+1]] <- p[[i+1]] + theme(axis.title.y=element_blank())
    }
    fdir <- paste0("../Figures/",dir1,t_sim[j],"days/vacc_eff_SA/")
    dir.create(fdir,recursive = T)
    ggsave(paste0(fdir,"pred_",by[b],"_by_strategy_vacc_eff_SA.pdf"),plot_grid(plotlist=p,nrow=1,rel_widths=c(1.05,1,1.15),labels=c("A","B","C","")),width = 15, height = 5)
  }
  
  
  ## Calculate and output impact of age-based targeting by age
  for (j in 1:length(t_sim)){
    x1 <- rbind(x_other,x_under20)
    setDT(x1)
    x1[,`:=`(deaths_averted=calc_outcome_averted(p_past_inf,lambda_adj,r,t_sim[j],population,v_e[[3]]),
             infections_averted=calc_outcome_averted(p_past_inf,lambda_inf,r_i,t_sim[j],population,v_ei[[3]]),
             cases_averted=calc_outcome_averted(p_past_inf,lambda_case,r_i,t_sim[j],population,v_ei[[3]]))]
    age_cat1 <- x1$age_cat
    age_cat1[age_cat1=="<10"] <- "0-9"
    age_cat1[age_cat1=="80+"] <- "80-89"
    age_lbs1 <- as.numeric(gsub("-.*|\\+|>","",age_cat1))
    age_ubs1 <- as.numeric(gsub(".*-|<","",age_cat1)) + 1
    age_mid1 <- (age_lbs1 + age_ubs1)/2
    x1[,`:=`(DALYs_averted=calc_DALYs_from_deaths(deaths_averted,life_expectancy,d[3],w[3])+
               calc_DALYs_from_infections_and_cases(infections_averted,cases_averted,d,w),
             QALYs_averted=calc_QALYs_from_deaths(deaths_averted,life_expectancy,age_mid1,Q_a)+
               calc_QALYs_from_infections_and_cases(infections_averted,cases_averted,u))]
    agg_x1 <- x1[,lapply(.SD,sum),.SDcols=c("population","deaths_averted","infections_averted","cases_averted","DALYs_averted","QALYs_averted"),by=.(age_cat)]
    agg_x1_long <- melt(agg_x1,id.vars = "age_cat")
    agg_x1_long[,prop_CA:=value/sum(value),by=.(variable)]
    agg_x1 <- dcast(agg_x1_long, age_cat ~ variable,value.var = c("value","prop_CA"))
    names(agg_x1) <- gsub("value_","",names(agg_x1))
    agg_x1 <- agg_x1[c((nrow(agg_x1)-1):1,nrow(agg_x1)),]
    write.csv(agg_x1,paste0("../Data/",dir1,"prop_CA_cases_deaths_DALYs_averted_by_age_",t_sim[j],"days",by[b],".csv"),row.names = F)
    tbl2 <- agg_x1[,.(`Age (years)`=age_cat,
                      `Population vaccinated, N (%)`=num_and_perc(population,prop_CA_population,d=0),
                      `Cases averted, N (%)`=num_and_perc(cases_averted,prop_CA_cases_averted,d=2,method="signif"),
                      `Deaths averted, N (%)`=num_and_perc(deaths_averted,prop_CA_deaths_averted,d=2,method="signif"),
                      `DALYs averted, N (%)`=num_and_perc(DALYs_averted,prop_CA_DALYs_averted,d=2,method="signif"))]
    write.csv(tbl2,paste0("../Data/",dir1,"tbl2_",t_sim[j],"days",by[b],".csv"),row.names = F)
  }
}

