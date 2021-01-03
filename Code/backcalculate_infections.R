rm(list=ls())

source("Code/calc_hosp_ICU_death_risk.R")
source("Code/analysis_functions.R")

# library(EpiNow2)
library(ggplot2)
library(actuar)
library(surveillance)
library(lubridate)
library(abind)
library(reshape2)
library(broom)

# Read in aggregated death data
agg_deaths <- read.csv("~/Downloads/processed_death_data_2020-01-19_1.csv",stringsAsFactors = F)
agg_deaths$date_of_death <- as.Date(agg_deaths$date_of_death)
# Assign group IDs
unq_agg_deaths <- agg_deaths[!duplicated(agg_deaths[,c("county_res","age_cat","sex","race_ethnicity")]),c("county_res","age_cat","sex","race_ethnicity")]
unq_agg_deaths$grp <- 1:nrow(unq_agg_deaths)
agg_deaths <- merge(agg_deaths,unq_agg_deaths,all.x=T)
# Plot
ggplot(agg_deaths[agg_deaths$grp %in% 1:100,],aes(x=date_of_death,y=n_deaths,group=as.factor(grp),color=as.factor(grp))) + geom_line() + theme(legend.position = "none")

# Aggregate and plot by age
agg_deaths$age_cat[agg_deaths$age_cat %in% c("<10","10-19","20-29","30-39")] <- "<40"
agg_deaths_age <- aggregate(n_deaths ~ date_of_death + age_cat,agg_deaths,sum)
agg_deaths_age$grp <- match(agg_deaths_age$age_cat,unique(agg_deaths_age$age_cat))
ggplot(agg_deaths_age,aes(x=date_of_death,y=n_deaths,group=as.factor(grp),color=as.factor(grp))) + geom_line()

####
# Read in and process data from O'Driscoll Nature 2020
IFR <- read.csv("Data/IFR_by_age_ODriscoll.csv",stringsAsFactors = F)
IFR$Age_group[IFR$Age_group=="80+"] <- "80-100"
IFR$age_low <- as.numeric(sub("-.*","",IFR$Age_group))
IFR$age_upp <- as.numeric(sub(".*-","",IFR$Age_group))
age_cat_lbs <- c(0,40,seq(50,80,by=10))
# lbls <- paste0(age_cat_lbs,"-",age_cat_lbs+9)
# lbls[lbls=="0-9"] <- "<10"
# lbls[lbls=="80-89"] <- "80+"
lbls <- c("<40","40-49","50-59","60-69","70-79","80+")
IFR$age_cat <- cut(IFR$age_low,c(age_cat_lbs,max(IFR$age_upp)+1),labels=lbls,right=F)

# Calculate means across the two 5-year age groups in each 10-year age group
agg_IFR <- aggregate(.~age_cat,IFR[,c(2:(ncol(IFR)-3),ncol(IFR))],mean)

# Reshape to long format for merging
grp <- c("male","female","mean")
IFR_long <-  reshape(agg_IFR,varying = list(paste0("Median_perc","_",grp),paste0("CI_95_LB","_",grp),paste0("CI_95_UB","_",grp)),v.names = c("median_perc","CI_95_LB","CI_95_UB"),idvar = "age_cat",drop = c(paste0("Sero_study_range_LB_",grp),paste0("Sero_study_range_UB_",grp)),direction = "long")
IFR_long$grp <- grp[IFR_long$time]
IFR_long$sex <- NA
IFR_long$sex[IFR_long$grp=="female"] <- 0
IFR_long$sex[IFR_long$grp=="male"] <- 1

# Read in California IFR data from CDC seroprevalence estimates
CA_IFR <- readRDS("Data/CA_IFR_by_age.RDS")

# Sum deaths in each risk group over time
cum_deaths <- backcalculate_infections(agg_deaths,IFR_long)
# # Multiply estimated number of infections from O'Driscoll IFR by ratio of total estimated from CDC CA seroprevalence estimates to that from O'Driscoll, so that total numbers of infections agree
# cum_deaths$n <- cum_deaths$n * sum(CA_IFR$infections)/sum(cum_deaths$n)
print(sum(cum_deaths$n)) # 6001684
# Plot distribution of numbers infected in each risk group 
hist(cum_deaths$n)

## Backcalculate infections using Becker back-projection algorithm from surveillance package
# Read in incubation period and onset-to-death delay data
incubation_period <- readRDS("Data/incubation_period.rds")
onset_to_death_delay <- readRDS("Data/onset_to_death_delay.rds")

# Convolve distributions to get infection-to-death delay distribution
incubation_period_pdf <- function(x) dlnorm(x,incubation_period$mean,incubation_period$sd)
onset_to_death_delay_pdf <- function(y) dlnorm(y,onset_to_death_delay$mean,onset_to_death_delay$sd)
infection_to_death_delay_pdf <- function(z) integrate(function(x,z) onset_to_death_delay_pdf(z-x)*incubation_period_pdf(x),-Inf,Inf,z)$value
infection_to_death_delay_pdf <- Vectorize(infection_to_death_delay_pdf)

# Discretised infection-to-death delay distribution
dmax <- incubation_period$max + onset_to_death_delay$max
infection_to_death_delay_pmf <- c(0,sapply(0:(dmax-1),function(x){integrate(function(z) infection_to_death_delay_pdf(z),x,x+1)$value}))/integrate(function(z) infection_to_death_delay_pdf(z),0,dmax)$value

# Plot to check distribution
# set.seed(1)
# X <- rlnorm(1000,incubation_period$mean,incubation_period$sd)
# Y <- rlnorm(1000,onset_to_death_delay$mean,onset_to_death_delay$sd)
# Z <- X + Y
# # compare the methods
# hist(Z,freq=F,breaks=50, xlim=c(0,30))
# z <- seq(0,50,0.01)
# lines(z,infection_to_death_delay_pdf(z),col="red")
# lines(0:dmax,infection_to_death_delay_pmf,col="green") 

# Cast aggregate death data to wide format (death counts by date by age group)
agg_deaths_wide <- dcast(agg_deaths_age,date_of_death ~ grp,value.var = "n_deaths")
dates <- agg_deaths_wide$date_of_death
row.names(agg_deaths_wide) <- dates
agg_deaths_wide$date_of_death <- NULL

# Remove all-zero time series for now
# agg_deaths_wide <- agg_deaths_wide[,!apply(agg_deaths_wide,2,function(x){all(x==0)})]
# Convert to sts object
deaths_sts <- sts(agg_deaths_wide,start=c(year(min(dates)),yday(min(dates))),frequency = 365)
# Plot death time series
plot(deaths_sts,type=observed~time)

# Deconvolve death curve to infection curve (for those that died)
mean_incubation_period <- exp(incubation_period$mean+incubation_period$sd^2/2)
mean_onset_to_death_delay <- exp(onset_to_death_delay$mean+onset_to_death_delay$sd^2/2)
mean_infection_to_death_delay <- round(mean_incubation_period + mean_onset_to_death_delay)
control <- list(Tmark = nrow(deaths_sts) - mean_infection_to_death_delay,eq3a.method = "C")
deaths_stsBP <- backprojNP(deaths_sts,infection_to_death_delay_pmf,control = control)
plot(deaths_stsBP,xaxis.labelFormat=NULL,legend=NULL,lwd=c(1,1,2),lty=c(1,1,1),ylim=c(0,100),main="")

# Scale by IFR to get infection curve (for all those infected)
infections_dead <- deaths_stsBP@upperbound
infections_dead <- infections_dead[1:(nrow(infections_dead)-mean_infection_to_death_delay),]
infections <- t(t(infections_dead)/(agg_IFR$Median_perc_mean/100))
infections_long <- melt(infections)
# Plot
ggplot(infections_long,aes(x=Var1,y=value,group=as.factor(Var2),color=as.factor(Var2))) + geom_line()

# Convolve to get unobserved case curve (i.e. the case curve assuming everyone had been tested)
# Placeholder values for reporting delay - [ ] UPDATE WITH ESTIMATES FROM CDPH DATA
rdmax <- 15
mean_log_reporting_delay <- log(3)
sd_log_reporting_delay <- 1
reporting_delay_pmf <- c(0,(plnorm(1:rdmax,mean_log_reporting_delay,sd_log_reporting_delay)-plnorm(0:(rdmax-1),mean_log_reporting_delay,sd_log_reporting_delay))/plnorm(rdmax,mean_log_reporting_delay,sd_log_reporting_delay))

cases <- matrix(nrow = nrow(infections),ncol = ncol(infections))
for (j in 1:ncol(infections)){
  for (i in 1:nrow(infections)){
    cases[i,j] <- sum(sapply(1:ifelse(i<=rdmax,i,rdmax+1),function(k){infections[i-k+1,j] * reporting_delay_pmf[k]}))
  }
  plot(infections[,j])
  points(cases[,j],col="red")
}

# Divide observed case curve by unobserved case curve
agg_cases_age <- read.csv("~/Downloads/agg_cases_age.csv",stringsAsFactors = F)
agg_cases_age$first_report_date <- as.Date(agg_cases_age$first_report_date)
x <- expand.grid(first_report_date=seq.Date(min(agg_cases_age$first_report_date),max(agg_cases_age$first_report_date),by=1),age_cat=unique(agg_cases_age$age_cat),stringsAsFactors = F)
agg_cases_age <- merge(agg_cases_age,x,all=T)
agg_cases_age$n[is.na(agg_cases_age$n)] <- 0
agg_cases_age$n_deaths[is.na(agg_cases_age$n_deaths)] <- 0
# Aggregate cases for age groups <40
agg_cases_age$age_cat[agg_cases_age$age_cat %in% c("<10","10-19","20-29","30-39")] <- "<40"
agg_cases_age <- aggregate(cbind(n,n_deaths) ~ first_report_date + age_cat,agg_cases_age,sum)
names(agg_cases_age)[names(agg_cases_age)=="first_report_date"] <- "date"
# Plot
pdf("Figures/Backcalculation/case_curves_by_age.pdf",width = 6,height=4)
print(ggplot(agg_cases_age,aes(x=date,y=n,group=age_cat,color=age_cat)) + geom_line() + xlab("Date") + ylab("Cases") + labs(color = "Age"))
dev.off()

# Melt data frame with estimated case counts to long format
cases_long <- melt(cases,varnames = c("date","age_cat"),value.name = "n_est")
cases_long$date <- rownames(agg_deaths_wide)[cases_long$date]
cases_long$age_cat <- lbls[cases_long$age_cat]

# Calculate the proportion of cases reported
agg_cases_age <- merge(agg_cases_age,cases_long,by=c("date","age_cat"),all.x=T)
agg_cases_age$prop_reported <- agg_cases_age$n/agg_cases_age$n_est
agg_cases_age$prop_reported[is.infinite(agg_cases_age$prop_reported)|is.nan(agg_cases_age$prop_reported)|agg_cases_age$prop_reported>1] <- 0
agg_cases_age$prop_reported[agg_cases_age$prop_reported>1] <- 1
# Plot
pdf("Figures/Backcalculation/prop_reported_over_time_by_age.pdf",width = 6,height=4)
ggplot(agg_cases_age[agg_cases_age$date>=as.Date("2020-04-01"),],aes(x=date,y=prop_reported,group=age_cat,color=age_cat)) + geom_line() + xlab("Date") + ylab("Proportion of infections reported") + labs(color = "Age")
dev.off()

# Read in aggregated test data
agg_tests_age_result <- read.csv("~/Downloads/agg_tests_age.csv",stringsAsFactors = F)
agg_tests_age_result$result_date <- as.Date(agg_tests_age_result$result_date)
# Combine age groups <40
agg_tests_age_result$age_cat[agg_tests_age_result$age_cat %in% c("<10","10-19","20-29","30-39")] <- "<40"
# Sum over all test results
agg_tests_age <- aggregate(n_tests ~ result_date + age_cat,agg_tests_age_result,sum)
names(agg_tests_age)[names(agg_tests_age)=="result_date"] <- "date"
# Plot tests over time by age
pdf("Figures/Backcalculation/tests_over_time_by_age.pdf",width = 6,height=4)
ggplot(agg_tests_age[agg_tests_age$date>=as.Date("2020-04-01"),],aes(x=date,y=n_tests,group=age_cat,color=age_cat)) + geom_line() + xlab("Date") + ylab("Tests") + labs(color = "Age")
dev.off()

# Merge with case counts data frame
agg_cases_age <- merge(agg_cases_age,agg_tests_age,by=c("date","age_cat"),all.x=T)
agg_cases_age$n_tests[is.na(agg_cases_age$n_tests)] <- 0
# Plot
pdf("Figures/Backcalculation/prop_reported_vs_tests_by_age.pdf",width = 6,height=4)
ggplot(agg_cases_age[agg_cases_age$date>=as.Date("2020-04-01"),],aes(x=n_tests,y=prop_reported,group=age_cat,color=age_cat)) + geom_point() + xlab("Tests") + ylab("Proportion of infections reported") + labs(color = "Age")
dev.off()

# Load population data
agg_pop <- read.csv("Data/agg_pop1.csv",stringsAsFactors = F)
# Combine age groups <40
agg_pop$age_cat[agg_pop$age_cat %in% c("<10","10-19","20-29","30-39")] <- "<40"
agg_pop <- aggregate(population ~ age_cat,agg_pop,sum)

# Merge with case counts data frame
agg_cases_age <- merge(agg_cases_age,agg_pop,by="age_cat",all.x=T)
# Calculate number of tests per head of population in each age group
agg_cases_age$n_tests_pp <- agg_cases_age$n_tests/agg_cases_age$population
# Plot
pdf("Figures/Backcalculation/prop_reported_vs_tests_per_person_by_age.pdf",width = 6,height=4)
ggplot(agg_cases_age[agg_cases_age$date>=as.Date("2020-04-01"),],aes(x=n_tests_pp,y=prop_reported,group=age_cat,color=age_cat)) + geom_point() + xlab("Tests per person") + ylab("Proportion of infections reported") + labs(color = "Age")
dev.off()

# # Merge with infections and deaths data frame
# inc <- merge(cum_deaths,agg_pop,by=c("county_res","age_cat","sex","race_ethnicity"),all.x = T)
# inc$exp_time <- inc$population * as.numeric(as.Date("2020-10-19")-as.Date("2020-01-19"))/1e5
# 
# # Fit Poisson regression model to estimated infections
# glm_fit <- glm(round(n) ~ offset(log(exp_time)) + county_res + age_cat + sex + race_ethnicity,family = poisson(link = log),data = inc)
# saveRDS(glm_fit,"regression_output_death_IFR_model.RDS")
# coeffs <- exp(coef(glm_fit))
# print(coeffs)
# 
# process_regression_output(dir,"regression_output_death_IFR_model.RDS")
# 
# # Load coefficient estimates from 3-month case model
# glm_fit1 <- readRDS("regression_output_2020-07-19_1.RDS")
# coeffs1 <- exp(coef(glm_fit1))
# 
# # Fit Poisson regression model to cumulative cases
# inc1 <- glm_fit1$data
# inc1$exp_time <- inc1$population * inc1$time[1]/1e5
# glm_fit2 <- glm(cum_cases ~ offset(log(exp_time)) + county_res + age_cat + sex + race_ethnicity,poisson(link = log),data = inc1)
# coeffs2 <- exp(coef(glm_fit2))
# print(coeffs2)
# 
# # Combine coefficient estimates into single data frame
# coeffs_comp <- data.frame(names = names(coeffs),death_IFR_model=coeffs,three_month_PH_model=coeffs1,three_month_Poiss_model=coeffs2)
# coeffs_comp_long <- melt(coeffs_comp)
# 
# # Plot coefficients to compare them
# pdf("../Figures/model_coefficients_comparison.pdf",width = 9,height = 5)
# ggplot(coeffs_comp_long,aes(x=names,y=value,group=variable,color=as.factor(variable))) + geom_point() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + xlab("coefficient") + labs(color="Model")
# dev.off()
