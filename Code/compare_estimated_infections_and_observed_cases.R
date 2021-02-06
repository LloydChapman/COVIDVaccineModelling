rm(list=ls())

library(reshape2)
library(data.table)
library(ggplot2)
library(broom)

source("calc_hosp_ICU_death_risk.R")
source("analysis_functions.R")
source("prediction_functions.R")

# Read in death data
# glm_fit <- readRDS("../Data/death_regression_output_2020-02-05_3.RDS")
# cum_deaths <- glm_fit$data
# cum_deaths$n_deaths <- NULL
# names(cum_deaths)[names(cum_deaths)=="cum_deaths"] <- "n_deaths"
cum_deaths <- read.csv("../Data/cum_deaths.csv",stringsAsFactors = F)

# Read in and process data from O'Driscoll Nature 2020
IFR <- read.csv("../Data/IFR_by_age_ODriscoll.csv",stringsAsFactors = F)
IFR$Age_group[IFR$Age_group=="80+"] <- "80-100"
IFR$age_low <- as.numeric(sub("-.*","",IFR$Age_group))
IFR$age_upp <- as.numeric(sub(".*-","",IFR$Age_group))
age_cat_lbs <- seq(0,80,by=10)
lbls <- paste0(age_cat_lbs,"-",age_cat_lbs+9)
lbls[lbls=="0-9"] <- "<10"
lbls[lbls=="80-89"] <- "80+"
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

# Save
saveRDS(IFR_long,"../Data/IFR_by_age_ODriscoll.RDS")

# Read in CDPH death age data with LTCF/non-LTCF breakdown
deaths_age <- read.csv("../Data/CDPHDeathAgeData/CDPHDeathAgeData.csv",stringsAsFactors = F)

# Extract lower and upper bounds of age groups
deaths_age$YEARS[deaths_age$YEARS==">100"] <- "100-105"
deaths_age$age_low <- as.numeric(sub("-.*","",deaths_age$YEARS))
deaths_age$age_upp <- as.numeric(sub(".*-","",deaths_age$YEARS)) - 1
# Correct age group labels
deaths_age$age_cat <- cut(deaths_age$age_low,c(age_cat_lbs,max(deaths_age$age_upp)+1),labels=lbls,right=F)

agg_deaths_age <- aggregate(cbind(SNF.LTC.residents,Not.SNF.LTC.residents) ~ age_cat,deaths_age,sum)
agg_deaths_age$LTCFprop <- agg_deaths_age$SNF.LTC.residents/(agg_deaths_age$SNF.LTC.residents + agg_deaths_age$Not.SNF.LTC.residents)
agg_deaths_age <- rbind(list("<10",NA,NA,0),agg_deaths_age)

# Calculate proportion of overall deaths that were LTCF residents
LTCFprop <- sum(deaths_age$SNF.LTC.residents)/sum(deaths_age$SNF.LTC.residents + deaths_age$Not.SNF.LTC.residents)
print(LTCFprop) #0.3240

# Merge with CDPH death data
cum_deaths <- merge(cum_deaths,agg_deaths_age[,c("age_cat","LTCFprop")],by="age_cat",all.x=T)
cum_deaths$n_LTCF_deaths <- cum_deaths$LTCFprop * cum_deaths$n_deaths
cum_deaths$total_deaths <- cum_deaths$n_deaths
cum_deaths$n_deaths <- cum_deaths$n_deaths - cum_deaths$n_LTCF_deaths

# Backcalculate infections in each risk group using O'Driscoll IFR
cum_deaths <- backcalculate_infections(cum_deaths,IFR_long)
# # Multiply estimated number of infections from O'Driscoll IFR by ratio of total estimated from CDC CA seroprevalence estimates to that from O'Driscoll, so that total numbers of infections agree
# cum_deaths$n <- cum_deaths$n * sum(CA_IFR$infections)/sum(cum_deaths$n)
print(sum(cum_deaths$n)) #6425783
# Plot distribution of numbers infected in each risk group 
hist(cum_deaths$n)

# Calculate cases in each risk group using Davies age-dependent clinical fraction
p_clin <- calc_hosp_ICU_death_risk()$p_clin_given_infctn
p_clin <- data.frame(age_cat=names(p_clin),p_clin=p_clin)
cum_deaths <- calc_cases_from_infections(cum_deaths,p_clin)

# Aggregate estimated infections by age
cum_inf_age <- aggregate(cbind(n,n_cases) ~ age_cat,cum_deaths,sum)

# Read in aggregated counts of confirmed cases up to Oct 22, 2020, in CDPH data
cum_cases <- read.csv("../Data/cum_cases.csv",stringsAsFactors = F)

# Aggregate over age
cum_cases_age <- aggregate(n ~ age_cat,cum_cases,sum)

# Merge estimated infections and observed cases
cum_inf_cases_age <- merge(cum_inf_age,cum_cases_age,by="age_cat",suffixes = c("","_cases_obs"))
cum_inf_cases_age$prop_inf <- cum_inf_cases_age$n/sum(cum_inf_cases_age$n)
cum_inf_cases_age$prop_cases <- cum_inf_cases_age$n_cases/sum(cum_inf_cases_age$n_cases)
cum_inf_cases_age$prop_cases_obs <- cum_inf_cases_age$n_cases_obs/sum(cum_inf_cases_age$n_cases_obs)

# Calculate multipler for IFR to match estimated total number of clinical cases with total number of confirmed cases
total_cases_obs <- sum(cum_inf_cases_age$n_cases_obs) #865563
total_cases_SNF <- 27019 # from LA Times rip of CDPH data
total_cases_ALF <- 10677 # from LA Times rip of CDSS data
total_cases_non_LTCF <- total_cases_obs - total_cases_SNF - total_cases_ALF
total_cases_est <- sum(cum_inf_cases_age$n_cases) #2186465
IFR_ratio <- total_cases_est/total_cases_non_LTCF #2.641083
saveRDS(IFR_ratio,"../Data/IFR_ratio2.RDS")

# Plot observed case counts vs estimated case counts by age
cum_inf_cases_age_long <- reshape2::melt(cum_inf_cases_age[,c("age_cat","n_cases_obs","n_cases")])
pdf("../Figures/Backcalculation/obsvd_vs_estd_case_age_counts2.pdf",width = 6,height = 4)
ggplot(cum_inf_cases_age_long,aes(x=age_cat,y=value,group=variable,fill=variable)) + geom_bar(stat="identity",position="dodge") + xlab("Age") + ylab("Cases") + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_fill_discrete(name="",labels=c("obsvd","estd"))
dev.off()

# Plot observed case counts vs corrected estimated case counts by age
cum_inf_cases_age_long$value[cum_inf_cases_age_long$variable=="n_cases"] <- cum_inf_cases_age_long$value[cum_inf_cases_age_long$variable=="n_cases"]/IFR_ratio
pdf("../Figures/Backcalculation/obsvd_vs_estd_case_age_counts_crrctd2.pdf",width = 6,height = 4)
ggplot(cum_inf_cases_age_long,aes(x=age_cat,y=value,group=variable,fill=variable)) + geom_bar(stat="identity",position="dodge") + xlab("Age") + ylab("Cases") + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_fill_discrete(name="",labels=c("obsvd","estd"))
dev.off()

# Plot observed and estimated age distributions of cases and infections
cum_inf_cases_age_long1 <- reshape2::melt(cum_inf_cases_age[,c("age_cat","prop_inf","prop_cases_obs","prop_cases")])
pdf("../Figures/Backcalculation/obsvd_vs_estd_case_age_distns2.pdf",width = 6,height = 4)
ggplot(cum_inf_cases_age_long1,aes(x=age_cat,y=value,group=variable,fill=variable)) + geom_bar(stat="identity",position="dodge") + xlab("Age") + ylab("Proportion") + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_fill_discrete(name="",labels=c("estd infctns","obsvd cases","estd cases"))
dev.off()

# Read in CDC seroprevalence data
seroprev_long <- read.csv("../Data/CA_seroprevalence1.csv",stringsAsFactors = F)
seroprev <- aggregate(seroprev ~ group,seroprev_long,mean)
seroprev <- seroprev[!(seroprev$group %in% c("Cumulative","Female","Male")),]
names(seroprev)[names(seroprev)=="group"] <- "age_cat_sero"
setDT(seroprev)

# Read in synthetic population
df <- readRDS("../Data/CA_pop_with_risk_ests_deaths3.RDS")
df$population <- 1
agg_df_age <- aggregate(population ~ age_cat_sero,df,sum)
rm(df)
gc()
agg_df_age <- merge(agg_df_age,seroprev,by="age_cat_sero")
# agg_df_age <- aggregate(seroprev ~ age_cat_sero,df,sum)

# Calculate expected number seropositive in each serological age category
agg_df_age$n_sero <- agg_df_age$population * agg_df_age$seroprev

# Calculate estimated cumulative number of infections from seroprevalence data and backcalculated from deaths using IFR
total_inf_sero <- sum(agg_df_age$n_sero) #1927303
total_inf_IFR <- sum(cum_inf_cases_age$n) #6425783
inf_ratio <- total_inf_IFR/total_inf_sero #3.33408
print(inf_ratio) 

# Estimate number of infections in each serological age category for backcalculated infection numbers 
agg_df_age$n <- NA
agg_df_age$n[agg_df_age$age_cat_sero=="0-17"] <- cum_inf_cases_age$n[cum_inf_cases_age$age_cat=="<10"] + 0.8*cum_inf_cases_age$n[cum_inf_cases_age$age_cat=="10-19"]
agg_df_age$n[agg_df_age$age_cat_sero=="18-49"] <- 0.2*cum_inf_cases_age$n[cum_inf_cases_age$age_cat=="10-19"] + cum_inf_cases_age$n[cum_inf_cases_age$age_cat=="20-29"] + cum_inf_cases_age$n[cum_inf_cases_age$age_cat=="30-39"] + cum_inf_cases_age$n[cum_inf_cases_age$age_cat=="40-49"]
agg_df_age$n[agg_df_age$age_cat_sero=="50-64"] <- cum_inf_cases_age$n[cum_inf_cases_age$age_cat=="50-59"] + 0.5*cum_inf_cases_age$n[cum_inf_cases_age$age_cat=="60-69"]
agg_df_age$n[agg_df_age$age_cat_sero=="65+"] <- 0.5*cum_inf_cases_age$n[cum_inf_cases_age$age_cat=="60-69"] + cum_inf_cases_age$n[cum_inf_cases_age$age_cat=="70-79"] + cum_inf_cases_age$n[cum_inf_cases_age$age_cat=="80+"]

# Plot observed vs estimated numbers of infections by age
agg_df_age_long <- reshape2::melt(agg_df_age,id.vars="age_cat_sero",measure.vars=c("n_sero","n"))
pdf("../Figures/Backcalculation/infections_seroprev_vs_IFR3.pdf",width = 6,height = 4)
ggplot(agg_df_age_long,aes(x=age_cat_sero,y=value,fill=variable)) + geom_bar(stat="identity",position="dodge") + xlab("Age") + ylab("Infections") + scale_fill_discrete(name="",labels=c("seroprev","IFR"))
dev.off()

# Plot observed vs corrected estimated numbers of infections by age
agg_df_age_long$value[agg_df_age_long$variable=="n"] <- agg_df_age_long$value[agg_df_age_long$variable=="n"]/IFR_ratio
pdf("../Figures/Backcalculation/infections_seroprev_vs_IFR_crrctd3.pdf",width = 6,height = 4)
ggplot(agg_df_age_long,aes(x=age_cat_sero,y=value,fill=variable)) + geom_bar(stat="identity",position="dodge") + xlab("Age") + ylab("Infections") + scale_fill_discrete(name="",labels=c("seroprev","IFR"))
dev.off()

# Apply IFR correction to estimated infections and cases
cum_deaths$n <- cum_deaths$n/IFR_ratio
cum_deaths$n_cases <- cum_deaths$n_cases/IFR_ratio

# Load population data
agg_pop <- read.csv("../Data/agg_pop1.csv",stringsAsFactors = F)
# Merge with infections and deaths data frame
x <- expand.grid(county_res=unique(cum_deaths$county_res),age_cat=unique(cum_deaths$age_cat),sex=c(0,1),race_ethnicity=unique(cum_deaths$race_ethnicity),stringsAsFactors = F)
x <- x[x$race_ethnicity!="Unknown",]
demogrphcs <- c("county_res","age_cat","sex","race_ethnicity")
cum_deaths <- merge(x,cum_deaths,by=demogrphcs,all.x = T)
cum_deaths$n_cases[is.na(cum_deaths$n_cases)] <- 0
inc <- merge(cum_deaths,agg_pop,by=demogrphcs,all.x = T)
inc$exp_time <- inc$population * as.numeric(as.Date("2020-10-22")-as.Date("2020-01-27"))/1e5 # assumes 9-day reporting to death lag with first death on 2020-02-05

# # Fit Poisson regression model to estimated infections
# glm_fit_inf <- glm(round(n) ~ offset(log(exp_time)) + county_res + age_cat + sex + race_ethnicity,family = poisson(link = log),data = inc)
# saveRDS(glm_fit_inf,"../Data/regression_output_death_IFR_model_infections.RDS")
# HR_inf <- data.frame(Estimate=exp(coef(glm_fit_inf)),exp(confint.default(glm_fit_inf)))
# # process_regression_output("../Data/","regression_output_death_IFR_model_infections.RDS")

# Fit Poisson regression model to estimated cumulative cases
inc$age_cat <- factor(inc$age_cat,levels = c("50-59",lbls[lbls!="50-59"]))
glm_fit1 <- glm(round(n_cases) ~ offset(log(exp_time)) + county_res + age_cat + sex + race_ethnicity,family = poisson(link = log),data = inc)
saveRDS(glm_fit1,"../Data/regression_output_death_IFR_model1.RDS")
HR <- data.frame(Estimate=exp(coef(glm_fit1)),exp(confint.default(glm_fit1)))
HR$mod <- "Estd"
HR$names <- row.names(HR)
# process_regression_output("../Data/","regression_output_death_IFR_model.RDS")

# Load data
# glm_fit2 <- readRDS("../Data/regression_output_2020-07-19_1.RDS")
# glm_fit2 <- readRDS("../Data/regression_output_2020-01-19_1.RDS")
cum_cases <- merge(x,cum_cases,by=demogrphcs,all.x = T)
cum_cases$n[is.na(cum_cases$n)] <- 0

# Fit Poisson regression model to cumulative cases
inc1 <- merge(cum_cases,agg_pop,by=demogrphcs,all.x = T)
inc1$age_cat <- factor(inc1$age_cat,levels = c("50-59",lbls[lbls!="50-59"]))
inc1$exp_time <- inc1$population * as.numeric(as.Date("2020-10-22")-as.Date("2020-01-27"))/1e5 # assumes 9-day reporting to death lag with first death on 2020-02-05
glm_fit3 <- glm(n ~ offset(log(exp_time)) + county_res + age_cat + sex + race_ethnicity,poisson(link = log),data = inc1)
HR1 <- data.frame(Estimate=exp(coef(glm_fit3)),exp(confint.default(glm_fit3)))
HR1$mod <- "Obsvd"
HR1$names <- row.names(HR1)

# Compare coefficients
HR_comp <- rbind(HR1,HR)
HR_comp$mod <- factor(HR_comp$mod,levels = c("Obsvd","Estd"))

# Plot
pdf("../Figures/Backcalculation/param_ests_comparison2.pdf",width = 9,height = 6)
ggplot(HR_comp,aes(x=names,y=Estimate,group=mod,color=as.factor(mod))) + geom_point() + geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..)) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + xlab("Parameter") + ylab("Value") + labs(color="Cases")
dev.off()
