rm(list=ls())

library(reshape2)
library(data.table)

source("calc_hosp_ICU_death_risk.R")
source("prediction_functions.R")

## Inputs
# Read in death regression output
glm_fit <- readRDS("../Data/death_regression_output_2020-02-05_2020-09-30_3.RDS")
alpha <- 0.05 # significance level for CI around death rate

# Load processed simulated population data frame
county_names <- read.csv("../Data/county_names.csv",stringsAsFactors = F)
df <- readRDS("../Data/CA_pop1.RDS")

# # Subset to Alameda County for testing
# df1 <- df
# df <- df1[df1$county_res=="Alameda",]

# Relative risks for HCWs, prisoners, SNF residents, teachers and homeless individuals
RR_essential <- 1.7 # relative risk for essential workers from Allen Nature Human Behaviour 2020
RR_essential_LB <- 1.1
RR_essential_UB <- 2.5
RR_f_vs_nf <- 1.7
N_f <- sum(df[,special.population]==6) # number of frontline essential workers
N_nf <- sum(df[,special.population]==7) # number of non-frontline essential workers
RR_nf <- RR_essential*(N_f+N_nf)/(RR_f_vs_nf*N_f+N_nf) # relative risk for non-frontline essential workers
RR_nf_LB <- RR_essential_LB*(N_f+N_nf)/(RR_f_vs_nf*N_f+N_nf)
RR_nf_UB <- RR_essential_UB*(N_f+N_nf)/(RR_f_vs_nf*N_f+N_nf)
RR <- c(HCW = 3.4,prisoner = 5.5,SNF = 6.4,educator = 1.8,homeless = 1.65,essential_f = RR_f_vs_nf*RR_nf,essential_nf = RR_nf,alf = 1.5)
RR_LB <- c(HCW = 3,prisoner = 5,SNF = 5,educator = 1.2,homeless = 1.5,essential_f = RR_f_vs_nf*RR_nf_LB,essential_nf = RR_nf_LB,alf = 1.4)
RR_UB <- c(HCW = 4,prisoner = 6,SNF = 8,educator = 2.8,homeless = 1.8,essential_f = RR_f_vs_nf*RR_nf_UB,essential_nf = RR_nf_UB,alf = 1.6)

# Relative risks of death for frontline and non-frontline essential workers from Chen medRxiv 2021
RR_death_essential <- data.table(special.population=c(6,7),RR=c(1.2,1.1),RR_LB=c(1,1),RR_UB=c(1.4,1.2))

# Get hazard ratios for death by comorbidity
HR <- calc_hosp_ICU_death_risk()$HR

# Frailty index for LTCF residents
frlty_idx <- 3

# Read in life table data frame
lt_long <- as.data.table(read.csv("../Data/US_life_table.csv",stringsAsFactors = F))

# Read in seroprevalence data
seroprev_long <- read.csv("../Data/CA_seroprevalence1.csv",stringsAsFactors = F)
seroprev <- aggregate(seroprev ~ group,seroprev_long,mean)
seroprev <- seroprev[!(seroprev$group %in% c("Cumulative","Female","Male")),]
names(seroprev)[names(seroprev)=="group"] <- "age_cat_sero"
setDT(seroprev)

# Read in IFR data
IFR_long <- readRDS("../Data/IFR_by_age_ODriscoll.RDS")
IFR_long$sex <- NULL
names(IFR_long)[names(IFR_long)=="grp"] <- "sex"
setDT(IFR_long)
IFR_ratio <- readRDS("../Data/IFR_ratio2.RDS") # estimated ratio of CA IFR to O'Driscoll IFR 

# Read in CDPH death data
deaths_age <- read.csv("../Data/CDPHDeathAgeData/CDPHDeathAgeData.csv",stringsAsFactors = F)

# Extract lower and upper bounds of age groups
age_cat_lbs <- seq(0,80,by=10)
lbls <- paste0(age_cat_lbs,"-",age_cat_lbs+9)
lbls[lbls=="0-9"] <- "<10"
lbls[lbls=="80-89"] <- "80+"
deaths_age$YEARS[deaths_age$YEARS==">100"] <- "100-105"
deaths_age$age_low <- as.numeric(sub("-.*","",deaths_age$YEARS))
deaths_age$age_upp <- as.numeric(sub(".*-","",deaths_age$YEARS)) - 1
# Correct age group labels
deaths_age$age_cat <- cut(deaths_age$age_low,c(age_cat_lbs,max(deaths_age$age_upp)+1),labels=lbls,right=F)

# Aggregate to age groups in CDPH line list data
agg_deaths_age <- aggregate(cbind(SNF.LTC.residents,Not.SNF.LTC.residents) ~ age_cat,deaths_age,sum)
agg_deaths_age <- rbind(list("<10",0,1,0),agg_deaths_age)
setDT(agg_deaths_age)

# Calculate LTCF and non-LTCF populations by age
df[,LTCF:=(special.population %in% c(3,8))]
agg_df_age <- df[,.(population=.N),by=.(age_cat,LTCF)]
agg_df_age_wide <- dcast(agg_df_age,age_cat ~ LTCF,value.var = "population")
setnames(agg_df_age_wide,c("FALSE","TRUE"),c("non_LTCF","LTCF"))
agg_df_age_wide[is.na(LTCF),LTCF:=0]

# Merge with population data
agg_deaths_age <- merge(agg_deaths_age,agg_df_age_wide,by="age_cat")

# Calculate RR of death for LTCF residents
agg_deaths_age[,RR_LTCF:=(SNF.LTC.residents/LTCF)/(Not.SNF.LTC.residents/non_LTCF)]
agg_deaths_age[is.na(RR_LTCF),RR_LTCF:=0]

# Calculate 95% CIs for RR of death for LTCF residents
agg_deaths_age[,`:=`(RR_LTCF_LB=calc_RR_CI(RR_LTCF,SNF.LTC.residents,LTCF,Not.SNF.LTC.residents,non_LTCF)$CI_low,RR_LTCF_UB=calc_RR_CI(RR_LTCF,SNF.LTC.residents,LTCF,Not.SNF.LTC.residents,non_LTCF)$CI_upp)]
agg_deaths_age[is.na(RR_LTCF_LB),RR_LTCF_LB:=0]
agg_deaths_age[is.na(RR_LTCF_UB),RR_LTCF_UB:=0]

# # Merge with simulated population data frame
# df <- merge(df,agg_deaths_age[,.(age_cat,RR_LTCF)],by="age_cat")

# Remove LTCF indicator from simulated population data frame
df[,LTCF:=NULL]

## Add risk estimates from regression of CDPH data and hospitalisation and death probabilities to synthetic population data frame
df <- add_risk_ests_deaths(glm_fit,alpha,df,RR,RR_LB,RR_UB,HR,lt_long,seroprev,IFR_long,IFR_ratio,agg_deaths_age,frlty_idx,RR_death_essential)#
saveRDS(df,"../Data/CA_pop_with_risk_ests_deaths_vldtn2.RDS")
