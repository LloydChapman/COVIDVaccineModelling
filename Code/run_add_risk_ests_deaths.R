rm(list=ls())

library(reshape2)
library(data.table)

source("calc_hosp_ICU_death_risk.R")
source("processing_functions.R")
source("prediction_functions.R")

## Inputs
# Read in death regression output
glm_fit <- readRDS("../Data/death_regression_output_2020-01-02_2.RDS")
alpha <- 0.05 # significance level for CI around death rate

# Load processed synthetic population data frame
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

# Get hazard ratios for death by comorbidity
HR <- calc_hosp_ICU_death_risk()$HR

# Read in life table data frame
lt_long <- as.data.table(read.csv("../Data/US_life_table.csv",stringsAsFactors = F))

# Read in seroprevalence data
seroprev_long <- read.csv("../Data/CA_seroprevalence.csv",stringsAsFactors = F)
seroprev <- aggregate(seroprev ~ group,seroprev_long,mean)
seroprev <- seroprev[!(seroprev$group %in% c("Cumulative","Female","Male")),]
names(seroprev)[names(seroprev)=="group"] <- "age_cat_sero"
setDT(seroprev)

# Read in IFR data
IFR_long <- readRDS("../Data/IFR_by_age_ODriscoll.RDS")
IFR_long$sex <- NULL
names(IFR_long)[names(IFR_long)=="grp"] <- "sex"
setDT(IFR_long)
IFR_ratio <- readRDS("../Data/IFR_ratio.RDS") # estimated ratio of CA IFR to O'Driscoll IFR 

## Add risk estimates from regression of CDPH data and hospitalisation and death probabilities to synthetic population data frame
df <- add_risk_ests_deaths(glm_fit,alpha,df,RR,RR_LB,RR_UB,HR,lt_long,seroprev,IFR_long,IFR_ratio)
saveRDS(df,"../Data/CA_pop_with_risk_ests_deaths1.RDS")
