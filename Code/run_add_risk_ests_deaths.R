rm(list=ls())

library(reshape2)
library(data.table)

source("calc_hosp_ICU_death_risk.R")
source("processing_functions.R")
source("prediction_functions.R")

## Inputs
fnm <- "death_regression_output_2020-01-19_1.RDS"

# Relative risks for HCWs, prisoners, SNF residents, teachers and homeless individuals
RR <- c(HCW = 3,prisoner = 5.5,SNF = 9,teacher = 1.8,homeless = 1.65,essential = 1.7)

# Get hazard ratios for death by comorbidity
HR <- calc_hosp_ICU_death_risk()$HR

# Load processed synthetic population data frame
county_names <- read.csv("../Data/county_names.csv",stringsAsFactors = F)
df <- readRDS("../Data/CA_pop.RDS")

# # Subset to Alameda County for testing
# df1 <- df
# df <- df1[df1$county_res=="Alameda",]

# Convert to data table
setDT(df)

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
df <- add_risk_ests_deaths(paste0("../Data/",fnm),df,RR,HR,lt_long,seroprev,IFR_long,IFR_ratio)
saveRDS(df,"../Data/CA_pop_with_risk_ests_deaths.RDS")
