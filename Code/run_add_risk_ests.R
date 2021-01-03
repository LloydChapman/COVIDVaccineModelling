rm(list=ls())

library(reshape2)

source("~/Dropbox/COVIDVaccineModelling/Code/processing_functions.R")
source("~/Dropbox/COVIDVaccineModelling/Code/prediction_functions.R")

# Set working directory
dir <- "~/Dropbox/COVIDVaccineModelling/Data/"
setwd(dir)

## Inputs
fnm <- "regression_output_2020-07-19_1.RDS"

# Relative risks for HCWs, prisoners, SNF residents, teachers and homeless individuals
# RR <- c(HCW = 3,prisoner = 5.5,SNF = 9,teacher = 1.8,homeless = 1.65)
RR <- c(HCW = 3,prisoner = 5.5,SNF = 9,teacher = 1.8,homeless = 1.65,essential = 1.7)

# Load processed synthetic population data frame
county_names <- read.csv("county_names.csv",stringsAsFactors = F)
# df <- load_synth_pop(paste0(dir,"Simulation1126/"),county_names)
df <- readRDS("CA_pop.RDS")

# # Subset to Alameda County for testing
# df1 <- df
# df <- df1[df1$county_res=="Alameda",]

# Read in death risk estimates from regression of CDPH data
p_death <- readRDS("prob_death_given_clin_by_age_and_sex_2020-07-19_1.RDS")

# Read in life table data frame
lt_long <- read.csv("US_life_table.csv",stringsAsFactors = F)

# Read in seroprevalence data
seroprev_long <- read.csv("CA_seroprevalence.csv",stringsAsFactors = F)
seroprev <- aggregate(seroprev ~ group,seroprev_long,mean)
seroprev <- seroprev[!(seroprev$group %in% c("Cumulative","Female","Male")),]
names(seroprev)[names(seroprev)=="group"] <- "age_cat_sero"

## Add risk estimates from regression of CDPH data and hospitalisation and death probabilities to synthetic population data frame
df <- add_risk_ests(fnm,RR,df,p_death,lt_long,seroprev)
saveRDS(df,"CA_pop_with_risk_ests1.RDS")