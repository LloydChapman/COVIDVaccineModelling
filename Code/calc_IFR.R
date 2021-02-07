rm(list=ls())

source("seroprevalence_functions.R")

glm_fit <- readRDS("../Data/death_regression_output_2020-01-19_1.RDS")
x <- glm_fit$data

agg_x <- aggregate(cbind(population,cum_deaths) ~ age_cat,x,sum)
agg_x$age_cat[agg_x$age_cat=="<10"] <- "0-9"
agg_x$age_cat[agg_x$age_cat=="80+"] <- "80-100"
agg_x$age_low <- as.numeric(sub("-.*","",agg_x$age_cat))
agg_x$age_upp <- as.numeric(sub(".*-","",agg_x$age_cat))

# Read in age seroprevalence curve parameters
res <- readRDS("../Data/seroprev_pars.RDS")
pars <- res$par

agg_x$seroprev <- seroprevalence(agg_x,pars[1],pars[2])
agg_x$infections <- agg_x$seroprev * agg_x$population
print(sum(agg_x$infections))

agg_x$IFR <- agg_x$cum_deaths/agg_x$infections
saveRDS(agg_x,"../Data/CA_IFR_by_age.RDS")
