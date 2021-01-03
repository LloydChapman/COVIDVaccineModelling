library(reshape2)
library(ggplot2)

setwd("~/Dropbox/COVIDVaccineModelling/Data/california-coronavirus-data/")

calc_new_cases <- function(x){c(0,diff(x))}

calc_exposure_time <- function(x){
  idx <- (x$date==max(x$date))
  sus_exp_time <- (x$population[idx]-x$confirmed_cases_total[idx])*x$time[idx] # exposure time of individuals still susceptible
  inf_exp_time <- sum(x$new_cases * x$time) # exposure time of individuals that have been infected
  exp_time <- (sus_exp_time + inf_exp_time)/1e5 # total exposure time in units of 100,000 person-days
}

# COVID-19 case data
# State total
cdph_state_total <- read.csv("cdph-state-totals.csv",stringsAsFactors = F)
cdph_state_total$date <- as.Date(cdph_state_total$date)
cdph_state_total <- cdph_state_total[order(cdph_state_total$date),]

# Correct error in 18-49 category
cdph_state_total$age_0_to_17[cdph_state_total$age_0_to_17==80869] <- 80069 # bit of a guess for what the correct number is!
cdph_state_total$age_18_to_49[cdph_state_total$age_18_to_49==203739] <- 303739
cdph_state_total$age_50_to_64[cdph_state_total$age_50_to_64==147760] <- 148760
cdph_state_total$age_65_and_up[cdph_state_total$age_65_and_up==52138] <- 53138

cdph_state_total$new_cases <- calc_new_cases(cdph_state_total$confirmed_cases)
cdph_state_total$age_0_to_17_new_cases <- calc_new_cases(cdph_state_total$age_0_to_17)
cdph_state_total$age_18_to_49_new_cases <- calc_new_cases(cdph_state_total$age_18_to_49)
cdph_state_total$age_50_to_64_new_cases <- calc_new_cases(cdph_state_total$age_50_to_64)
cdph_state_total$age_65_and_up_new_cases <- calc_new_cases(cdph_state_total$age_65_and_up)

# which(cdph_state_total$age_0_to_17_new_cases<0)
# which(cdph_state_total$age_18_to_49_new_cases<0)
# which(cdph_state_total$age_50_to_64_new_cases<0)
# which(cdph_state_total$age_65_and_up_new_cases<0)

# clrs <- c("All"="black","0-17"="red","18-49"="green","50-64"="blue","65+"="orange")
ggplot(cdph_state_total) + geom_line(aes(x=date,y=new_cases,col="black")) + geom_line(aes(x=date,y=age_0_to_17_new_cases,col="blue")) +
  geom_line(aes(x=date,y=age_18_to_49_new_cases,col="green")) + geom_line(aes(x=date,y=age_50_to_64_new_cases,col="orange")) +
  geom_line(aes(x=date,y=age_65_and_up_new_cases,col="red")) + scale_color_manual(name = "Age",values=c("black","blue","green","orange","red"),labels=c("All","0-17","18-49","50-64","65+"))

# Age groups - data only goes up to the end of August and seems to have quite a few errors, so not that useful
cdph_age <- read.csv("cdph-age.csv",stringsAsFactors = F)
cdph_age$date <- as.Date(cdph_age$date)
cdph_age$age_grp <- "age_0_to_17"
cdph_age$age_grp[cdph_age$age %in% c("18-34","35-49")] <- "age_18_to_49"
cdph_age$age_grp[cdph_age$age %in% c("50-59","60-64")] <- "age_50_to_64"
cdph_age$age_grp[cdph_age$age %in% c("65-69","70-74","75-79","80+")] <- "age_65_and_up"

# Aggregate age groups to cross-check with state dataset
agg_cdph_age <- aggregate(cbind(confirmed_cases_total,deaths_total) ~ date + age_grp,cdph_age,sum)
agg_cdph_age_wide <- dcast(agg_cdph_age,date ~ age_grp, value.var = "confirmed_cases_total")
agg_cdph_age_wide$age_0_to_17_new_cases <- calc_new_cases(agg_cdph_age_wide$age_0_to_17)
agg_cdph_age_wide$age_18_to_49_new_cases <- calc_new_cases(agg_cdph_age_wide$age_18_to_49)
agg_cdph_age_wide$age_50_to_64_new_cases <- calc_new_cases(agg_cdph_age_wide$age_50_to_64)
agg_cdph_age_wide$age_65_and_up_new_cases <- calc_new_cases(agg_cdph_age_wide$age_65_and_up)
agg_cdph_age_wide$new_cases <- apply(agg_cdph_age_wide[,6:ncol(agg_cdph_age_wide)],1,sum)
  
ggplot(agg_cdph_age_wide[apply(agg_cdph_age_wide[,6:ncol(agg_cdph_age_wide)],1,function(x) all(x>=0 & x<7.5e4)),]) + geom_line(aes(x=date,y=new_cases,col="black")) + geom_line(aes(x=date,y=age_0_to_17_new_cases,col="blue")) +
  geom_line(aes(x=date,y=age_18_to_49_new_cases,col="green")) + geom_line(aes(x=date,y=age_50_to_64_new_cases,col="orange")) +
  geom_line(aes(x=date,y=age_65_and_up_new_cases,col="red")) + scale_color_manual(name = "Age",values=c("black","blue","green","orange","red"),labels=c("All","0-17","18-49","50-64","65+"))

# Race/ethnicity
cdph_race_ethnicity <- read.csv("cdph-race-ethnicity.csv",stringsAsFactors = F)
# Remove extra aggregated age categories
cdph_race_ethnicity <- cdph_race_ethnicity[!(cdph_race_ethnicity$age %in% c("18+","all")),]
cdph_race_ethnicity$date <- as.Date(cdph_race_ethnicity$date)
cdph_race_ethnicity <- cdph_race_ethnicity[order(cdph_race_ethnicity$date),]
# sum(cdph_race_ethnicity$population_percent[cdph_race_ethnicity$date=="2020-10-10"])

agg_age_df <- read.csv("../CA_age_distn.csv",stringsAsFactors = F)
cdph_race_ethnicity$population <- cdph_race_ethnicity$population_percent*agg_age_df$population[match(cdph_race_ethnicity$age,agg_age_df$age_grp)]
cdph_race_ethnicity$time <- as.numeric(cdph_race_ethnicity$date - min(cdph_race_ethnicity$date)) 

races <- unique(cdph_race_ethnicity$race)
age_grps <- unique(cdph_race_ethnicity$age)
cdph_race_ethnicity$new_cases <- NA
for (i in 1:length(races)){
  for (j in 1:length(age_grps)){
    idx <- (cdph_race_ethnicity$race==races[i] & cdph_race_ethnicity$age==age_grps[j])
    cdph_race_ethnicity$new_cases[idx] <- calc_new_cases(cdph_race_ethnicity$confirmed_cases_total[idx])
  }
}

inc_df <- cdph_race_ethnicity[cdph_race_ethnicity$date==max(cdph_race_ethnicity$date) & cdph_race_ethnicity$race!="cdph-other",]
inc_df$exp_time <- NA
for (i in 1:nrow(inc_df)){
    idx <- (cdph_race_ethnicity$race==inc_df$race[i] & cdph_race_ethnicity$age==inc_df$age[i])
    inc_df$exp_time[i] <- calc_exposure_time(cdph_race_ethnicity[idx,])
}
summary(inc_df)
tapply(inc_df$confirmed_cases_total,inc_df$race,function(x){c(mean(x),var(x))})

lm <- glm(confirmed_cases_total ~ offset(log(exp_time)) + race + age,family = poisson(link=log),data = inc_df)
summary(lm)

# Testing
cdph_pstve_test_rate <- read.csv("cdph-positive-test-rate.csv",stringsAsFactors = F)

# Hospitalizations
cdph_hospital_patient_county <- read.csv("cdph-hospital-patient-county-totals.csv",stringsAsFactors = F)

# Prison counts (by prison and total)
cdcr_prison_totals <- read.csv("cdcr_prison_totals.csv",stringsAsFactors = F)
cdcr_state_totals <- read.csv("cdcr_prison_totals.csv",stringsAsFactors = F)

# Nursing home totals by county
cdph_nursing_home_county_total <- read.csv("cdph-nursing-home-county-totals.csv",stringsAsFactors = F)

# Adult and senior care counts (by facility and total)
cdph_ASC_facility <- read.csv("cdph-adult-and-senior-care-facilities.csv",stringsAsFactors = F)
cdph_ASC_total <- read.csv("cdph-adult-and-senior-care-totals.csv",stringsAsFactors = F)

# SNF counts (by facility and total)
cdph_SNF_facility <- read.csv("cdph-skilled-nursing-facilities.csv",stringsAsFactors = F)
cdph_SNF_totals <- read.csv("cdph-skilled-nursing-totals.csv",stringsAsFactors = F)





