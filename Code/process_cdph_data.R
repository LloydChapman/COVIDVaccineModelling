rm(list=ls())

library(lubridate)

source("/mnt/nlo_shared/code/processing_functions.R")
# source("~/Dropbox/COVIDVaccineModelling/Code/processing_functions.R")

setwd("/mnt/nlo/cdph_lane6_dua_tables_20201019/")

#################
### Case data ###
#################

# Read in data
cases <- read.csv("tbl1_case_linelist.csv",stringsAsFactors=F)
print(dim(cases))
print(head(cases))
print(summary(cases))
print(summary(as.factor(cases$race_ethnicity)))

# Convert date strings to dates
date_cols <- grep("date",names(cases))
for (i in 1:length(date_cols)){
  cases[,date_cols[i]] <- as.Date(cases[,date_cols[i]])
}
print(sum(is.na(cases$onset_date)))
print(sum(is.na(cases$lab_result_date)))
print(sum(is.na(cases$first_report_date)))
print(sum(is.na(cases$episode_date)))
print(cases[is.na(cases$first_report_date),])
print(summary(as.numeric(cases$episode_date-cases$first_report_date)))

# Correct incorrect date of death
print(cases[cases$date_of_death=="1975-06-06" & !is.na(cases$date_of_death),])
print(head(sort(cases$date_of_death[!is.na(cases$date_of_death)])))
cases$date_of_death[cases$date_of_death=="1975-06-06"] <- NA
cases$episode_date[cases$episode_date=="1975-06-06"] <- NA                                 

# Correct lab result dates
print(head(sort(cases$lab_result_date[!is.na(cases$lab_result_date)])))
print(tail(sort(cases$lab_result_date[!is.na(cases$lab_result_date)])))
print(head(cases[order(cases$lab_result_date),]))
print(sort(cases$lab_result_date[!is.na(cases$lab_result_date)])[1:40])

# Remove completely incorrect lab result dates
cases$lab_result_date[cases$lab_result_date=="1899-12-30"] <- NA
cases$lab_result_date[cases$lab_result_date=="1999-12-31"] <- NA  

# Correct years of lab result dates with wrong year
cases$lab_result_year <- as.numeric(format(cases$lab_result_date,"%Y"))
print(cases$lab_result_year[cases$lab_result_year<2020 & !is.na(cases$lab_result_year)])
cases$lab_result_year[cases$lab_result_year<2020 & !is.na(cases$lab_result_year)] <- 2020
print(summary(cases$lab_result_year))
cases$lab_result_date <- as.Date(paste0(cases$lab_result_year,"-",format(cases$lab_result_date,"%m-%d")))
print(summary(cases$lab_result_date))

# Make numeric age-group variable
cases$age_grp <- match(cases$age_cat,unique(cases$age_cat))

# Create binary sex variable based on gender variable
print(summary(as.factor(cases$gender)))
cases$sex <- NA
cases$sex[cases$gender=="F"] <- 0
cases$sex[cases$gender=="M"] <- 1
print(summary(as.factor(cases$sex)))

# Count cases with all demographic variables and reporting date
print(sum(cases$county_res!="UNASSIGNED" & !is.na(cases$age_cat) & !is.na(cases$sex) & cases$race_ethnicity!="Unknown" & !is.na(cases$first_report_date)))
# 606469

# Add count variable for aggregation
cases$n <- 1

# Add count variable for aggregation of deaths
cases$n_deaths <- 0
cases$n_deaths[!is.na(cases$date_of_death)] <- 1

# # Save
# write.csv(cases,"/mnt/nlo_shared/data/cases.csv",row.names = F) 

####################
### Testing data ###
####################

tests <- read.csv("tbl3_testing_linelist.csv",stringsAsFactors=F)
print(dim(tests))
print(summary(tests))
print(summary(as.factor(tests$result)))

# Deduplicated positive test data
tests1 <- read.csv("tbl3_a_testing_linelist_pos_only.csv",stringsAsFactors=F)                    
print(dim(tests1))
print(summary(tests1))

tests$n_tests <- 1
agg_tests_age <- aggregate(n_tests ~ result_date + age_cat + result,tests,sum)
write.csv(agg_tests_age,"/mnt/nlo_shared/data/agg_tests_age.csv",row.names = F)

###########################
### Aggregate case data ###
###########################

# Aggregate
agg_cases <- aggregate(cbind(n,n_deaths) ~ first_report_date + county_res + age_cat + sex + race_ethnicity,cases,sum)
agg_cases <- agg_cases[do.call(order,agg_cases),]
print(head(agg_cases))

agg_cases_age <- aggregate(cbind(n,n_deaths) ~ first_report_date + age_cat,cases,sum)
write.csv(agg_cases_age,"/mnt/nlo_shared/data/agg_cases_age.csv",row.names = F)

# # Create exposure time variable
# # Find earliest event date
# apply(cases[,c("lab_result_date","first_report_date")],2,function(x) min(x,na.rm = T))
# start_date <- as.Date(min(apply(cases[,c("lab_result_date","first_report_date")],2,function(x) min(x,na.rm = T))))
# agg_cases$time <- as.numeric(agg_cases$first_report_date - start_date)
# print(head(agg_cases))
# # Save
# write.csv(agg_cases,"/mnt/nlo_shared/data/agg_cases.csv",row.names = F)
# 
# ## Subset data to last 3 months ##
# start_date1 <- as.Date("2020-07-15")
# cases1 <- cases[cases$first_report_date > start_date1,]
# cases1 <- cases[cases$first_report_date >= start_date1 & !is.na(cases$first_report_date),]
# print(dim(cases1))
# # Save
# write.csv(cases1,"/mnt/nlo_shared/data/cases1.csv",row.names = F)
# 
# # Aggregate
# agg_cases1 <- aggregate(n ~ first_report_date + county_res + age_cat + sex + race_ethnicity,cases1,sum)
# agg_cases1 <- agg_cases1[do.call(order,agg_cases1),]
# print(dim(agg_cases1))
# agg_cases1$time <- as.numeric(agg_cases1$first_report_date - start_date1)
# print(head(agg_cases1))
# 
# # Calculate cumulative cases by date for each risk factor group
# agg_cases1$cum_cases <- ave(agg_cases1$n,list(agg_cases1$county_res,agg_cases1$age_cat,agg_cases1$sex,agg_cases1$race_ethnicity),FUN=cumsum)
# 
# # Check
# print(agg_cases1[order(agg_cases1$race_ethnicity,agg_cases1$sex,agg_cases1$age_cat,agg_cases1$county_res,agg_cases1$first_report_date),][1:40,])
# 
# # Save
# write.csv(agg_cases1,"/mnt/nlo_shared/data/agg_cases1.csv",row.names = F)

############################
### Aggregate death data ###
############################

# Aggregate
agg_deaths <- aggregate(n_deaths ~ date_of_death + county_res + age_cat + sex + race_ethnicity,cases,sum)
agg_deaths <- agg_deaths[do.call(order,agg_deaths),]
print(head(agg_deaths))

# # Save
# write.csv(agg_deaths,"/mnt/nlo_shared/data/agg_deaths1.csv",row.names = F)

####################################
### Calculate analysis variables ###
####################################

# Create dummy data
agg_cases_dummy <- agg_cases
set.seed(123)
agg_cases_dummy$n <- rpois(nrow(agg_cases),mean(agg_cases$n))
write.csv(agg_cases_dummy,"/mnt/nlo_shared/data/agg_cases_dummy1.csv",row.names=F)

agg_deaths_dummy <- agg_deaths
agg_deaths_dummy$n_deaths <- rpois(nrow(agg_deaths),mean(agg_deaths$n_deaths))
write.csv(agg_deaths_dummy,"/mnt/nlo_shared/data/agg_deaths_dummy1.csv",row.names=F)

# # Read in dummy data
# setwd("~/Dropbox/COVIDVaccineModelling/Code")
# agg_cases <- read.csv("../Data/agg_cases_dummy.csv")
# agg_cases$first_report_date <- as.Date(agg_cases$first_report_date)                                            

# Subset and process the data for the last 1, 3, 6 and 9 months 
start_dates <- max(agg_cases$first_report_date) %m-% months(c(1,3,6,9))
agg_pop <- read.csv("/mnt/nlo_shared/data/agg_pop1.csv",stringsAsFactors=F)
for (i in 1:length(start_dates)){
  y <- calc_analysis_vars(agg_cases,start_dates[i],agg_pop)
  # Save
  write.csv(y,paste0("/mnt/nlo_shared/data/processed_data_",start_dates[i],"_1.csv"),row.names = F)
  # write.csv(y,paste0("../Data/processed_data_",start_date,"_1.csv"),row.names = F)
  y_deaths <- calc_death_analysis_vars(agg_deaths,start_dates[i],agg_pop)
  # Save
  write.csv(y_deaths,paste0("/mnt/nlo_shared/data/processed_death_data_",start_dates[i],"_1.csv"),row.names = F)
}
