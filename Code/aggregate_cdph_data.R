rm(list=ls())

library(lubridate)

setwd("/mnt/nlo/cdph_lane6_dua_tables_20210102/")

y <- read.csv("tbl1_case_linelist.csv")
y$sex <- NA
y$sex[y$gender=="F"] <- 0
y$sex[y$gender=="M"] <- 1

y$cum_cases <- 1
y$date_of_death <- as.Date(y$date_of_death)
y$cum_deaths <- 0
y$cum_deaths[!is.na(y$date_of_death) & year(y$date_of_death)>=2020] <- 1

# Read in aggregated simulated population
agg_pop <- read.csv("/mnt/nlo_shared/data/agg_pop2.csv",stringsAsFactors = F)

# Function for plotting numbers and cumulative incidence
aggregate_by_demogrphcs <- function(y,df,demogrphcs,vrble){
  lagg_y <- vector("list",length(demogrphcs))
  lagg_mrg <- vector("list",length(demogrphcs))
  for (i in 1:length(demogrphcs)){
    agg_y <- aggregate(y[,vrble],by=list(y[,demogrphcs[i]]),FUN=sum)
    names(agg_y)[names(agg_y)=="Group.1"] <- demogrphcs[i]
    agg_df <- aggregate(df$population,by=list(df[,demogrphcs[i]]),FUN=sum)
    names(agg_df)[names(agg_df)=="Group.1"] <- demogrphcs[i]
    names(agg_df)[names(agg_df)=="x"] <- "population"
    agg_mrg <- merge(agg_df,agg_y,by=demogrphcs[i],all.x=T)
    agg_mrg$cum_inc <- 100*agg_mrg$x/agg_mrg$population
    
    lagg_y[[i]] <- agg_y
    lagg_mrg[[i]] <- agg_mrg
  }
  return(list(lagg_y=lagg_y,lagg_mrg=lagg_mrg))
}

# Aggregate cases and deaths by demographic factors
demogrphcs <- c("county_res","age_cat","sex","race_ethnicity")
l <- aggregate_by_demogrphcs(y,agg_pop,demogrphcs,"cum_cases")
l_deaths <- aggregate_by_demogrphcs(y,agg_pop,demogrphcs,"cum_deaths")

# Save
save(l,l_deaths,file="/mnt/nlo_shared/data/case_and_death_demogrphcs.RData")




