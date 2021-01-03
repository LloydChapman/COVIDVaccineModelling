rm(list=ls())

source("~/Dropbox/COVIDVaccineModelling/Code/processing_functions.R")

# Set working directory
dir <- "~/Dropbox/COVIDVaccineModelling/Data/"
setwd(dir)

# Change county names to those in CDPH data
county_names <- read.csv("county_names.csv",stringsAsFactors = F)

# Load synthetic CA population
df <- load_synth_pop(paste0(dir,"Simulation1126/"),county_names)

# Add essential workers to special populations
set.seed(123)
df <- add_essential_workers(df)

# Save
saveRDS(df,"CA_pop.RDS")

### Aggregation for working with CDPH data ###
# Count variable for aggregation
df$population <- 1

# Aggregate by age and sex
agg_age_sex_df <- aggregate(population ~ sex + age_cat,df,sum)
agg_age_sex_df <- dcast(agg_age_sex_df,age_cat ~ sex)
write.csv(agg_age_sex_df,"CA_age_sex_distn1.csv",row.names = F)

# Aggregate with sub-totals for each special population in each demographic risk group
df$non_spec_pop <- ifelse(df$special.population==0,1,0)
df$hcw <- ifelse(df$special.population==1,1,0)
df$prisoner <- ifelse(df$special.population==2,1,0)
df$snf <- ifelse(df$special.population==3,1,0)
df$teacher <- ifelse(df$special.population==4,1,0)
df$homeless <- ifelse(df$special.population==5,1,0)
df$essential <- ifelse(df$special.population==6,1,0)
  
# Aggregate (without special populations)
# agg_df <- aggregate(cbind(population,non_spec_pop,hcw,prisoner,snf,teacher,homeless) ~ county_res + age_cat + sex + race_ethnicity,df,sum)
agg_df <- aggregate(cbind(population,non_spec_pop,hcw,prisoner,snf,teacher,homeless,essential) ~ county_res + age_cat + sex + race_ethnicity,df,sum)
agg_df <- agg_df[do.call(order,agg_df),]
View(agg_df)

# Look at subgroup sizes
summary(agg_df$population)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4    1712    6969   16999   17957  421616 
hist(agg_df$population,100)

# Save
write.csv(agg_df,"agg_pop1.csv",row.names = F)

