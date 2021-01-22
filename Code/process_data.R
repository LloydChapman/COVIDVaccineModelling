rm(list=ls())

source("processing_functions.R")

# Change county names to those in CDPH data
county_names <- read.csv("../Data/county_names.csv",stringsAsFactors = F)

# Load synthetic CA population
df <- load_synth_pop("../Data/Simulation1226/",county_names)

# Convert to data table
setDT(df)

# Save
saveRDS(df,"../Data/CA_pop1.RDS")

### Aggregation for working with CDPH data ###
# Count variable for aggregation
df$population <- 1

# Aggregate by age and sex
agg_age_sex_df <- aggregate(population ~ sex + age_cat,df,sum)
agg_age_sex_df <- dcast(agg_age_sex_df,age_cat ~ sex)
write.csv(agg_age_sex_df,"../Data/CA_age_sex_distn2.csv",row.names = F)

# Aggregate with sub-totals for each special population in each demographic risk group
df$non_spec_pop <- ifelse(df$special.population==0,1,0)
df$hcw <- ifelse(df$special.population==1,1,0)
df$prisoner <- ifelse(df$special.population==2,1,0)
df$snf <- ifelse(df$special.population==3,1,0)
df$educator <- ifelse(df$special.population==4,1,0)
df$homeless <- ifelse(df$special.population==5,1,0)
df$essential_f <- ifelse(df$special.population==6,1,0)
df$essential_nf <- ifelse(df$special.population==7,1,0)
df$alf <- ifelse(df$special.population==8,1,0)
  
# Aggregate (without special populations)
agg_df <- aggregate(cbind(population,non_spec_pop,hcw,prisoner,snf,teacher,homeless,essential_f,essential_nf,alf) ~ county_res + age_cat + sex + race_ethnicity,df,sum)
cols <- c("population","non_spec_pop","hcw","prisoner","snf","educator","homeless","essential_f","essential_nf","alf")
agg_df <- df[,lapply(.SD,sum),.SDcols=cols,by=.(county_res,age_cat,sex,race_ethnicity)]
agg_df <- agg_df[do.call(order,agg_df),]
View(agg_df)

# Look at subgroup sizes
summary(agg_df$population)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4    1724    6991   16999   17968  421616 
hist(agg_df$population,100)

# Save
write.csv(agg_df,"../Data/agg_pop2.csv",row.names = F)

