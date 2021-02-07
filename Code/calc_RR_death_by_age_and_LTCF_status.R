rm(list=ls())

library(reshape2)

# Read in simulated population
df <- readRDS("../Data/CA_pop.RDS")

# Add variable for aggregation
df$population <- 1
df$SNF <- ifelse(df$special.population==3,1,0)

age_cat1_lbs <- c(seq(0,100,by=5),max(df$age)+1)
lbls <- paste0(age_cat1_lbs[1:(length(age_cat1_lbs)-1)],"-",age_cat1_lbs[1:(length(age_cat1_lbs)-1)]+4)
df$age_cat1 <- cut(df$age,age_cat1_lbs,right = F,labels = lbls)

# Aggregate over age and special population
agg_age_SNF <- aggregate(population ~ age_cat1 + SNF,df,sum)
# Cast to wide format
agg_age_SNF_wide <- dcast(agg_age_SNF,age_cat1 ~ SNF)
names(agg_age_SNF_wide)[names(agg_age_SNF_wide) %in% c(0,1)] <- c("non_SNF","SNF")
agg_age_SNF_wide$SNF[is.na(agg_age_SNF_wide$SNF)] <- 0

# Read in CDPH age death data
deaths_age <- read.csv("../Data/CDPHDeathAgeData/CDPHDeathAgeData.csv",stringsAsFactors = F)

deaths_age$YEARS[deaths_age$YEARS==">100"] <- "100-105"

# Extract lower and upper bounds of age groups
deaths_age$age_low <- as.numeric(sub("-.*","",deaths_age$YEARS))
deaths_age$age_upp <- as.numeric(sub(".*-","",deaths_age$YEARS)) - 1
# Correct age group labels
deaths_age$age_cat1 <- paste0(deaths_age$age_low,"-",deaths_age$age_low+4)

# Merge with population totals from simulated population
deaths_age_mrg <- merge(deaths_age,agg_age_SNF_wide,by="age_cat1",all.x=T)

# Calculate cumulative incidence of death for LTCF residents and general population
deaths_age_mrg$cum_inc_SNF <- deaths_age_mrg$SNF.LTC.residents/deaths_age_mrg$SNF
deaths_age_mrg$cum_inc_non_SNF <- deaths_age_mrg$Not.SNF.LTC.residents/deaths_age_mrg$non_SNF

calc_RR_CI <- function(RR,a,b,c,d,alpha = 0.05){
  z <- qnorm(1-alpha/2)
  CI_low <- exp(log(RR) - z * sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d)))
  CI_upp <- exp(log(RR) + z * sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d)))
  list(CI_low=CI_low,CI_upp=CI_upp)
}

# Calculate RR of death for LTCF residents compared to the general population
deaths_age_mrg$RR_SNF <- deaths_age_mrg$cum_inc_SNF/deaths_age_mrg$cum_inc_non_SNF
res <- calc_RR_CI(deaths_age_mrg$RR_SNF,deaths_age_mrg$SNF.LTC.residents,deaths_age_mrg$SNF,deaths_age_mrg$Not.SNF.LTC.residents,deaths_age_mrg$non_SNF)
deaths_age_mrg$RR_SNF_CI_low <- res$CI_low
deaths_age_mrg$RR_SNF_CI_upp <- res$CI_upp

# Calculate RR of death by age compared to 15-19-year-olds in LTCF residents and the general population
deaths_age_mrg$RR_age_SNF <- deaths_age_mrg$cum_inc_SNF/deaths_age_mrg$cum_inc_SNF[deaths_age_mrg$age_cat1=="20-24"]
res <- calc_RR_CI(deaths_age_mrg$RR_age_SNF,deaths_age_mrg$SNF.LTC.residents,deaths_age_mrg$SNF,deaths_age_mrg$SNF.LTC.residents[deaths_age_mrg$age_cat1=="20-24"],deaths_age_mrg$SNF[deaths_age_mrg$age_cat1=="20-24"])
deaths_age_mrg$RR_age_SNF_CI_low <- res$CI_low
deaths_age_mrg$RR_age_SNF_CI_upp <- res$CI_upp
deaths_age_mrg$RR_age_non_SNF <- deaths_age_mrg$cum_inc_non_SNF/deaths_age_mrg$cum_inc_non_SNF[deaths_age_mrg$age_cat1=="20-24"]
res <- calc_RR_CI(deaths_age_mrg$RR_age_non_SNF,deaths_age_mrg$Not.SNF.LTC.residents,deaths_age_mrg$non_SNF,deaths_age_mrg$Not.SNF.LTC.residents[deaths_age_mrg$age_cat1=="20-24"],deaths_age_mrg$non_SNF[deaths_age_mrg$age_cat1=="20-24"])
deaths_age_mrg$RR_age_non_SNF_CI_low <- res$CI_low
deaths_age_mrg$RR_age_non_SNF_CI_upp <- res$CI_upp

# Reorder
deaths_age_mrg <- deaths_age_mrg[c(2:nrow(deaths_age_mrg),1),]
write.csv(deaths_age_mrg,"../Data/CDPHDeathAgeData/RR_death_by_age_and_LTCF_status.csv",row.names = F)
