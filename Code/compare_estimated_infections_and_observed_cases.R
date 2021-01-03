rm(list=ls())

library(ggplot2)
library(broom)

source("calc_hosp_ICU_death_risk.R")
source("analysis_functions.R")
source("prediction_functions.R")

# Read in aggregated death data
agg_deaths <- read.csv("~/Downloads/processed_death_data_2020-01-19_1.csv",stringsAsFactors = F)
agg_deaths$date_of_death <- as.Date(agg_deaths$date_of_death)

# Read in and process data from O'Driscoll Nature 2020
IFR <- read.csv("../Data/IFR_by_age_ODriscoll.csv",stringsAsFactors = F)
IFR$Age_group[IFR$Age_group=="80+"] <- "80-100"
IFR$age_low <- as.numeric(sub("-.*","",IFR$Age_group))
IFR$age_upp <- as.numeric(sub(".*-","",IFR$Age_group))
age_cat_lbs <- seq(0,80,by=10)
lbls <- paste0(age_cat_lbs,"-",age_cat_lbs+9)
lbls[lbls=="0-9"] <- "<10"
lbls[lbls=="80-89"] <- "80+"
IFR$age_cat <- cut(IFR$age_low,c(age_cat_lbs,max(IFR$age_upp)+1),labels=lbls,right=F)

# Calculate means across the two 5-year age groups in each 10-year age group
agg_IFR <- aggregate(.~age_cat,IFR[,c(2:(ncol(IFR)-3),ncol(IFR))],mean)

# Reshape to long format for merging
grp <- c("male","female","mean")
IFR_long <-  reshape(agg_IFR,varying = list(paste0("Median_perc","_",grp),paste0("CI_95_LB","_",grp),paste0("CI_95_UB","_",grp)),v.names = c("median_perc","CI_95_LB","CI_95_UB"),idvar = "age_cat",drop = c(paste0("Sero_study_range_LB_",grp),paste0("Sero_study_range_UB_",grp)),direction = "long")
IFR_long$grp <- grp[IFR_long$time]
IFR_long$sex <- NA
IFR_long$sex[IFR_long$grp=="female"] <- 0
IFR_long$sex[IFR_long$grp=="male"] <- 1

# Save
saveRDS(IFR_long,"../Data/IFR_by_age_ODriscoll.RDS")

# Backcalculate infections in each risk group using O'Driscoll IFR
cum_deaths <- backcalculate_infections(agg_deaths,IFR_long)
# # Multiply estimated number of infections from O'Driscoll IFR by ratio of total estimated from CDC CA seroprevalence estimates to that from O'Driscoll, so that total numbers of infections agree
# cum_deaths$n <- cum_deaths$n * sum(CA_IFR$infections)/sum(cum_deaths$n)
print(sum(cum_deaths$n)) # 6001684
# Plot distribution of numbers infected in each risk group 
hist(cum_deaths$n)

# Calculate cases in each risk group using Davies age-dependent clinical fraction
p_clin <- calc_hosp_ICU_death_risk()$p_clin_given_infctn
p_clin <- data.frame(age_cat=names(p_clin),p_clin=p_clin)
cum_deaths <- calc_cases_from_infections(cum_deaths,p_clin)


# Read in aggregated case data
agg_cases_age <- read.csv("~/Downloads/agg_cases_age.csv",stringsAsFactors = F)

cum_inf_age <- aggregate(n ~ age_cat,cum_deaths,sum)
cum_cases_age <- aggregate(n ~ age_cat,agg_cases_age,sum)

cum_inf_cases_age <- merge(cum_inf_age,cum_cases_age,by="age_cat")
cum_inf_cases_age$prop.x <- cum_inf_cases_age$n.x/sum(cum_inf_cases_age$n.x)
p_clin <- calc_hosp_ICU_death_risk()$p_clin_given_infctn
cum_inf_cases_age$n.x_clin <- cum_inf_cases_age$n.x*p_clin
cum_inf_cases_age$prop.x_clin <- cum_inf_cases_age$n.x_clin/sum(cum_inf_cases_age$n.x_clin)
cum_inf_cases_age$prop.y <- cum_inf_cases_age$n.y/sum(cum_inf_cases_age$n.y)

obsvd_total_cases <- sum(cum_inf_cases_age$n.y) #865203
estd_total_cases <- sum(cum_inf_cases_age$n.x_clin) #2140575
IFR_ratio <- estd_total_cases/obsvd_total_cases #2.474072
saveRDS(IFR_ratio,"../Data/IFR_ratio.RDS")

# Plot observed case counts vs estimated case counts by age
cum_inf_cases_age_long <- melt(cum_inf_cases_age[,c("age_cat","n.y","n.x_clin")])
pdf("../Figures/Backcalculation/obsvd_vs_estd_case_age_counts.pdf",width = 6,height = 4)
ggplot(cum_inf_cases_age_long,aes(x=age_cat,y=value,group=variable,fill=variable)) + geom_bar(stat="identity",position="dodge") + xlab("Age") + ylab("Cases") + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_fill_discrete(name="",labels=c("obsvd","estd"))
dev.off()

# Plot observed case counts vs corrected estimated case counts by age
cum_inf_cases_age_long$value[cum_inf_cases_age_long$variable=="n.x_clin"] <- cum_inf_cases_age_long$value[cum_inf_cases_age_long$variable=="n.x_clin"]/IFR_ratio
pdf("../Figures/Backcalculation/obsvd_vs_estd_case_age_counts_crrctd.pdf",width = 6,height = 4)
ggplot(cum_inf_cases_age_long,aes(x=age_cat,y=value,group=variable,fill=variable)) + geom_bar(stat="identity",position="dodge") + xlab("Age") + ylab("Cases") + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_fill_discrete(name="",labels=c("obsvd","estd"))
dev.off()

# Plot observed and estimated age distributions of cases
cum_inf_cases_age_long1 <- melt(cum_inf_cases_age[,c("age_cat","prop.x","prop.y","prop.x_clin")])
pdf("../Figures/Backcalculation/obsvd_vs_estd_case_age_distns.pdf",width = 6,height = 4)
ggplot(cum_inf_cases_age_long1,aes(x=age_cat,y=value,group=variable,fill=variable)) + geom_bar(stat="identity",position="dodge") + xlab("Age") + ylab("Proportion") + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_fill_discrete(name="",labels=c("estd infctns","obsvd cases","estd cases"))
dev.off()

# Read in CDC seroprevalence data
seroprev_long <- read.csv("../Data/CA_seroprevalence.csv",stringsAsFactors = F)
seroprev <- aggregate(seroprev ~ group,seroprev_long,mean)
seroprev <- seroprev[!(seroprev$group %in% c("Cumulative","Female","Male")),]
names(seroprev)[names(seroprev)=="group"] <- "age_cat_sero"
setDT(seroprev)

df <- readRDS("../Data/CA_pop_with_risk_ests_deaths.RDS")
agg_df_age <- aggregate(seroprev ~ age_cat_sero,df,sum)
agg_df_age$age_low <- as.numeric(sub("-.*|\\+","",agg_df_age$age_cat_sero))
agg_df_age$age_upp <- as.numeric(sub(".*-|\\+","",agg_df_age$age_cat_sero))
agg_df_age$age_upp[nrow(agg_df_age)] <- 100
agg_df_age$age_mid <- (agg_df_age$age_low + agg_df_age$age_upp)/2

cum_inf_cases_age$age_low <- as.numeric(sub("<|-.*|\\+","",cum_inf_cases_age$age_cat))
cum_inf_cases_age$age_low[1] <- 0
cum_inf_cases_age$age_upp <- c(cum_inf_cases_age$age_low[1:(nrow(cum_inf_cases_age)-1)]+9,100)
cum_inf_cases_age$age_mid <- (cum_inf_cases_age$age_low + cum_inf_cases_age$age_upp)/2

print(sum(agg_df_age$seroprev)) #2047013
print(sum(cum_inf_cases_age$n.x)) #6001684
print(sum(cum_inf_cases_age$n.x)/sum(agg_df_age$seroprev)) #2.931922

agg_df_age$n.x <- NA
agg_df_age$n.x[1] <- cum_inf_cases_age$n.x[cum_inf_cases_age$age_cat=="<10"] + 0.8*cum_inf_cases_age$n.x[cum_inf_cases_age$age_cat=="10-19"]
agg_df_age$n.x[2] <- 0.2*cum_inf_cases_age$n.x[cum_inf_cases_age$age_cat=="10-19"] + cum_inf_cases_age$n.x[cum_inf_cases_age$age_cat=="20-29"] + cum_inf_cases_age$n.x[cum_inf_cases_age$age_cat=="30-39"] + cum_inf_cases_age$n.x[cum_inf_cases_age$age_cat=="40-49"]
agg_df_age$n.x[3] <- cum_inf_cases_age$n.x[cum_inf_cases_age$age_cat=="50-59"] + 0.5*cum_inf_cases_age$n.x[cum_inf_cases_age$age_cat=="60-69"]
agg_df_age$n.x[4] <- 0.5*cum_inf_cases_age$n.x[cum_inf_cases_age$age_cat=="60-69"] + cum_inf_cases_age$n.x[cum_inf_cases_age$age_cat=="70-79"] + cum_inf_cases_age$n.x[cum_inf_cases_age$age_cat=="80+"]

# Plot observed vs estimated numbers of infections by age
pdf("../Figures/Backcalculation/infections_seroprev_vs_IFR.pdf",width = 6,height = 4)
ggplot() + geom_line(aes(x=age_mid,y=seroprev,color="seroprev"),agg_df_age) + geom_line(aes(x=age_mid,y=n.x,color="IFR"),data=agg_df_age) + xlab("Age") + ylab("Infections") + scale_color_manual(name="",values=c("seroprev"="black","IFR"="red"))
dev.off()

# Plot observed vs corrected estimated numbers of infections by age
pdf("../Figures/Backcalculation/infections_seroprev_vs_IFR_crrctd.pdf",width = 6,height = 4)
ggplot() + geom_line(aes(x=age_mid,y=seroprev,color="seroprev"),agg_df_age) + geom_line(aes(x=age_mid,y=n.x/IFR_ratio,color="IFR"),data=agg_df_age) + xlab("Age") + ylab("Infections") + scale_color_manual(name="",values=c("seroprev"="black","IFR"="red"))
dev.off()

# Apply IFR correction to estimated infections and cases
cum_deaths$n <- cum_deaths$n/IFR_ratio
cum_deaths$n_cases <- cum_deaths$n_cases/IFR_ratio

# Load population data
agg_pop <- read.csv("../Data/agg_pop1.csv",stringsAsFactors = F)
# Merge with infections and deaths data frame
inc <- merge(cum_deaths,agg_pop,by=c("county_res","age_cat","sex","race_ethnicity"),all.x = T)
inc$exp_time <- inc$population * as.numeric(as.Date("2020-10-19")-as.Date("2020-01-19"))/1e5

# Fit Poisson regression model to estimated infections
glm_fit <- glm(round(n) ~ offset(log(exp_time)) + county_res + age_cat + sex + race_ethnicity,family = poisson(link = log),data = inc)
saveRDS(glm_fit,"../Data/regression_output_death_IFR_model_infections.RDS")
coeffs <- exp(coef(glm_fit))
print(coeffs)

# process_regression_output("../Data/","regression_output_death_IFR_model_infections.RDS")

# Fit Poisson regression model to estimated infections
glm_fit3 <- glm(round(n_cases) ~ offset(log(exp_time)) + county_res + age_cat + sex + race_ethnicity,family = poisson(link = log),data = inc)
saveRDS(glm_fit3,"../Data/regression_output_death_IFR_model.RDS")
coeffs3 <- exp(coef(glm_fit3))
print(coeffs3)

process_regression_output("../Data/","regression_output_death_IFR_model.RDS")

# Load coefficient estimates from 3-month case model
glm_fit1 <- readRDS("../Data/regression_output_2020-07-19_1.RDS")
coeffs1 <- exp(coef(glm_fit1))

# Fit Poisson regression model to cumulative cases
inc1 <- glm_fit1$data
inc1$exp_time <- inc1$population * inc1$time[1]/1e5
glm_fit2 <- glm(cum_cases ~ offset(log(exp_time)) + county_res + age_cat + sex + race_ethnicity,poisson(link = log),data = inc1)
coeffs2 <- exp(coef(glm_fit2))
print(coeffs2)

normalise_coeffs <- function(x,coeffs){
  x[grep("age",names(x))] <- x[grep("age",names(x))] * x["(Intercept)"]/coeffs["(Intercept)"]
  x["(Intercept)"] <- x["(Intercept)"] * coeffs["(Intercept)"]/x["(Intercept)"]
  return(x)
}

# Normalise RRs for different age categories for regression of estimated cases under death IFR model by baseline risk
coeffs3 <- normalise_coeffs(coeffs3,coeffs2)

# Combine coefficient estimates into single data frame
coeffs_comp <- data.frame(names = names(coeffs),death_IFR_model_inf=coeffs,death_IFR_model_cases=coeffs3,three_month_PH_model=coeffs1,three_month_Poiss_model=coeffs2)
coeffs_comp_long <- melt(coeffs_comp)

# Plot coefficients to compare them - [ ] RERUN REGRESSIONS WITH 50-59 AS REFERENCE AGE CATEGORY
pdf("../Figures/model_coefficients_comparison.pdf",width = 9,height = 5)
ggplot(coeffs_comp_long,aes(x=names,y=value,group=variable,color=as.factor(variable))) + geom_point() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + xlab("coefficient") + labs(color="Model")
dev.off()
