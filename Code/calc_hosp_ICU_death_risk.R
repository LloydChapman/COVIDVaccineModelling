calc_hosp_ICU_death_risk <- function(){

# age groups: 0-9, 10-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+
# sex: M/F

age_cats <- c("<10","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")
n_age_cats <- length(age_cats)
comorbs <- c("asthma","diabetes","smoker","heart.disease","heart.failure","hypertension","obesity")
n_comorbs <- length(comorbs)
comorbs1 <- as.vector(rbind(paste0(comorbs,0),paste0(comorbs,1)))
sexes <- c("female","male")
n_sexes <- length(sexes)

## Probability of clinical symptoms given infection
# Davies Nat Med 2020
p_clin_given_infctn <- c(0.29,0.21,0.27,0.33,0.4,0.49,0.63,0.69,0.69) # value for 70+ repeated for 80+
names(p_clin_given_infctn) <- age_cats

## Probability of hospitalisation given infection
# Salje Science 2020
p_hosp_given_infctn <- matrix(nrow=n_age_cats,ncol=n_sexes,dimnames = list(age_cats,sexes))
# females
p_hosp_given_infctn[,"female"] <- c(0.001,0.001,0.005,0.009,0.013,0.026,0.051,0.078,0.193) # value for 0-19 age group repeated for 0-9 and 10-19
# males
p_hosp_given_infctn[,"male"] <- c(0.002,0.002,0.006,0.012,0.016,0.032,0.067,0.11,0.376) # value for 0-19 age group repeated for 0-9 and 10-19

# # Verity LID 2020
# p_hosp_given_infctn <- c(0,0.0004,0.010,0.034,0.043,0.082,0.118,0.166,0.184) 

## Probability of hospitalisation given clinical symptoms
p_hosp_given_clin <- p_hosp_given_infctn/p_clin_given_infctn

# p_comorb_given_hosp <- 261/1099 # Guan NEJM 2020

## Probability of ICU given hospitalisation
p_ICU_given_hosp <- matrix(nrow=n_age_cats,ncol=n_sexes,dimnames = list(age_cats,sexes))

# # Salje Science 2020
# # females
# p_ICU_given_hosp[,"female"] <- c(0.167,0.167,0.087,0.119,0.166,0.207,0.231,0.187,0.042) # value for 0-19 age group repeated for 0-9 and 10-19
# # males
# p_ICU_given_hosp[,"male"] <- c(0.269,0.269,0.140,0.192,0.269,0.334,0.373,0.302,0.068) # value for 0-19 age group repeated for 0-9 and 10-19

# Lewnard BMJ 2020
# females
p_ICU_given_hosp[,"female"] <- c(0.331,0.331,0.311,0.207,0.323,0.346,0.339,0.332,0.307) # value for 10-19 age group used for 0-9 age group
# males
p_ICU_given_hosp[,"male"] <- c(0.373,0.373,0.351,0.433,0.525,0.489,0.464,0.571,0.436) # value for 10-19 age group used for 0-9 age group

## Probability of death given hospitalisation
p_death_given_hosp_by_age_and_sex <- matrix(nrow=n_age_cats,ncol=2,dimnames = list(age_cats,sexes))

# # Salje Science 2020
# # females
# p_death_given_hosp_by_age_and_sex[,"female"] <- c(0.005,0.005,0.009,0.015,0.026,0.052,0.101,0.167,0.252)
# # males
# p_death_given_hosp_by_age_and_sex[,"male"] <- c(0.007,0.007,0.013,0.022,0.038,0.076,0.148,0.246,0.371)

# Lewnard BMJ 2020
# females
p_death_given_hosp_by_age_and_sex[,"female"] <- c(0.017,0.017,0.026,0.036,0.062,0.099,0.151,0.238,0.365) # value for 10-19 age group used for 0-9 age group
# males
p_death_given_hosp_by_age_and_sex[,"male"] <- c(0.024,0.024,0.036,0.059,0.095,0.143,0.221,0.363,0.538) # value for 10-19 age group used for 0-9 age group

# Hazard ratios for co-morbidities in comorbs
# HR_ashtma <- (2454403 * 0.99  + 291670 * 1.13)/(2454403 + 291670) # population-weighted HR for asthma from Williamson Nature 2020, almost exactly 1 so just use 1
HR_DM <- (1038082 * 1.31 + 486491 * 1.95)/(1038082 + 486491) # population-weighted HR for diabetes from Williamson Nature 2020
HR <- c(1,HR_DM,1,1.16,1.77,1,1.33) # Docherty BMJ 2020 (some significant RFs from paper not included here, HRs for asthma and smoking taken as 1 as they were not significant in multivariable model, HR for heart failure taken from Petrilli BMJ 2020)
names(HR) <- comorbs

# Possibly include relative risk of death with hypertension by age from Grasselli JAMA 2020? Would be unadjusted for other co-morbidities though

rltve_probs <- numeric(2*n_comorbs)
rltve_probs[seq(1,2*n_comorbs-1,by=2)] <- 1/(1+HR)
rltve_probs[seq(2,2*n_comorbs,by=2)] <- HR/(1+HR)

p_death_given_hosp <- array(dim=c(n_age_cats,n_sexes,2*n_comorbs),dimnames = list(age_cats,sexes,comorbs1))
for (k in 1:length(rltve_probs)){
  p_death_given_hosp[,,k] <- p_death_given_hosp_by_age_and_sex*rltve_probs[k]
}

## Probability of death given ICU
p_death_given_ICU <- array(dim=c(n_age_cats,n_sexes,2*n_comorbs),dimnames = list(age_cats,sexes,comorbs1))
for (k in 1:(2*n_comorbs)){
  p_death_given_ICU[,,k] <- p_death_given_hosp[,,k]/p_ICU_given_hosp
}
  
# Calculate probabilities for each of the different end-points
p_ICU_given_clin <- p_ICU_given_hosp*p_hosp_given_clin
p_death_given_clin <- array(dim=c(n_age_cats,n_sexes,2*n_comorbs),dimnames = list(age_cats,sexes,comorbs1))
for (k in 1:(2*n_comorbs)){
  p_death_given_clin[,,k] <- p_death_given_hosp[,,k]*p_hosp_given_clin
}

# Calculate probabilities of death given clinical infection by age and sex
p_death_given_clin_by_age_and_sex <- p_death_given_hosp_by_age_and_sex*p_hosp_given_clin

return(list(p_clin_given_infctn=p_clin_given_infctn,p_hosp_given_clin=p_hosp_given_clin,p_ICU_given_clin=p_ICU_given_clin,p_death_given_clin=p_death_given_clin,p_death_given_hosp_by_age_and_sex=p_death_given_hosp_by_age_and_sex,p_death_given_clin_by_age_and_sex=p_death_given_clin_by_age_and_sex,HR=HR))

}

