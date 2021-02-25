# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Author: Poojan Shukla
# Goal: Simulate populations of all counties in CA

rm(list=ls())

#Load Packages
require(ggplot2)
library(triangle)
library(dplyr)
library(tidycensus) # the one you care about 
library(acs)
library(censusapi)
library(devtools)
if (!require(hashmap)){
  install_github("nathan-russell/hashmap")
}
library(hashmap)
library(foreach)
library(doParallel)
library(MASS)
library(readxl)
library(truncnorm)

chis_data_validation <- function(county_name, county_df) {
  #for each co morbidiity, examine prevalence by age
  #now check by race (this should be flat because of the way its programmed)
  #now by eth
  #now by sex
  #county_df = Yuba_df
  #county_name = "Yuba County, California"

  pop = nrow(county_df)
  total_ab40 = sum(county_df$Asthma)
  total_ab34 = sum(county_df$Heart.Disease)
  total_diabetes = sum(county_df$Diabetes)
  total_smoking = sum(county_df$Smoker)
  total_ab29 = sum(county_df$Hypertension)
  total_ab52 = sum(county_df$Heart.Failure)
  avg_bmi = sum(county_df$BMI)
  
  filt_ab34_age_county = filter(chis_age_county_ab34, srcnty_code == county_name)
  filt_ab40_age_county = filter(chis_age_county_ab40, srcnty_code == county_name)
  filt_ab29_age_county = filter(chis_age_county_ab29, srcnty_code == county_name)
  filt_diabetes_age_county = filter(chis_age_county_diabetes, srcnty_code == county_name)
  filt_smkcur_age_county = filter(chis_age_county_smkcur, srcnty_code == county_name)
  filt_bmi_age_county = filter(chis_age_county_bmi, srcnty_code == county_name)
  filt_bmi_age_county_totals = filter(chis_age_county_bmi_totals, srcnty_code == county_name)
  
  filt_ab34_age_county[filt_ab34_age_county == 'na'] = 0
  filt_ab34_age_county[is.na(filt_ab34_age_county)] = 0
  
  filt_ab40_age_county[filt_ab40_age_county == 'na'] = 0
  filt_ab40_age_county[is.na(filt_ab40_age_county)] = 0
  
  filt_ab29_age_county[filt_ab29_age_county == 'na'] = 0
  filt_ab29_age_county[is.na(filt_ab29_age_county)] = 0
  
  filt_diabetes_age_county[filt_diabetes_age_county == 'na'] = 0
  filt_diabetes_age_county[is.na(filt_diabetes_age_county)] = 0
  
  filt_smkcur_age_county[filt_smkcur_age_county == 'na'] = 0
  filt_smkcur_age_county[is.na(filt_smkcur_age_county)] = 0
  
  filt_bmi_age_county[filt_bmi_age_county == 'na'] = 0
  filt_bmi_age_county[is.na(filt_bmi_age_county)] = 0
  
  filt_bmi_age_county_totals[filt_bmi_age_county_totals == 'na'] = 0
  filt_bmi_age_county_totals[is.na(filt_bmi_age_county_totals)] = 0
  
  print(paste("Ab40 prevalence", total_ab40/pop))
  print(paste("Ab34 prevalence", total_ab34/pop))
  print(paste("Ab29 prevalence", total_ab29/pop))
  print(paste("Smoking prevalence", total_smoking/pop))
  print(paste("Diabetes prevalence", total_diabetes/pop))
  print(paste("Ab52 prevalence", total_ab52/pop))
  print(paste("Avg BMI", avg_bmi/pop))
  
  CA_wide = c(ab40_overall, ab34_overall, ab29_overall,
              ab52_overall, diabetes_overall, smkcur_overall, CA_bmi_mean)
  county_overall = c(total_ab40/pop, total_ab34/pop, total_ab29/pop,
                     total_ab52/pop, total_diabetes/pop, total_smoking/pop, avg_bmi)
  comparison_df = data.frame(county_overall, CA_wide)
  colnames(comparison_df) = c("Simulation", "Real Data")
  rownames = c("Asthma-ab40", "Heart Disease-ab34", "HTN-ab29", "HF-ab52", "Diabetes",
               "Smoking", "BMI")
  
  pop_18_24 = nrow(filter(county_df, Age>=18, Age<25))
  pop_25_34 = nrow(filter(county_df, Age>=25, Age<35))
  pop_35_44 = nrow(filter(county_df, Age>=35, Age<45))
  pop_45_54 = nrow(filter(county_df, Age>=45, Age<55))
  pop_55_64 = nrow(filter(county_df, Age>=55, Age<65))
  pop_65_74 = nrow(filter(county_df, Age>=65, Age<75))
  pop_75_84 = nrow(filter(county_df, Age>=75, Age<85))
  pop_85 = nrow(filter(county_df, Age>=85))
  
  bin_18_24 = filter(county_df, Age>=18, Age<25)
  bin_25_34 = filter(county_df, Age>=25, Age<35)
  bin_35_44 = filter(county_df, Age>=35, Age<45)
  bin_45_54 = filter(county_df, Age>=45, Age<55)
  bin_55_64 = filter(county_df, Age>=55, Age<65)
  bin_65_74 = filter(county_df, Age>=65, Age<75)
  bin_75_84 = filter(county_df, Age>=75, Age<85)
  bin_85 = filter(county_df, Age>=85)
  
  ab40_prev_18_24 = sum(bin_18_24$Asthma)/pop_18_24
  ab34_prev_18_24 = sum(bin_18_24$Heart.Disease)/pop_18_24
  ab29_prev_18_24 = sum(bin_18_24$Hypertension)/pop_18_24
  ab52_prev_18_24 = sum(bin_18_24$Heart.Failure)/pop_18_24
  smoking_prev_18_24 = sum(bin_18_24$Smoker)/pop_18_24
  diabetes_prev_18_24 = sum(bin_18_24$Diabetes)/pop_18_24
  bmi_18_24 = sum(bin_18_24$BMI)/pop_18_24
  
  ab40_prev_25_34 = sum(bin_25_34$Asthma)/pop_25_34
  ab34_prev_25_34 = sum(bin_25_34$Heart.Disease)/pop_25_34
  ab29_prev_25_34 = sum(bin_25_34$Hypertension)/pop_25_34
  ab52_prev_25_34 = sum(bin_25_34$Heart.Failure)/pop_25_34
  smoking_prev_25_34 = sum(bin_25_34$Smoker)/pop_25_34
  diabetes_prev_25_34 = sum(bin_25_34$Diabetes)/pop_25_34
  bmi_25_34 = sum(bin_25_34$BMI)/pop_25_34
  
  ab40_prev_35_44 = sum(bin_35_44$Asthma)/pop_35_44
  ab34_prev_35_44 = sum(bin_35_44$Heart.Disease)/pop_35_44
  ab29_prev_35_44 = sum(bin_35_44$Hypertension)/pop_35_44
  ab52_prev_35_44 = sum(bin_35_44$Heart.Failure)/pop_35_44
  smoking_prev_35_44 = sum(bin_35_44$Smoker)/pop_35_44
  diabetes_prev_35_44 = sum(bin_35_44$Diabetes)/pop_35_44
  bmi_35_44 = sum(bin_35_44$BMI)/pop_35_44
  
  ab40_prev_45_54 = sum(bin_45_54$Asthma)/pop_45_54
  ab34_prev_45_54 = sum(bin_45_54$Heart.Disease)/pop_45_54
  ab29_prev_45_54 = sum(bin_45_54$Hypertension)/pop_45_54
  ab52_prev_45_54 = sum(bin_45_54$Heart.Failure)/pop_45_54
  smoking_prev_45_54 = sum(bin_45_54$Smoker)/pop_45_54
  diabetes_prev_45_54 = sum(bin_45_54$Diabetes)/pop_45_54
  bmi_45_54 = sum(bin_45_54$BMI)/pop_45_54
  
  ab40_prev_55_64 = sum(bin_55_64$Asthma)/pop_55_64
  ab34_prev_55_64 = sum(bin_55_64$Heart.Disease)/pop_55_64
  ab29_prev_55_64 = sum(bin_55_64$Hypertension)/pop_55_64
  ab52_prev_55_64 = sum(bin_55_64$Heart.Failure)/pop_55_64
  smoking_prev_55_64 = sum(bin_55_64$Smoker)/pop_55_64
  diabetes_prev_55_64 = sum(bin_55_64$Diabetes)/pop_55_64
  bmi_55_64 = sum(bin_55_64$BMI)/pop_55_64
  
  ab40_prev_65_74 = sum(bin_65_74$Asthma)/pop_65_74
  ab34_prev_65_74 = sum(bin_65_74$Heart.Disease)/pop_65_74
  ab29_prev_65_74 = sum(bin_65_74$Hypertension)/pop_65_74
  ab52_prev_65_74 = sum(bin_65_74$Heart.Failure)/pop_65_74
  smoking_prev_65_74 = sum(bin_65_74$Smoker)/pop_65_74
  diabetes_prev_65_74 = sum(bin_65_74$Diabetes)/pop_65_74
  bmi_65_74 = sum(bin_65_74$BMI)/pop_65_74
  
  ab40_prev_75_84 = sum(bin_75_84$Asthma)/pop_75_84
  ab34_prev_75_84 = sum(bin_75_84$Heart.Disease)/pop_75_84
  ab29_prev_75_84 = sum(bin_75_84$Hypertension)/pop_75_84
  ab52_prev_75_84 = sum(bin_75_84$Heart.Failure)/pop_75_84
  smoking_prev_75_84 = sum(bin_75_84$Smoker)/pop_75_84
  diabetes_prev_75_84 = sum(bin_75_84$Diabetes)/pop_75_84
  bmi_75_84 = sum(bin_75_84$BMI)/pop_75_84
  
  ab40_prev_85 = sum(bin_85$Asthma)/pop_85
  ab34_prev_85 = sum(bin_85$Heart.Disease)/pop_85
  ab29_prev_85 = sum(bin_85$Hypertension)/pop_85
  ab52_prev_85 = sum(bin_85$Heart.Failure)/pop_85
  smoking_prev_85 = sum(bin_85$Smoker)/pop_85
  diabetes_prev_85 = sum(bin_85$Diabetes)/pop_85
  bmi_85 = sum(bin_85$BMI)/pop_85
  
  ab40_CA = chis_ab40[5:12,4]
  ab34_CA = chis_ab34[5:12,4]
  ab29_CA = chis_ab29[5:12,4]
  ab52_CA = chis_ab52[5:12,4]
  smkcur_CA = chis_smkcur[5:12,3]
  diabetes_CA = chis_diabetes[5:12,4]
  
  bmi_CA = chis_bmi_age[1:8,3]
  
  ab40_county = as.numeric(filt_ab40_age_county[1:8,5])
  ab34_county = as.numeric(filt_ab34_age_county[1:8,5])
  ab29_county = as.numeric(filt_ab29_age_county[1:8,4])
  smkcur_county = as.numeric(filt_smkcur_age_county[1:8,5])
  diabetes_county = as.numeric(filt_diabetes_age_county[1:8,5])
  
  bmi_county = filt_bmi_age_county[1:8,4]
  
  for (i in c(1:8)) {
    if (as.numeric(filt_ab29_age_county[i,14])+as.numeric(filt_ab29_age_county[i,15]) < 10) {
      ab29_county[i] = ab29_CA[i]
    }
  }
  for (i in c(1:8)) {
    if (as.numeric(filt_ab34_age_county[i,14])+as.numeric(filt_ab34_age_county[i,15]) < 10) {
      ab34_county[i] = ab34_CA[i]
    }
  }
  for (i in c(1:8)) {
    if (as.numeric(filt_ab40_age_county[i,14])+as.numeric(filt_ab40_age_county[i,15]) < 10) {
      ab40_county[i] = ab40_CA[i]
    }
  }
  for (i in c(1:8)) {
    if (as.numeric(filt_diabetes_age_county[i,14])+as.numeric(filt_diabetes_age_county[i,15]) < 10) {
      diabetes_county[i] = diabetes_CA[i]
    }
  }
  for (i in c(1:8)) {
    if (as.numeric(filt_smkcur_age_county[i,14])+as.numeric(filt_smkcur_age_county[i,15]) < 10) {
      smkcur_county[i] = smkcur_CA[i]
    }
  }
  for (i in c(1:8)) {
    if (as.numeric(filt_bmi_age_county_totals[i,5] < 10)) {
      bmi_county[i] = bmi_CA[i]
    }
  }
  
  ab40_validation_df = data.frame("Simulated" = c(ab40_prev_18_24, ab40_prev_25_34,
                                                  ab40_prev_35_44, ab40_prev_45_54,
                                                  ab40_prev_55_64, ab40_prev_65_74,
                                                  ab40_prev_75_84, ab40_prev_85),
                                  #"Real Data" = ab40_CA)
                                  "Real Data" = ab40_county)
  ab40_validation_df = as.data.frame(t(as.matrix(ab40_validation_df)))
  colnames(ab40_validation_df) = c("18-24", "25-34", "35-44", "45-54", "
                                   55-64", "65-74", "75-84", "85+")
  barplot(as.matrix(ab40_validation_df),
          main = paste(county_name, "Asthma Prevalence by Age, Simulation vs Real Data"),
          xlab = "Age",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("bottom",
         c("Simulated","CHIS"),
         fill = c("blue","yellow")
  )
  
  ab34_validation_df = data.frame("Simulated" = c(ab34_prev_18_24, ab34_prev_25_34,
                                                  ab34_prev_35_44, ab34_prev_45_54,
                                                  ab34_prev_55_64, ab34_prev_65_74,
                                                  ab34_prev_75_84, ab34_prev_85),
                                  #"Real Data" = ab34_CA)
                                  "Real Data" = ab34_county)
  ab34_validation_df = as.data.frame(t(as.matrix(ab34_validation_df)))
  colnames(ab34_validation_df) = c("18-24", "25-34", "35-44", "45-54", "
                                   55-64", "65-74", "75-84", "85+")
  barplot(as.matrix(ab34_validation_df),
          main = paste(county_name, "Heart Disease Prevalence by Age, Simulation vs Real Data"),
          xlab = "Age",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("bottom",
         c("Simulated","CHIS Data"),
         fill = c("blue","yellow")
  )
  
  ab29_validation_df = data.frame("Simulated" = c(ab29_prev_18_24, ab29_prev_25_34,
                                                  ab29_prev_35_44, ab29_prev_45_54,
                                                  ab29_prev_55_64, ab29_prev_65_74,
                                                  ab29_prev_75_84, ab29_prev_85),
                                  #"Real Data" = ab29_CA)
                                  "Real Data" = ab29_county)
  ab29_validation_df = as.data.frame(t(as.matrix(ab29_validation_df)))
  colnames(ab29_validation_df) = c("18-24", "25-34", "35-44", "45-54", "
                                   55-64", "65-74", "75-84", "85+")
  barplot(as.matrix(ab29_validation_df),
          main = paste(county_name, "HTN Prevalence by Age, Simulation vs Real Data"),
          xlab = "Age",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("bottom",
         c("Simulated","CHIS Data"),
         fill = c("blue","yellow")
  )
  
  ab52_validation_df = data.frame("Simulated" = c(ab52_prev_18_24, ab52_prev_25_34,
                                                  ab52_prev_35_44, ab52_prev_45_54,
                                                  ab52_prev_55_64, ab52_prev_65_74,
                                                  ab52_prev_75_84, ab52_prev_85),
                                  "Real Data" = ab52_CA)
  ab52_validation_df = as.data.frame(t(as.matrix(ab52_validation_df)))
  colnames(ab52_validation_df) = c("18-24", "25-34", "35-44", "45-54", "
                                   55-64", "65-74", "75-84", "85+")
  barplot(as.matrix(ab52_validation_df),
          main = paste(county_name, "HF Prevalence by Age, Simulation vs Real Data"),
          xlab = "Age",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("bottom",
         c("Simulated","CHIS Data"),
         fill = c("blue","yellow")
  )
  
  diabetes_validation_df = data.frame("Simulated" = c(diabetes_prev_18_24, diabetes_prev_25_34,
                                                  diabetes_prev_35_44, diabetes_prev_45_54,
                                                  diabetes_prev_55_64, diabetes_prev_65_74,
                                                  diabetes_prev_75_84, diabetes_prev_85),
                                  #"Real Data" = diabetes_CA)
                                  "Real Data" = diabetes_county)
  diabetes_validation_df = as.data.frame(t(as.matrix(diabetes_validation_df)))
  colnames(diabetes_validation_df) = c("18-24", "25-34", "35-44", "45-54", "
                                   55-64", "65-74", "75-84", "85+")
  barplot(as.matrix(diabetes_validation_df),
          main = paste(county_name, "Diabetes Prevalence by Age, Simulation vs Real Data"),
          xlab = "Age",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("bottom",
         c("Simulated","CHIS Data"),
         fill = c("blue","yellow")
  )
  
  smkcur_validation_df = data.frame("Simulated" = c(smoking_prev_18_24, smoking_prev_25_34,
                                                  smoking_prev_35_44, smoking_prev_45_54,
                                                  smoking_prev_55_64, smoking_prev_65_74,
                                                  smoking_prev_75_84, smoking_prev_85),
                                  #"Real Data" = smkcur_CA)
                                  "Real Data" = smkcur_county)
  smkcur_validation_df = as.data.frame(t(as.matrix(smkcur_validation_df)))
  colnames(smkcur_validation_df) = c("18-24", "25-34", "35-44", "45-54", "
                                   55-64", "65-74", "75-84", "85+")
  barplot(as.matrix(smkcur_validation_df),
          main = paste(county_name, "Smoking Prevalence by Age, Simulation vs Real Data"),
          xlab = "Age",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("bottom",
         c("Simulated","CHIS Data"),
         fill = c("blue","yellow")
  )
  
  bmi_validation_df = data.frame("Simulated" = c(bmi_18_24, bmi_25_34,
                                                  bmi_35_44, bmi_45_54,
                                                  bmi_55_64, bmi_65_74,
                                                  bmi_75_84, bmi_85),
                                  #"Real Data" = bmi_CA)
                                 "Real Data" = bmi_county)
  bmi_validation_df = as.data.frame(t(as.matrix(bmi_validation_df)))
  colnames(ab29_validation_df) = c("18-24", "25-34", "35-44", "45-54", "
                                   55-64", "65-74", "75-84", "85+")
  barplot(as.matrix(bmi_validation_df),
          main = paste(county_name, "Avg BMI Age, Simulation vs Real Data"),
          xlab = "Age",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("bottom",
         c("Simulated","CHIS Data"),
         fill = c("blue","yellow")
  )
  return (comparison_df)
}

essential_workers_validation <- function(county_name, county_df) {
  #Expected numbers of frontline workers
  county_pop_census = filter(census, NAME == county_name)
  total_16_19 = county_pop_census[8,4] + county_pop_census[32,4]
  total_18_24 = county_pop_census[9,4]+county_pop_census[10,4]+county_pop_census[11,4]+county_pop_census[12,4]+
                county_pop_census[33,4]+county_pop_census[34,4]+county_pop_census[35,4]+county_pop_census[36,4]
  total_25_34 = county_pop_census[13,4] + county_pop_census[14,4] + county_pop_census[37,4] + county_pop_census[38,4]
  total_35_44 = county_pop_census[15,4] + county_pop_census[16,4] + county_pop_census[39,4] + county_pop_census[40,4]
  total_45_54 = county_pop_census[17,4] + county_pop_census[18,4] + county_pop_census[41,4] + county_pop_census[42,4]
  total_55_64 = county_pop_census[19,4] + county_pop_census[20,4] + county_pop_census[21,4] + county_pop_census[43,4] + county_pop_census[44,4]+county_pop_census[45,4]
  total_65_plus = county_pop_census[22,4] + county_pop_census[23,4] +county_pop_census[24,4]+county_pop_census[25,4]+
                county_pop_census[26,4]+county_pop_census[27,4]+county_pop_census[46,4]+county_pop_census[47,4]+
                county_pop_census[48,4]+county_pop_census[49,4]+county_pop_census[50,4]+county_pop_census[51,4]
  
  normalization_constant = nrow(county_df)/total_CA_pop
  expected_frontline_16_19 = CA_workers_16_19*frontline_16_19_percentage*normalization_constant
  expected_frontline_20_24 = CA_workers_20_24*frontline_20_24_percentage*normalization_constant
  expected_frontline_25_34 = CA_workers_25_34*frontline_25_34_percentage*normalization_constant
  expected_frontline_35_44 = CA_workers_35_44*frontline_35_44_percentage*normalization_constant
  expected_frontline_45_54 = CA_workers_45_54*frontline_45_54_percentage*normalization_constant
  expected_frontline_55_64 = CA_workers_55_64*frontline_55_64_percentage*normalization_constant
  expected_frontline_65_plus = CA_workers_65_plus*frontline_65_plus_percentage*normalization_constant
  
  expected_non_frontline_16_19 = CA_workers_16_19*non_frontline_16_19_percentage*normalization_constant
  expected_non_frontline_20_24 = CA_workers_20_24*non_frontline_20_24_percentage*normalization_constant
  expected_non_frontline_25_34 = CA_workers_25_34*non_frontline_25_34_percentage*normalization_constant
  expected_non_frontline_35_44 = CA_workers_35_44*non_frontline_35_44_percentage*normalization_constant
  expected_non_frontline_45_54 = CA_workers_45_54*non_frontline_45_54_percentage*normalization_constant
  expected_non_frontline_55_64 = CA_workers_55_64*non_frontline_55_64_percentage*normalization_constant
  expected_non_frontline_65_plus = CA_workers_65_plus*non_frontline_65_plus_percentage*normalization_constant
  
  
  actual_frontline_16_19 = nrow(filter(county_df, Age>14, Age<20, Special.Population == 6))
  actual_frontline_20_24 = nrow(filter(county_df, Age>19, Age<25, Special.Population == 6))
  actual_frontline_25_34 = nrow(filter(county_df, Age>24, Age<35, Special.Population == 6))
  actual_frontline_35_44 = nrow(filter(county_df, Age>34, Age<45, Special.Population == 6))
  actual_frontline_45_54 = nrow(filter(county_df, Age>44, Age<55, Special.Population == 6))
  actual_frontline_55_64 = nrow(filter(county_df, Age>54, Age<65, Special.Population == 6))
  actual_frontline_65_plus = nrow(filter(county_df, Age>64, Age<101, Special.Population == 6))
  
  actual_non_frontline_16_19 = nrow(filter(county_df, Age>14, Age<20, Special.Population == 7))
  actual_non_frontline_20_24 = nrow(filter(county_df, Age>19, Age<25, Special.Population == 7))
  actual_non_frontline_25_34 = nrow(filter(county_df, Age>24, Age<35, Special.Population == 7))
  actual_non_frontline_35_44 = nrow(filter(county_df, Age>34, Age<45, Special.Population == 7))
  actual_non_frontline_45_54 = nrow(filter(county_df, Age>44, Age<55, Special.Population == 7))
  actual_non_frontline_55_64 = nrow(filter(county_df, Age>54, Age<65, Special.Population == 7))
  actual_non_frontline_65_plus = nrow(filter(county_df, Age>64, Age<101, Special.Population == 7))
  
  
  frontline_validation_df = data.frame("Simulated" = c(actual_frontline_16_19, actual_frontline_20_24,
                                                       actual_frontline_25_34, actual_frontline_35_44,
                                                       actual_frontline_45_54, actual_frontline_55_64,
                                                       actual_frontline_65_plus),
                                  #"Real Data" = ab40_CA)
                                  "Real Data" = c(expected_frontline_16_19, expected_frontline_20_24,
                                                  expected_frontline_25_34, expected_frontline_35_44,
                                                  expected_frontline_45_54, expected_frontline_55_64,
                                                  expected_frontline_65_plus))
  frontline_validation_df = as.data.frame(t(as.matrix(frontline_validation_df)))
  colnames(frontline_validation_df) = c("16-19", "18-24", "25-34", "35-44", "
                                   45-54", "55-64", "65+")
  barplot(as.matrix(frontline_validation_df),
          main = paste(county_name, "Frontline Workers, Simulation vs Real Data"),
          xlab = "Age",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("bottom",
         c("Simulated","BLS"),
         fill = c("blue","yellow")
  )
  non_frontline_validation_df = data.frame("Simulated" = c(actual_non_frontline_16_19, actual_non_frontline_20_24,
                                                       actual_non_frontline_25_34, actual_non_frontline_35_44,
                                                       actual_non_frontline_45_54, actual_non_frontline_55_64,
                                                       actual_non_frontline_65_plus),
                                       #"Real Data" = ab40_CA)
                                       "Real Data" = c(expected_non_frontline_16_19, expected_non_frontline_20_24,
                                                       expected_non_frontline_25_34, expected_non_frontline_35_44,
                                                       expected_non_frontline_45_54, expected_non_frontline_55_64,
                                                       expected_non_frontline_65_plus))
  non_frontline_validation_df = as.data.frame(t(as.matrix(non_frontline_validation_df)))
  colnames(non_frontline_validation_df) = c("16-19", "18-24", "25-34", "35-44", "
                                   45-54", "55-64", "65+")
  barplot(as.matrix(non_frontline_validation_df),
          main = paste(county_name, "non_frontline Workers, Simulation vs Real Data"),
          xlab = "Age",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("bottom",
         c("Simulated","BLS"),
         fill = c("blue","yellow")
  )
}

ALF_validation <- function(county_name, county_df){
  normalization_constant = nrow(county_df)/total_CA_pop
  expected_ALF_55_64 = 0.03*150000*normalization_constant
  expected_ALF_65_74 = 0.109*150000*normalization_constant
  expected_ALF_75_84 = 0.319*150000*normalization_constant
  expected_ALF_85_plus = 0.548*150000*normalization_constant
  
  actual_ALF_55_64 = nrow(filter(county_df, Special.Population == 8, Age >=55, Age <= 64))
  actual_ALF_65_74 = nrow(filter(county_df, Special.Population == 8, Age >=65, Age <= 74))
  actual_ALF_75_84 = nrow(filter(county_df, Special.Population == 8, Age >=75, Age <= 84))
  actual_ALF_85_plus = nrow(filter(county_df, Special.Population == 8, Age >=85))
  
  frontline_validation_df = data.frame("Simulated" = c(actual_ALF_55_64, actual_ALF_65_74,
                                                       actual_ALF_75_84, actual_ALF_85_plus),
                                       #"Real Data" = ab40_CA)
                                       "Real Data" = c(expected_ALF_55_64, expected_ALF_65_74,
                                                       expected_ALF_75_84, expected_ALF_85_plus))
  frontline_validation_df = as.data.frame(t(as.matrix(frontline_validation_df)))
  colnames(frontline_validation_df) = c("55-64", "65-74", "75-84", "85+")
  barplot(as.matrix(frontline_validation_df),
          main = paste(county_name, "ALF Residents, Simulation vs Real Data"),
          xlab = "Age",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("top",
         c("Simulated","Data"),
         fill = c("blue","yellow")
  )
}

special_populations_validation <- function(county_name, county_df) {
  # 0 for no special population
  # 1 for healthcare workers
  # 2 for incarcerated
  # 3 for SNF
  # 4 for educational occupations
  # 5 for homeless
  
  #county_name  ="Yuba County, California"
  #county_df = Yuba_df
  
  num_hcw = nrow(filter(county_df, Special.Population == 1))
  num_prisoners = nrow(filter(county_df, Special.Population == 2))
  num_SNF = nrow(filter(county_df, Special.Population == 3))
  num_educ = nrow(filter(county_df, Special.Population == 4))
  num_homeless = nrow(filter(county_df, Special.Population == 5))
  
  occ_filtered = filter(census_occ, NAME == county_name)
  pop_filtered = filter(census, NAME == county_name)
  num_male_educ = occ_filtered[1,4] * frac_without_archivists
  num_male_hcw = occ_filtered[2,4] + occ_filtered[5,4]
  num_female_educ = occ_filtered[6,4] * frac_without_archivists
  num_female_hcw = occ_filtered[7,4] + occ_filtered[10,4]
  
  county_homeless_df = filter(county_df, Special.Population == 5)
  num_homeless_white = nrow(filter(county_homeless_df, Race == "White"))
  num_homeless_african_american = nrow(filter(county_homeless_df, Race == "African American"))
  num_homeless_AIAN = nrow(filter(county_homeless_df, Race == "AIAN"))
  num_homeless_asian = nrow(filter(county_homeless_df, Race == "Asian Alone"))
  num_homeless_NHPI = nrow(filter(county_homeless_df, Race == "Native Hawaiian And Other Pacific Islander Alone"))
  num_homeless_other = nrow(filter(county_homeless_df, Race == "Some other race alone"))
  num_homeless_multi = nrow(filter(county_homeless_df, Race == "Two or more races"))
  
  county_homeless_by_race = c(num_homeless_white, num_homeless_african_american, 
                              num_homeless_AIAN, num_homeless_asian, num_homeless_NHPI,
                              num_homeless_other+num_homeless_multi)
  
  homeless_county_white = filter(homeless_data, CA_county == county_name)[1, 28]
  homeless_county_african_american = filter(homeless_data, CA_county == county_name)[1, 29]
  homeless_county_asian = filter(homeless_data, CA_county == county_name)[1, 30]
  homeless_county_AIAN = filter(homeless_data, CA_county == county_name)[1, 31]
  homeless_county_NHPI = filter(homeless_data, CA_county == county_name)[1, 32]
  homeless_county_multiple = filter(homeless_data, CA_county == county_name)[1, 33]

  county_homeless_CoC = c(homeless_county_white, homeless_county_african_american, homeless_county_AIAN,
                          homeless_county_asian, homeless_county_NHPI, homeless_county_multiple)
  
  homeless_validation_df = data.frame("Simulated" = county_homeless_by_race,
                                  "Real Data" = county_homeless_CoC)
  homeless_validation_df = as.data.frame(t(as.matrix(homeless_validation_df)))
  
  colnames(homeless_validation_df) = c("Wh", "Af_Am", "AIAN", "Asian", "NHPI", "Multi")
  barplot(as.matrix(homeless_validation_df),
          main = paste(county_name, "Homelessness by Race, Simulation vs Real Data"),
          xlab = "Race",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("top",
         c("Simulated","CoC Data"),
         fill = c("blue","yellow")
  )
  county_SNF_df = filter(county_df, Special.Population == 3)
  num_SNF_white = nrow(filter(county_SNF_df, Race == "White"))
  num_SNF_african_american = nrow(filter(county_SNF_df, Race == "African American"))
  num_SNF_AIAN = nrow(filter(county_SNF_df, Race == "AIAN"))
  num_SNF_asian = nrow(filter(county_SNF_df, Race == "Asian Alone"))
  num_SNF_NHPI = nrow(filter(county_SNF_df, Race == "Native Hawaiian And Other Pacific Islander Alone"))
  num_SNF_other = nrow(filter(county_SNF_df, Race == "Some other race alone"))
  num_SNF_multi = nrow(filter(county_SNF_df, Race == "Two or more races"))
  num_SNF_hispanic = nrow(filter(county_SNF_df, Ethnicity == "Hispanic or Latino"))
  num_SNF_not_hispanic = nrow(filter(county_SNF_df, Ethnicity == "Not Hispanic or Latino"))
  
  county_SNF_race = c(num_SNF_white, num_SNF_african_american, num_SNF_AIAN,
                      num_SNF_asian, num_SNF_NHPI, num_SNF_other, num_SNF_multi,
                      num_SNF_hispanic, num_SNF_not_hispanic)
  
  
  county_CF_df = filter(county_df, Special.Population == 2)
  num_CF_white = nrow(filter(county_CF_df, Race == "White"))
  num_CF_african_american = nrow(filter(county_CF_df, Race == "African American"))
  num_CF_AIAN = nrow(filter(county_CF_df, Race == "AIAN"))
  num_CF_asian = nrow(filter(county_CF_df, Race == "Asian Alone"))
  num_CF_NHPI = nrow(filter(county_CF_df, Race == "Native Hawaiian And Other Pacific Islander Alone"))
  num_CF_other = nrow(filter(county_CF_df, Race == "Some other race alone"))
  num_CF_multi = nrow(filter(county_CF_df, Race == "Two or more races"))
  num_CF_hispanic = nrow(filter(county_CF_df, Ethnicity == "Hispanic or Latino"))
  num_CF_not_hispanic = nrow(filter(county_CF_df, Ethnicity == "Not Hispanic or Latino"))
  
  county_CF_race = c(num_CF_white, num_CF_african_american, num_CF_AIAN,
                      num_CF_asian, num_CF_NHPI, num_CF_other, num_CF_multi,
                      num_CF_hispanic, num_CF_not_hispanic)
  
  pop_weight = nrow(county_df)/total_CA_pop
  
  SNF_18_24 = census_SNF[2,4]*pop_weight
  SNF_25_34 = census_SNF[3,4]*pop_weight
  SNF_35_44 = census_SNF[4,4]*pop_weight
  SNF_45_54 = census_SNF[5,4]*pop_weight
  SNF_55_64 = census_SNF[6,4]*pop_weight
  SNF_65_74 = census_SNF[7,4]*pop_weight
  SNF_75_84 = census_SNF[8,4]*pop_weight
  SNF_85_plus = census_SNF[9,4]*pop_weight
  
  SNF_ACS_age = c(SNF_18_24, SNF_25_34, SNF_35_44, SNF_45_54, SNF_55_64, SNF_65_74,
                  SNF_75_84, SNF_85_plus)
  
  num_SNF_18_24 = nrow(filter(county_df, Special.Population == 3, Age>=18 & Age <=24))
  num_SNF_25_34 = nrow(filter(county_df, Special.Population == 3, Age>=25 & Age <=34))
  num_SNF_35_44 = nrow(filter(county_df, Special.Population == 3, Age>=35 & Age <=44))
  num_SNF_45_54 = nrow(filter(county_df, Special.Population == 3, Age>=45 & Age <=54))
  num_SNF_55_64 = nrow(filter(county_df, Special.Population == 3, Age>=55 & Age <=64))
  num_SNF_65_74 = nrow(filter(county_df, Special.Population == 3, Age>=65 & Age <=74))
  num_SNF_75_84 = nrow(filter(county_df, Special.Population == 3, Age>=75 & Age <=84))
  num_SNF_85_plus = nrow(filter(county_df, Special.Population == 3, Age>=85))
  
  county_SNF_age = c(num_SNF_18_24, num_SNF_25_34, num_SNF_35_44, num_SNF_45_54, 
                     num_SNF_55_64, num_SNF_65_74, num_SNF_75_84, num_SNF_85_plus)
  
  CF_15_17 = census_CF[2,4]*pop_weight
  CF_18_24 = census_CF[3,4]*pop_weight
  CF_25_34 = census_CF[4,4]*pop_weight
  CF_35_44 = census_CF[5,4]*pop_weight
  CF_45_54 = census_CF[6,4]*pop_weight
  CF_55_64 = census_CF[7,4]*pop_weight
  CF_65_74 = census_CF[8,4]*pop_weight
  CF_75_84 = census_CF[9,4]*pop_weight
  CF_85_plus = census_CF[10,4]*pop_weight
  
  CF_ACS_age = c(CF_15_17, CF_18_24, CF_25_34, CF_35_44, CF_45_54, CF_55_64,
                  CF_65_74, CF_75_84, CF_85_plus)
  
  num_CF_15_17 = nrow(filter(county_df, Special.Population == 2, Age>=15 & Age <=17))
  num_CF_18_24 = nrow(filter(county_df, Special.Population == 2, Age>=18 & Age <=24))
  num_CF_25_34 = nrow(filter(county_df, Special.Population == 2, Age>=25 & Age <=34))
  num_CF_35_44 = nrow(filter(county_df, Special.Population == 2, Age>=35 & Age <=44))
  num_CF_45_54 = nrow(filter(county_df, Special.Population == 2, Age>=45 & Age <=54))
  num_CF_55_64 = nrow(filter(county_df, Special.Population == 2, Age>=55 & Age <=64))
  num_CF_65_74 = nrow(filter(county_df, Special.Population == 2, Age>=65 & Age <=74))
  num_CF_75_84 = nrow(filter(county_df, Special.Population == 2, Age>=75 & Age <=84))
  num_CF_85_plus = nrow(filter(county_df, Special.Population == 2, Age>=85))
  
  county_CF_age = c(num_CF_15_17, num_CF_18_24, num_CF_25_34,
                     num_CF_35_44, num_CF_45_54, num_CF_55_64, num_CF_65_74,
                     num_CF_75_84, num_CF_85_plus)
  
  CF_white = census_SNF[35,4]*pop_weight
  CF_african_american = census_SNF[43,4]*pop_weight
  CF_AIAN = census_SNF[51,4]*pop_weight
  CF_asian = census_SNF[59,4]*pop_weight
  CF_NHPI = census_SNF[67,4]*pop_weight
  CF_other = census_SNF[75,4]*pop_weight
  CF_multi = census_SNF[83,4]*pop_weight
  
  
  SNF_white = census_SNF[36,4]*pop_weight
  SNF_african_american = census_SNF[44,4]*pop_weight
  SNF_AIAN = census_SNF[52,4]*pop_weight
  SNF_asian = census_SNF[60,4]*pop_weight
  SNF_NHPI = census_SNF[68,4]*pop_weight
  SNF_other = census_SNF[76,4]*pop_weight
  SNF_multi = census_SNF[84,4]*pop_weight
  
  SNF_hispanic = census_SNF[100,4]*pop_weight
  SNF_not_hispanic = (census_SNF[1,4] - census_SNF[100,4])*pop_weight
  
  SNF_ACS_race = c(SNF_white, SNF_african_american, SNF_AIAN, SNF_asian, 
                   SNF_NHPI, SNF_other, SNF_multi, SNF_hispanic, SNF_not_hispanic)
  
  CF_hispanic = census_CF[103,4]*pop_weight
  CF_not_hispanic = census_CF[1,4]*pop_weight - CF_hispanic
  
  CF_ACS_race = c(CF_white, CF_african_american, CF_AIAN, CF_asian, CF_NHPI,
                  CF_other, CF_multi, CF_hispanic, CF_not_hispanic)
  #TODO
  #need to make graphs of healthcare workers by age and sex
  #do the same for education workers
  #need to make graphs of prisoners by age, sex, race
  #need to measure number of homeless ppl
  #need to determine if racial disparities are reflected in simulated homeless pop
  #need to make graphs of SNF by age, sex, race
  
  SNF_race_validation_df = data.frame("Simulated" = county_SNF_race,
                                      "Real Data" = SNF_ACS_race)
  SNF_race_validation_df = as.data.frame(t(as.matrix(SNF_race_validation_df)))
  
  colnames(SNF_race_validation_df) = c("Wh", "Af_Am", "AIAN", "Asn", "NHPI", "Oth","Multi",
                                       "H", "NH")
  barplot(as.matrix(SNF_race_validation_df),
          main = paste(county_name, "SNF by Race, Simulation vs Real Data"),
          xlab = "Race",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("top",
         c("Simulated","ACS Data"),
         fill = c("blue","yellow")
  )
  
  SNF_age_validation_df = data.frame("Simulated" = county_SNF_age,
                                      "Real Data" = SNF_ACS_age)
  SNF_age_validation_df = as.data.frame(t(as.matrix(SNF_age_validation_df)))
  
  colnames(SNF_age_validation_df) = c("<25", "<35", "<45", "<55", "<65", "<75", "<85", "<100")
  barplot(as.matrix(SNF_age_validation_df),
          main = paste(county_name, "SNF by Age, Simulation vs Real Data"),
          xlab = "Age",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("top",
         c("Simulated","ACS Data"),
         fill = c("blue","yellow")
  )
  CF_race_validation_df = data.frame("Simulated" = county_CF_race,
                                      "Real Data" = CF_ACS_race)
  CF_race_validation_df = as.data.frame(t(as.matrix(CF_race_validation_df)))
  
  colnames(CF_race_validation_df) = c("Wh", "Af_Am", "AIAN", "Asn", "NHPI", "Oth","Multi",
                                       "H", "NH")
  barplot(as.matrix(CF_race_validation_df),
          main = paste(county_name, "CF by Race, Simulation vs Real Data"),
          xlab = "Race",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("top",
         c("Simulated","ACS Data"),
         fill = c("blue","yellow")
  )
  CF_age_validation_df = data.frame("Simulated" = county_CF_age,
                                     "Real Data" = CF_ACS_age)
  CF_age_validation_df = as.data.frame(t(as.matrix(CF_age_validation_df)))
  
  colnames(CF_age_validation_df) = c("<18", "<25", "<35", "<45", "<55", "<65", "<75", "<85", "<100")
  barplot(as.matrix(CF_age_validation_df),
          main = paste(county_name, "CF by Age, Simulation vs Real Data"),
          xlab = "Age",
          col = c("blue","yellow"), beside = TRUE
  )
  legend("top",
         c("Simulated","ACS Data"),
         fill = c("blue","yellow")
  )
  print(num_hcw -  num_male_hcw -  num_female_hcw)
  print(num_educ - num_male_educ - num_female_educ)
  print(num_SNF)
  print(num_homeless)
  
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Load census data 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

my_key <-"665e6ac52698b6e2518d187bc9865bdc1c4192eb"
census_api_key(my_key, overwrite=T, install=T)
readRenviron("~/.Renviron")

cpsaat11b <- read_excel("../Data/cpsaat11b.xlsx")
hcw_age_data = data.frame(cpsaat11b)

load_variables <- load_variables(year=2018,dataset="acs5")
test <- load_variables
var <- test$name[c(1:330, 590:610)] # Extract section of possibly relevant variables 
var_label1 <- rep(test$label[c(1:330, 590:610)], 58) # 58 counties in CA 
var_label2 <- rep(test$concept[c(1:330, 590:610)],58) # 58 counties in CA 

census <- get_acs(geography = "county",
                  variables =  var,
                  state = "CA",
                  year = 2018, key = my_key) # obtain data 
census <- data.frame(census, label=var_label1, concept=var_label2) # label variables 

var_label1 <- rep(test$label[c(1:330, 590:610)], 1) # 58 counties in CA 
var_label2 <- rep(test$concept[c(1:330, 590:610)],1) # 58 counties in CA 
census_CA = get_acs(geography = "state",
                    variables =  var,
                    state = "CA",
                    year = 2018, key = my_key) # obtain data 
census_CA <- data.frame(census_CA, label=var_label1, concept=var_label2) # label variables 

load_variables2017 <- load_variables(year=2017,dataset="acs5")
load_variables2016 <- load_variables(year=2016,dataset="acs5")

SNF_var <- test$name[c(22367:22393, 22464:22538, 23254:23280, 23392:23482)]
SNF_label1 <- rep(test$label[c(22367:22393, 22464:22538, 23254:23280, 23392:23482)], 1) # 58 counties in CA 
SNF_label2 <- rep(test$concept[c(22367:22393, 22464:22538, 23254:23280, 23392:23482)],1) # 58 counties in CA 

census_SNF <- get_acs(geography = "state",
                      variables =  SNF_var,
                      state = "CA",
                      year = 2018, key = my_key) # obtain data 
census_SNF <- data.frame(census_SNF, label=SNF_label1, concept=SNF_label2) # label variables 

CF_var <- test$name[c(22337:22366, 22463:22536, 23224:23253, 23391:23481)]
CF_label1 <- rep(test$label[c(22337:22366, 22463:22536, 23224:23253, 23391:23481)], 1) # 58 counties in CA 
CF_label2 <- rep(test$concept[c(22337:22366, 22463:22536, 23224:23253, 23391:23481)],1) # 58 counties in CA 


census_CF <- get_acs(geography = "state",
                     variables =  CF_var,
                     state = "CA",
                     year = 2018, key = my_key) # obtain data 
census_CF <- data.frame(census_CF, label=CF_label1, concept=CF_label2) # label variables 



load_variables2018 <- load_variables(year=2018,dataset="acs1")

var_occ <- test$name[c(26007,26009,26010,26011,26013,26043,26045,26046,26047, 26049)] # Extract section of possibly relevant variables 
var_occ_label1 <- rep(test$label[c(26007,26009,26010,26011,26013,26043,26045,26046,26047, 26049)], 58) # 58 counties in CA 
var_occ_label2 <- rep(test$concept[c(26007,26009,26010,26011,26013,26043,26045,26046,26047, 26049)],58) # 58 counties in CA 

census_occ <- get_acs(geography = "county",
                      variables =  var_occ,
                      state = "CA",
                      year = 2018, survey = "acs5") # obtain data 
census_occ <- data.frame(census_occ, label=var_occ_label1, concept=var_occ_label2) # label variables 

CA_occ_totals <- get_acs(geography = "state",
                         variables =  var_occ,
                         state = "CA",
                         year = 2018, survey = "acs5")
CA_occ_totals <- data.frame(CA_occ_totals, label=var_occ_label1, concept=var_occ_label2)

CA_wide_hcw_total = CA_occ_totals[2,4] + CA_occ_totals[5,4] + CA_occ_totals[17,4] + CA_occ_totals[20,4]
CA_wide_educ_total = CA_occ_totals[1,4] + CA_occ_totals[6,4]
#B24010_043 to B24010_051 for teachers, librarians Males, 16119 to 16127      
#B24010_056 to B24010_063 for healthcare Males, 16132 to 16138
#B24010_065 to B24010_068 healthcare support Males 16141 to 16144
#B24010_194 to B24010_202 teachers, librarians, female 16270 to 16278
#B24010_207 to B24010_214 for healthcare, female 16283 to 16290
#B24010_216 to B24010_219 for healthcare support, female 16292 to 16295

#B24010A_014 white alone, total education and librarian occupation, male 16393 
#B24010A_017, B24010A_018 healthcare diagnosing and tech, white alone, male 16396,16397 16398
#B24010A_020 healtchare support, white alone, male 16399
#B24010A_050 education, white alone, female 16329
#B24010A_052 to B24010A_054 white alone, healthcare, female 16431 to 16433
#B24010A_056 healthcare support, white alone, female 16435

#apply offset of 73 for all subsequent groups

#same as above, but B for african americans
#C for AIAN
#D for asian alone
#E for native hawaiian
#F for other
#G for multi
var_occ_sex_alone = load_variables2018$name[c(16119:16127,16132:16138,16141:16144,16270:16278,16283:16290,16292:16295)]
var_occ_sex_and_race_1 = c(16393,16396,16397,16398,16399,16329,16431:16433,16435)

var_occ_sex_and_race = c(var_occ_sex_and_race_1, var_occ_sex_and_race_1+73) #white, african american
var_occ_sex_and_race = c(var_occ_sex_and_race, var_occ_sex_and_race_1+2*73) #AIAN
var_occ_sex_and_race = c(var_occ_sex_and_race, var_occ_sex_and_race_1+3*73) #asian
var_occ_sex_and_race = c(var_occ_sex_and_race, var_occ_sex_and_race_1+4*73) #native hawaiian
var_occ_sex_and_race = c(var_occ_sex_and_race, var_occ_sex_and_race_1+5*73) #other
var_occ_sex_and_race = c(var_occ_sex_and_race, var_occ_sex_and_race_1+6*73) #multi
var_occ_sex_and_race = c(var_occ_sex_and_race, var_occ_sex_and_race_1+8*73) #hispanic or latino

var_occ_sex_and_race_names = load_variables2018$name[var_occ_sex_and_race]

var_occ_sex_alone_l1 <- rep(load_variables2018$label[c(16119:16127,16132:16138,16141:16144,16270:16278,16283:16290,16292:16295)], 1) 
var_occ_sex_alone_l2 <- rep(load_variables2018$concept[c(16119:16127,16132:16138,16141:16144,16270:16278,16283:16290,16292:16295)],1) 

var_occ_sex_and_race_l1 <- rep(load_variables2018$label[var_occ_sex_and_race], 40) # 58 counties in CA 
var_occ_sex_and_race_l2 <- rep(load_variables2018$concept[var_occ_sex_and_race],40) # 58 counties in CA 


census_occ_sex_alone <- get_acs(geography = "state",
                                variables =  var_occ_sex_alone,
                                state = "CA",
                                year = 2018, survey = "acs1") # obtain data 
census_occ_sex_alone <- data.frame(census_occ_sex_alone, label=var_occ_sex_alone_l1, concept=var_occ_sex_alone_l2) # label variables 

census_occ_sex_and_race <- get_acs(geography = "state",
                                   variables =  var_occ_sex_and_race_names,
                                   state = "CA",
                                   year = 2018, survey = "acs1") # obtain data 
census_occ_sex_and_race <- data.frame(census_occ_sex_and_race, label=var_occ_sex_and_race_l1, concept=var_occ_sex_and_race_l2) # label variables 


numCores = detectCores()
registerDoParallel(numCores) 
################################################################################
#LA Almanac for homelessness data
CA_prop_male = census_CA[4,4]/census_CA[3,4]
CA_prop_female = census_CA[28,4]/census_CA[3,4]
num_homeless_CA = 151278 #from HUD estimate
num_homeless_LA = 66436
homeless_data = data.frame(read_excel("../Data/CoC by county.xlsx"))


prop_male_homeless_CA = 98404/(98404+50467)
prop_female_homeless_CA = 50467/(98404+50467)

prop_homeless_white = 0.254
prop_homeless_african_american = 0.338
prop_homeless_AIAN = 0.011
prop_homeless_asian = 0.012
prop_homeless_native_hawaiian_pacific_islander = 0.003
prop_homeless_other_race= 0.361 #same as prop hispanic/latino
prop_homeless_multi = 0.021

prop_hispanic_homeless = 0.361
prop_not_hispanic_homeless = 1- prop_hispanic_homeless

prop_homeless_under_18 = 7491/num_homeless_LA #this seems really high
prop_homeless_18_24 = 4181/num_homeless_LA
prop_homeless_25_54 = 37138/num_homeless_LA
prop_homeless_55_61 = 8606/num_homeless_LA
prop_homeless_62_plus = 6290/num_homeless_LA
prop_sums_age = prop_homeless_under_18 + prop_homeless_18_24+prop_homeless_25_54+prop_homeless_55_61+prop_homeless_62_plus

CA_homeless_props = c(0, 0, 0, prop_homeless_under_18, prop_homeless_18_24,
                      prop_homeless_25_54/3, prop_homeless_25_54/3, prop_homeless_25_54/3,
                      prop_homeless_55_61, prop_homeless_62_plus*2/3, prop_homeless_62_plus*1/3, 0)

CA_homeless_props_race = c(prop_homeless_white, prop_homeless_african_american,
                           prop_homeless_AIAN, prop_homeless_asian, prop_homeless_native_hawaiian_pacific_islander,
                           prop_homeless_other_race, prop_homeless_multi)

#Frontline and non-frontline essential workers
total_CA_workers = 18309012
management_and_science_CA = 7076526
service_CA = 3393502
sales_CA = 3968498
natural_resources_CA = 1653964
production_CA = 2216522


BLS_worker_total = strtoi(hcw_age_data[6,2])
BLS_16_19 = strtoi(hcw_age_data[6,3])
BLS_20_24 = strtoi(hcw_age_data[6,4])
BLS_25_34 = strtoi(hcw_age_data[6,5])
BLS_35_44 = strtoi(hcw_age_data[6,6])
BLS_45_54 = strtoi(hcw_age_data[6,7])
BLS_55_64 = strtoi(hcw_age_data[6,8])
BLS_65_plus = strtoi(hcw_age_data[6,9])

#Entire production category is essential and frontline
prod_and_transport = strtoi(hcw_age_data[456,2])

#Under service occupations (entire category minus personal care, but plus funerary services)
service = strtoi(hcw_age_data[216,2]) - strtoi(hcw_age_data[269,2]) + strtoi(hcw_age_data[278,2]) + strtoi(hcw_age_data[279,2])
frontline_BLS_service_percentage = service/strtoi(hcw_age_data[216,2])
#Management,business, science, art
healthcare_technical = strtoi(hcw_age_data[181,2])
community_and_social = strtoi(hcw_age_data[134,2])
life_physical = strtoi(hcw_age_data[110,2])
frontline_BLS_management_percentage = (healthcare_technical+community_and_social+life_physical)/strtoi(hcw_age_data[8,2])

#Sales and office related
sales = strtoi(hcw_age_data[292,2])
frontline_BLS_sales_percentage = sales/strtoi(hcw_age_data[291,2])

#Natural resources (this entire category is essential and frontline)
nat_resources = strtoi(hcw_age_data[365,2])

#Frontline total
frontline_BLS = prod_and_transport + service + healthcare_technical + community_and_social+ life_physical + sales + nat_resources
frontline_BLS_percentage = frontline_BLS/BLS_worker_total
CA_total_frontline = (natural_resources_CA+frontline_BLS_sales_percentage*sales_CA+
                        frontline_BLS_management_percentage*management_and_science_CA+production_CA+
                        frontline_BLS_service_percentage*service_CA)
CA_frontline_percentage = CA_total_frontline/total_CA_workers
##Non-frontline essential workers
#management, business, science
educators = strtoi(hcw_age_data[149,2])
management = strtoi(hcw_age_data[10,2])
computers = strtoi(hcw_age_data[71,2])
architecture_engineering = strtoi(hcw_age_data[88,2])
finance = strtoi(hcw_age_data[41,2])
entertainment = strtoi(hcw_age_data[161,2])
law = strtoi(hcw_age_data[143,2])

BLS_essential_management_percentage = (educators+management+computers+architecture_engineering+finance+
                                         entertainment+law)/strtoi(hcw_age_data[8,2])

#Sales and office occupations
office_admin = strtoi(hcw_age_data[311,2])
BLS_essential_sales_percentage = office_admin/strtoi(hcw_age_data[291,2])

#Total Essential, non-frontline workers
non_frontline_BLS = educators + management + computers + architecture_engineering + finance + entertainment+law
non_frontline_BLS_percentage = non_frontline_BLS/BLS_worker_total

BLS_essential_percentage = non_frontline_BLS_percentage + frontline_BLS_percentage
CA_essential_total = CA_total_frontline + BLS_essential_sales_percentage*sales_CA+BLS_essential_management_percentage*management_and_science_CA
CA_essential_percentage = CA_essential_total/total_CA_workers



#Reclassification of workers under CDC definitions
#CDC frontline
#Service Occupations
first_responders = strtoi(hcw_age_data[229,2])
frontline_service_16_19 = strtoi(hcw_age_data[229,3])
frontline_service_20_24 = strtoi(hcw_age_data[229,4])
frontline_service_25_34 = strtoi(hcw_age_data[229,5])
frontline_service_35_44 = strtoi(hcw_age_data[229,6])
frontline_service_45_54 = strtoi(hcw_age_data[229,7])
frontline_service_55_64 = strtoi(hcw_age_data[229,8])
frontline_service_65_plus = strtoi(hcw_age_data[229,9])

frontline_service_percentage = first_responders/strtoi(hcw_age_data[216,2])

#Management and professional
education = strtoi(hcw_age_data[149,2])
frontline_management_16_19 = strtoi(hcw_age_data[149,3])
frontline_management_20_24 = strtoi(hcw_age_data[149,4])
frontline_management_25_34 = strtoi(hcw_age_data[149,5])
frontline_management_35_44 = strtoi(hcw_age_data[149,6])
frontline_management_45_54 = strtoi(hcw_age_data[149,7])
frontline_management_55_64 = strtoi(hcw_age_data[149,8])
frontline_management_65_plus = strtoi(hcw_age_data[149,9])
frontline_management_percentage = education/strtoi(hcw_age_data[8,2])

#Zero out education workers because they're already their own separate category
frontline_management_16_19 = 0
frontline_management_20_24 = 0
frontline_management_25_34 = 0
frontline_management_35_44 = 0
frontline_management_45_54 = 0
frontline_management_55_64 = 0
frontline_management_65_plus = 0
frontline_management_percentage = 0

#Natural resources
food_ag = strtoi(hcw_age_data[366,2])
frontline_nat_res_16_19 = strtoi(hcw_age_data[366,3])
frontline_nat_res_20_24 = strtoi(hcw_age_data[366,4])
frontline_nat_res_25_34 = strtoi(hcw_age_data[366,5])
frontline_nat_res_35_44 = strtoi(hcw_age_data[366,6])
frontline_nat_res_45_54 = strtoi(hcw_age_data[366,7])
frontline_nat_res_55_64 = strtoi(hcw_age_data[366,8])
frontline_nat_res_65_plus = strtoi(hcw_age_data[366,9])
frontline_nat_resources_percentage = food_ag/strtoi(hcw_age_data[365,2])

#production and transportation
manufacturing = strtoi(hcw_age_data[457,2])
public_transit = strtoi(hcw_age_data[545,2]) + strtoi(hcw_age_data[552,2]) + strtoi(hcw_age_data[547,2])

frontline_production_16_19 = strtoi(hcw_age_data[457,3])+strtoi(hcw_age_data[545,3]) + strtoi(hcw_age_data[552,3]) + strtoi(hcw_age_data[547,3])
frontline_production_20_24 = strtoi(hcw_age_data[457,4])+strtoi(hcw_age_data[545,4]) + strtoi(hcw_age_data[552,4]) + strtoi(hcw_age_data[547,4])
frontline_production_25_34 = strtoi(hcw_age_data[457,5])+strtoi(hcw_age_data[545,5]) + strtoi(hcw_age_data[552,5]) + strtoi(hcw_age_data[547,5])
frontline_production_35_44 = strtoi(hcw_age_data[457,6])+strtoi(hcw_age_data[545,6]) + strtoi(hcw_age_data[552,6]) + strtoi(hcw_age_data[547,6])
frontline_production_45_54 = strtoi(hcw_age_data[457,7])+strtoi(hcw_age_data[545,7]) + strtoi(hcw_age_data[552,7]) + strtoi(hcw_age_data[547,7])
frontline_production_55_64 = strtoi(hcw_age_data[457,8])+strtoi(hcw_age_data[545,8]) + strtoi(hcw_age_data[552,8]) + strtoi(hcw_age_data[547,8])
frontline_production_65_plus = strtoi(hcw_age_data[457,9])+strtoi(hcw_age_data[545,9]) + strtoi(hcw_age_data[552,9]) + strtoi(hcw_age_data[547,9])

frontline_production_percentage = (manufacturing+public_transit)/strtoi(hcw_age_data[456,2])

#Sales and office
USPS = strtoi(hcw_age_data[345,2])+strtoi(hcw_age_data[346,2])+strtoi(hcw_age_data[347,2])
grocery_clerks = strtoi(hcw_age_data[293,2])+strtoi(hcw_age_data[294,2])+strtoi(hcw_age_data[295,2])+strtoi(hcw_age_data[296,2])

frontline_sales_16_19 = strtoi(hcw_age_data[345,3])+strtoi(hcw_age_data[346,3])+strtoi(hcw_age_data[347,3])+
                        strtoi(hcw_age_data[293,3])+strtoi(hcw_age_data[294,3])+strtoi(hcw_age_data[295,3])+
                        strtoi(hcw_age_data[296,3])
frontline_sales_20_24 = strtoi(hcw_age_data[345,4])+strtoi(hcw_age_data[346,4])+strtoi(hcw_age_data[347,4])+
                        strtoi(hcw_age_data[293,4])+strtoi(hcw_age_data[294,4])+strtoi(hcw_age_data[295,4])+
                        strtoi(hcw_age_data[296,4])
frontline_sales_25_34 = strtoi(hcw_age_data[345,5])+strtoi(hcw_age_data[346,5])+strtoi(hcw_age_data[347,5])+
                        strtoi(hcw_age_data[293,5])+strtoi(hcw_age_data[294,5])+strtoi(hcw_age_data[295,5])+
                        strtoi(hcw_age_data[296,5])
frontline_sales_35_44 = strtoi(hcw_age_data[345,6])+strtoi(hcw_age_data[346,6])+strtoi(hcw_age_data[347,6])+
                        strtoi(hcw_age_data[293,6])+strtoi(hcw_age_data[294,6])+strtoi(hcw_age_data[295,6])+
                        strtoi(hcw_age_data[296,6])
frontline_sales_45_54 = strtoi(hcw_age_data[345,7])+strtoi(hcw_age_data[346,7])+strtoi(hcw_age_data[347,7])+
                        strtoi(hcw_age_data[293,7])+strtoi(hcw_age_data[294,7])+strtoi(hcw_age_data[295,7])+
                        strtoi(hcw_age_data[296,7])
frontline_sales_55_64 = strtoi(hcw_age_data[345,8])+strtoi(hcw_age_data[346,8])+strtoi(hcw_age_data[347,8])+
                        strtoi(hcw_age_data[293,8])+strtoi(hcw_age_data[294,8])+strtoi(hcw_age_data[295,8])+
                        strtoi(hcw_age_data[296,8])
frontline_sales_65_plus = strtoi(hcw_age_data[345,9])+strtoi(hcw_age_data[346,9])+strtoi(hcw_age_data[347,9])+
                          strtoi(hcw_age_data[293,9])+strtoi(hcw_age_data[294,9])+strtoi(hcw_age_data[295,9])+
                          strtoi(hcw_age_data[296,9])
frontline_sales_percentage = (USPS+grocery_clerks)/strtoi(hcw_age_data[291,2])

CDC_frontline_total = first_responders+education+food_ag+manufacturing+USPS+grocery_clerks+public_transit
CDC_frontline_percentage = CDC_frontline_total/strtoi(hcw_age_data[6,2])

CA_frontline_percentage = (service_CA*frontline_service_percentage + production_CA*frontline_production_percentage+
                          management_and_science_CA*frontline_management_percentage +frontline_nat_resources_percentage*natural_resources_CA+
                          sales_CA*frontline_sales_percentage)/total_CA_workers

frontline_16_19_percentage = (frontline_management_16_19+frontline_sales_16_19+frontline_nat_res_16_19+
                                frontline_production_16_19+frontline_service_16_19)/strtoi(hcw_age_data[6,3])
frontline_20_24_percentage = (frontline_management_20_24+frontline_sales_20_24+frontline_nat_res_20_24+
                                frontline_production_20_24+frontline_service_20_24)/strtoi(hcw_age_data[6,4])
frontline_25_34_percentage = (frontline_management_25_34+frontline_sales_25_34+frontline_nat_res_25_34+
                                frontline_production_25_34+frontline_service_25_34)/strtoi(hcw_age_data[6,5])
frontline_35_44_percentage = (frontline_management_35_44+frontline_sales_35_44+frontline_nat_res_35_44+
                                frontline_production_35_44+frontline_service_35_44)/strtoi(hcw_age_data[6,6])
frontline_45_54_percentage = (frontline_management_45_54+frontline_sales_45_54+frontline_nat_res_45_54+
                                frontline_production_45_54+frontline_service_45_54)/strtoi(hcw_age_data[6,7])
frontline_55_64_percentage = (frontline_management_55_64+frontline_sales_55_64+frontline_nat_res_55_64+
                                frontline_production_55_64+frontline_service_55_64)/strtoi(hcw_age_data[6,8])
frontline_65_plus_percentage = (frontline_management_65_plus+frontline_sales_65_plus+frontline_nat_res_65_plus+
                                frontline_production_65_plus+frontline_service_65_plus)/strtoi(hcw_age_data[6,9])
#CDC essential non-frontline
#Production and Transportation
transportation = strtoi(hcw_age_data[539,2]) - strtoi(hcw_age_data[552,2]) - strtoi(hcw_age_data[547,2])

non_frontline_production_16_19 = strtoi(hcw_age_data[539,3]) - strtoi(hcw_age_data[552,3]) - strtoi(hcw_age_data[547,3])
non_frontline_production_20_24 = strtoi(hcw_age_data[539,4]) - strtoi(hcw_age_data[552,4]) - strtoi(hcw_age_data[547,4])
non_frontline_production_25_34 = strtoi(hcw_age_data[539,5]) - strtoi(hcw_age_data[552,5]) - strtoi(hcw_age_data[547,5])
non_frontline_production_35_44 = strtoi(hcw_age_data[539,6]) - strtoi(hcw_age_data[552,6]) - strtoi(hcw_age_data[547,6])
non_frontline_production_45_54 = strtoi(hcw_age_data[539,7]) - strtoi(hcw_age_data[552,7]) - strtoi(hcw_age_data[547,7])
non_frontline_production_55_64 = strtoi(hcw_age_data[539,8]) - strtoi(hcw_age_data[552,8]) - strtoi(hcw_age_data[547,8])
non_frontline_production_65_plus = strtoi(hcw_age_data[539,9]) - strtoi(hcw_age_data[552,9]) - strtoi(hcw_age_data[547,9])

non_frontline_production_percentage = transportation/strtoi(hcw_age_data[456,2])

#Service
food_service = strtoi(hcw_age_data[248,2])
grounds_maintenance = strtoi(hcw_age_data[262,2])
non_frontline_service_percentage = (food_service+grounds_maintenance)/strtoi(hcw_age_data[216,2])

non_frontline_service_16_19 = strtoi(hcw_age_data[248,3])+strtoi(hcw_age_data[262,3])
non_frontline_service_20_24 = strtoi(hcw_age_data[248,4])+strtoi(hcw_age_data[262,4])
non_frontline_service_25_34 = strtoi(hcw_age_data[248,5])+strtoi(hcw_age_data[262,5])
non_frontline_service_35_44 = strtoi(hcw_age_data[248,6])+strtoi(hcw_age_data[262,6])
non_frontline_service_45_54 = strtoi(hcw_age_data[248,7])+strtoi(hcw_age_data[262,7])
non_frontline_service_55_64 = strtoi(hcw_age_data[248,8])+strtoi(hcw_age_data[262,8])
non_frontline_service_65_plus = strtoi(hcw_age_data[248,9])+strtoi(hcw_age_data[262,9])

#Sales and office
office_admin_support = 0.5*strtoi(hcw_age_data[311,2])

non_frontline_sales_16_19 = 0.5*strtoi(hcw_age_data[311,3])
non_frontline_sales_20_24 = 0.5*strtoi(hcw_age_data[311,4])
non_frontline_sales_25_34 = 0.5*strtoi(hcw_age_data[311,5])
non_frontline_sales_35_44 = 0.5*strtoi(hcw_age_data[311,6])
non_frontline_sales_45_54 = 0.5*strtoi(hcw_age_data[311,7])
non_frontline_sales_55_64 = 0.5*strtoi(hcw_age_data[311,8])
non_frontline_sales_65_plus = 0.5*strtoi(hcw_age_data[311,9])

non_frontline_sales_percentage = office_admin_support/strtoi(hcw_age_data[291,2])

#Natural Resources
construction = strtoi(hcw_age_data[376,2])
IT_comms = strtoi(hcw_age_data[420,2]) + strtoi(hcw_age_data[444,2])
installation_repair = strtoi(hcw_age_data[417,2])

non_frontline_nat_res_16_19 = strtoi(hcw_age_data[376,3])+strtoi(hcw_age_data[420,3]) + strtoi(hcw_age_data[444,3])+strtoi(hcw_age_data[417,3])
non_frontline_nat_res_20_24 = strtoi(hcw_age_data[376,4])+strtoi(hcw_age_data[420,4]) + strtoi(hcw_age_data[444,4])+strtoi(hcw_age_data[417,4])
non_frontline_nat_res_25_34 = strtoi(hcw_age_data[376,5])+strtoi(hcw_age_data[420,5]) + strtoi(hcw_age_data[444,5])+strtoi(hcw_age_data[417,5])
non_frontline_nat_res_35_44 = strtoi(hcw_age_data[376,6])+strtoi(hcw_age_data[420,6]) + strtoi(hcw_age_data[444,6])+strtoi(hcw_age_data[417,6])
non_frontline_nat_res_45_54 = strtoi(hcw_age_data[376,7])+strtoi(hcw_age_data[420,7]) + strtoi(hcw_age_data[444,7])+strtoi(hcw_age_data[417,7])
non_frontline_nat_res_55_64 = strtoi(hcw_age_data[376,8])+strtoi(hcw_age_data[420,8]) + strtoi(hcw_age_data[444,8])+strtoi(hcw_age_data[417,8])
non_frontline_nat_res_65_plus = strtoi(hcw_age_data[376,9])+strtoi(hcw_age_data[420,9]) + strtoi(hcw_age_data[444,9])+strtoi(hcw_age_data[417,9])

non_frontline_nat_resouces_percentage = (construction+IT_comms+installation_repair)/strtoi(hcw_age_data[365,2])

#Management and professional
finance = strtoi(hcw_age_data[41,2])
media = strtoi(hcw_age_data[161,2])
legal = strtoi(hcw_age_data[143,2])
engineers = strtoi(hcw_age_data[88,2])

non_frontline_management_16_19 = strtoi(hcw_age_data[41,3])+strtoi(hcw_age_data[161,3])+strtoi(hcw_age_data[143,3])+strtoi(hcw_age_data[88,3])
non_frontline_management_20_24 = strtoi(hcw_age_data[41,4])+strtoi(hcw_age_data[161,4])+strtoi(hcw_age_data[143,4])+strtoi(hcw_age_data[88,4])
non_frontline_management_25_34 = strtoi(hcw_age_data[41,5])+strtoi(hcw_age_data[161,5])+strtoi(hcw_age_data[143,5])+strtoi(hcw_age_data[88,5])
non_frontline_management_35_44 = strtoi(hcw_age_data[41,6])+strtoi(hcw_age_data[161,6])+strtoi(hcw_age_data[143,6])+strtoi(hcw_age_data[88,6])
non_frontline_management_45_54 = strtoi(hcw_age_data[41,7])+strtoi(hcw_age_data[161,7])+strtoi(hcw_age_data[143,7])+strtoi(hcw_age_data[88,7])
non_frontline_management_55_64 = strtoi(hcw_age_data[41,8])+strtoi(hcw_age_data[161,8])+strtoi(hcw_age_data[143,8])+strtoi(hcw_age_data[88,8])
non_frontline_management_65_plus = strtoi(hcw_age_data[41,9])+strtoi(hcw_age_data[161,9])+strtoi(hcw_age_data[143,9])+strtoi(hcw_age_data[88,9])

non_frontline_management_percentage = (finance+media+legal+engineers)/strtoi(hcw_age_data[8,2])

CDC_non_frontline_total = transportation+food_service+construction+finance+
                          IT_comms+media+legal+engineers+grounds_maintenance+installation_repair+office_admin_support
CDC_non_frontline_percentage = CDC_non_frontline_total/strtoi(hcw_age_data[6,2])

CA_non_frontline_percentage = (sales_CA*non_frontline_sales_percentage+production_CA*non_frontline_production_percentage+
                                 natural_resources_CA*non_frontline_nat_resouces_percentage+service_CA*non_frontline_service_percentage+
                                 management_and_science_CA*non_frontline_management_percentage)/total_CA_workers

non_frontline_16_19_percentage = (non_frontline_management_16_19+non_frontline_sales_16_19+non_frontline_nat_res_16_19+
                                non_frontline_production_16_19+non_frontline_service_16_19)/strtoi(hcw_age_data[6,3])
non_frontline_20_24_percentage = (non_frontline_management_20_24+non_frontline_sales_20_24+non_frontline_nat_res_20_24+
                                non_frontline_production_20_24+non_frontline_service_20_24)/strtoi(hcw_age_data[6,4])
non_frontline_25_34_percentage = (non_frontline_management_25_34+non_frontline_sales_25_34+non_frontline_nat_res_25_34+
                                non_frontline_production_25_34+non_frontline_service_25_34)/strtoi(hcw_age_data[6,5])
non_frontline_35_44_percentage = (non_frontline_management_35_44+non_frontline_sales_35_44+non_frontline_nat_res_35_44+
                                non_frontline_production_35_44+non_frontline_service_35_44)/strtoi(hcw_age_data[6,6])
non_frontline_45_54_percentage = (non_frontline_management_45_54+non_frontline_sales_45_54+non_frontline_nat_res_45_54+
                                non_frontline_production_45_54+non_frontline_service_45_54)/strtoi(hcw_age_data[6,7])
non_frontline_55_64_percentage = (non_frontline_management_55_64+non_frontline_sales_55_64+non_frontline_nat_res_55_64+
                                non_frontline_production_55_64+non_frontline_service_55_64)/strtoi(hcw_age_data[6,8])
non_frontline_65_plus_percentage = (non_frontline_management_65_plus+non_frontline_sales_65_plus+non_frontline_nat_res_65_plus+
                                  non_frontline_production_65_plus+non_frontline_service_65_plus)/strtoi(hcw_age_data[6,9])

CA_workers_16_19 = total_CA_workers * strtoi(hcw_age_data[6,3])/strtoi(hcw_age_data[6,2])
CA_workers_20_24 = total_CA_workers * strtoi(hcw_age_data[6,4])/strtoi(hcw_age_data[6,2])
CA_workers_25_34 = total_CA_workers * strtoi(hcw_age_data[6,5])/strtoi(hcw_age_data[6,2])
CA_workers_35_44 = total_CA_workers * strtoi(hcw_age_data[6,6])/strtoi(hcw_age_data[6,2])
CA_workers_45_54 = total_CA_workers * strtoi(hcw_age_data[6,7])/strtoi(hcw_age_data[6,2])
CA_workers_55_64 = total_CA_workers * strtoi(hcw_age_data[6,8])/strtoi(hcw_age_data[6,2])
CA_workers_65_plus = total_CA_workers * strtoi(hcw_age_data[6,9])/strtoi(hcw_age_data[6,2])

#HCW's and education occupations
# Need to subtract archivists and curators from total 
frac_without_archivists = (census_occ_sex_alone[1,4] - census_occ_sex_alone[8,4])/census_occ_sex_alone[1,4]
educ_total = strtoi(hcw_age_data[149,2]) - strtoi(hcw_age_data[156,2]) - strtoi(hcw_age_data[157,2]) - strtoi(hcw_age_data[158,2])
educ_16_19 = (strtoi(hcw_age_data[149,3]) - strtoi(hcw_age_data[156,3]) - strtoi(hcw_age_data[157,3]) - strtoi(hcw_age_data[158,3]))/educ_total
educ_20_24 = (strtoi(hcw_age_data[149,4]) - strtoi(hcw_age_data[156,4]) - strtoi(hcw_age_data[157,4]) - strtoi(hcw_age_data[158,4]))/educ_total
educ_25_34 = (strtoi(hcw_age_data[149,5]) - strtoi(hcw_age_data[156,5]) - strtoi(hcw_age_data[157,5]) - strtoi(hcw_age_data[158,5]))/educ_total
educ_35_44 = (strtoi(hcw_age_data[149,6]) - strtoi(hcw_age_data[156,6]) - strtoi(hcw_age_data[157,6]) - strtoi(hcw_age_data[158,6]))/educ_total
educ_45_54 = (strtoi(hcw_age_data[149,7]) - strtoi(hcw_age_data[156,7]) - strtoi(hcw_age_data[157,7]) - strtoi(hcw_age_data[158,7]))/educ_total
educ_55_64 = (strtoi(hcw_age_data[149,8]) - strtoi(hcw_age_data[156,8]) - strtoi(hcw_age_data[157,8]) - strtoi(hcw_age_data[158,8]))/educ_total
educ_65_plus = (strtoi(hcw_age_data[149,9]) - strtoi(hcw_age_data[156,9]) - strtoi(hcw_age_data[157,9]) - strtoi(hcw_age_data[158,9]))/educ_total

hcw_total = strtoi(hcw_age_data[181,2]) + strtoi(hcw_age_data[217,2])
hcw_16_19 = (strtoi(hcw_age_data[181,3]) + strtoi(hcw_age_data[217,3]))/hcw_total
hcw_20_24 = (strtoi(hcw_age_data[181,4]) + strtoi(hcw_age_data[217,4]))/hcw_total
hcw_25_34 = (strtoi(hcw_age_data[181,5]) + strtoi(hcw_age_data[217,5]))/hcw_total
hcw_35_44 = (strtoi(hcw_age_data[181,6]) + strtoi(hcw_age_data[217,6]))/hcw_total
hcw_45_54 = (strtoi(hcw_age_data[181,7]) + strtoi(hcw_age_data[217,7]))/hcw_total
hcw_55_64 = (strtoi(hcw_age_data[181,8]) + strtoi(hcw_age_data[217,8]))/hcw_total
hcw_65_plus = (strtoi(hcw_age_data[181,9]) + strtoi(hcw_age_data[217,9]))/hcw_total


educ_ages = c(educ_16_19, educ_20_24, educ_25_34, educ_35_44, educ_45_54,
              educ_55_64, educ_65_plus)

hcw_ages = c(hcw_16_19, hcw_20_24, hcw_25_34, hcw_35_44, hcw_45_54, hcw_55_64,
             hcw_65_plus)
#SNF's, jails

counties_in_CA =  unique(census$NAME)
race_groups = c("White", "African American", "AIAN", "Asian Alone", 
                "Native Hawaiian And Other Pacific Islander Alone",
                "Some other race alone", "Two or more races")
ethnicity_list = c("Hispanic or Latino", "Not Hispanic or Latino")
sex_list = c("Male", "Female")

total_CA_pop = census_CA[3,4]
total_SNF_CA = census_SNF[1,4]
total_CA_prison_pop = census_CF[1,4]

prop_male_SNF = census_SNF[10,4]/total_SNF_CA
prop_female_SNF = census_SNF[19,4]/total_SNF_CA
male_18_24 = census_SNF[11,4]
male_25_34 = census_SNF[12,4]
male_35_44 = census_SNF[13,4]
male_45_54 = census_SNF[14,4]
male_55_64 = census_SNF[15,4]
male_65_74 = census_SNF[16,4]
male_75_84 = census_SNF[17,4]
male_85_plus = census_SNF[18,4]
female_18_24 = census_SNF[20,4]
female_25_34 = census_SNF[21,4]
female_35_44 = census_SNF[22,4]
female_45_54 = census_SNF[23,4]
female_55_64 = census_SNF[24,4]
female_65_74 = census_SNF[25,4]
female_75_84 = census_SNF[26,4]
female_85_plus = census_SNF[27,4]
num_white_SNF = census_SNF[36,4]
num_african_american_SNF = census_SNF[44,4]
num_AIAN_SNF = census_SNF[52,4]
num_asian_SNF = census_SNF[60,4]
num_native_hawaiian_pac_islander_SNF = census_SNF[68,4]
num_other_SNF = census_SNF[76,4]
num_multi_SNF = census_SNF[84,4]
num_hispanic_SNF = census_SNF[100,4]

male_15_17_CF = census_CF[12,4] 
male_18_24_CF = census_CF[13,4] 
male_25_34_CF = census_CF[14,4] 
male_35_44_CF = census_CF[15,4] 
male_45_54_CF = census_CF[16,4] 
male_55_64_CF = census_CF[17,4] 
male_65_74_CF = census_CF[18,4]
male_75_84_CF = census_CF[19,4]
male_85_plus_CF = census_CF[20,4]

female_15_17_CF = census_CF[22,4] 
female_18_24_CF = census_CF[23,4] 
female_25_34_CF = census_CF[24,4] 
female_35_44_CF = census_CF[25,4] 
female_45_54_CF = census_CF[26,4] 
female_55_64_CF = census_CF[27,4] 
female_65_74_CF = census_CF[28,4]
female_75_84_CF = census_CF[29,4]
female_85_plus_CF = census_CF[30,4]

num_white_CF = census_CF[39,4]
num_african_american_CF = census_CF[47,4]
num_AIAN_CF = census_CF[55,4]
num_asian_CF = census_CF[63,4]
num_native_hawaiian_pac_islander_CF = census_CF[71,4]
num_other_CF = census_CF[79,4]
num_multi_CF= census_CF[87,4]
num_hispanic_CF = census_CF[103,4]


num_white_CA = census_CA[333,4] + census_CA[343,4]
num_african_american_CA = census_CA[334,4] + census_CA[344,4]
num_AIAN_CA = census_CA[335,4] + census_CA[345,4]
num_asian_CA = census_CA[336,4] + census_CA[346,4]
num_native_hawaiian_pac_islander_CA = census_CA[337,4] + census_CA[347,4]
num_other_CA = census_CA[338,4] + census_CA[348,4]
num_multi_CA = census_CA[339,4] + census_CA[349,4]
num_hispanic_CA = census_CA[342,4]
num_not_hispanic_CA = total_CA_pop - num_hispanic_CA

male_0_4_CA = census_CA[5,4]
male_5_9_CA = census_CA[6,4]
male_10_14_CA = census_CA[7,4]
male_15_17_CA =  census_CA[8,4]
male_18_24_CA =  census_CA[9,4] + census_CA[10,4] + census_CA[11,4] + census_CA[12,4]
male_25_34_CA =  census_CA[13,4] + census_CA[14,4]
male_35_44_CA =  census_CA[15,4] + census_CA[16,4]
male_45_54_CA =  census_CA[17,4] + census_CA[18,4]
male_55_64_CA =  census_CA[19,4] + census_CA[20,4] + census_CA[21,4]
male_65_74_CA =  census_CA[22,4] + census_CA[23,4] + census_CA[24,4]
male_75_84_CA =  census_CA[25,4] + census_CA[26,4]
male_85_plus_CA =  census_CA[27,4]

female_0_4_CA = census_CA[24+5,4]
female_5_9_CA = census_CA[24+6,4]
female_10_14_CA = census_CA[24+7,4]
female_15_17_CA =  census_CA[24+8,4]
female_18_24_CA =  census_CA[24+9,4] + census_CA[24+10,4] + census_CA[24+11,4] + census_CA[24+12,4]
female_25_34_CA =  census_CA[24+13,4] + census_CA[24+14,4]
female_35_44_CA =  census_CA[24+15,4] + census_CA[24+16,4]
female_45_54_CA =  census_CA[24+17,4] + census_CA[24+18,4]
female_55_64_CA =  census_CA[24+19,4] + census_CA[24+20,4] + census_CA[24+21,4]
female_65_74_CA =  census_CA[24+22,4] + census_CA[24+23,4] + census_CA[24+24,4]
female_75_84_CA =  census_CA[24+25,4] + census_CA[24+26,4]
female_85_plus_CA =  census_CA[24+27,4]




CA_age_totals = c(male_0_4_CA, male_5_9_CA, male_10_14_CA, male_15_17_CA, male_18_24_CA,
                  male_25_34_CA, male_35_44_CA, male_45_54_CA, male_55_64_CA, male_65_74_CA,
                  male_75_84_CA, male_85_plus_CA, female_0_4_CA, female_5_9_CA, female_10_14_CA,
                  female_15_17_CA, female_18_24_CA, female_25_34_CA, female_35_44_CA,
                  female_45_54_CA, female_55_64_CA, female_65_74_CA, female_75_84_CA,
                  female_85_plus_CA)

frontline_hash = hashmap('a',1)
frontline_hash$insert(15, 0.5*CA_workers_16_19*0.5*frontline_16_19_percentage/(male_15_17_CA+female_15_17_CA))
frontline_hash$insert(18, 0.5*CA_workers_16_19*0.5*frontline_16_19_percentage/(male_15_17_CA+female_15_17_CA)+
                            CA_workers_20_24*frontline_20_24_percentage/(male_18_24_CA+female_18_24_CA))
frontline_hash$insert(25, CA_workers_25_34*frontline_25_34_percentage/(male_25_34_CA+female_25_34_CA))
frontline_hash$insert(35, CA_workers_35_44*frontline_35_44_percentage/(male_35_44_CA+female_35_44_CA))
frontline_hash$insert(45, CA_workers_45_54*frontline_45_54_percentage/(male_45_54_CA+female_45_54_CA))
frontline_hash$insert(55, CA_workers_55_64*frontline_55_64_percentage/(male_55_64_CA+female_55_64_CA))
frontline_hash$insert(65, 0.75*CA_workers_65_plus*0.75*frontline_65_plus_percentage/(male_65_74_CA+female_65_74_CA))
frontline_hash$insert(75, 0.25*CA_workers_65_plus*0.25*frontline_65_plus_percentage/(male_75_84_CA+female_75_84_CA))
frontline_hash$insert(85, 0)

non_frontline_hash = hashmap('a', 1)
non_frontline_hash$insert(15, 0.5*CA_workers_16_19*0.5*non_frontline_16_19_percentage/(male_15_17_CA+female_15_17_CA))
non_frontline_hash$insert(18, 0.5*CA_workers_16_19*0.5*non_frontline_16_19_percentage/(male_15_17_CA+female_15_17_CA)+
                        CA_workers_20_24*non_frontline_20_24_percentage/(male_18_24_CA+female_18_24_CA))
non_frontline_hash$insert(25, CA_workers_25_34*non_frontline_25_34_percentage/(male_25_34_CA+female_25_34_CA))
non_frontline_hash$insert(35, CA_workers_35_44*non_frontline_35_44_percentage/(male_35_44_CA+female_35_44_CA))
non_frontline_hash$insert(45, CA_workers_45_54*non_frontline_45_54_percentage/(male_45_54_CA+female_45_54_CA))
non_frontline_hash$insert(55, CA_workers_55_64*non_frontline_55_64_percentage/(male_55_64_CA+female_55_64_CA))
non_frontline_hash$insert(65, 0.75*CA_workers_65_plus*0.75*non_frontline_65_plus_percentage/(male_65_74_CA+female_65_74_CA))
non_frontline_hash$insert(75, 0.25*CA_workers_65_plus*0.25*non_frontline_65_plus_percentage/(male_75_84_CA+female_75_84_CA))
non_frontline_hash$insert(85, 0)

ALF_hash = hashmap('a',1)
ALF_hash$insert(15, 0)
ALF_hash$insert(18, 0)
ALF_hash$insert(25, 0)
ALF_hash$insert(35, 0)
ALF_hash$insert(45, 0)
ALF_hash$insert(55, 0.03*150000/(male_55_64_CA +female_55_64_CA))
ALF_hash$insert(65, 0.104*150000/(male_65_74_CA +female_65_74_CA))
ALF_hash$insert(75, 0.319*150000/(male_75_84_CA +female_75_84_CA))
ALF_hash$insert(85, 0.548*150000/(male_85_plus_CA +female_85_plus_CA))


CA_race_totals = c(num_white_CA, num_african_american_CA, num_AIAN_CA,
                   num_asian_CA, num_native_hawaiian_pac_islander_CA, 
                   num_other_CA, num_multi_CA)

SNF_race_totals = c(num_white_SNF, num_african_american_SNF, num_AIAN_SNF,
                    num_asian_SNF, num_native_hawaiian_pac_islander_SNF, 
                    num_other_SNF, num_multi_SNF)

CF_race_totals = c(num_white_CF, num_african_american_CF, num_AIAN_CF,
                    num_asian_CF, num_native_hawaiian_pac_islander_CF, 
                    num_other_CF, num_multi_CF)




set.seed(123)
stdev_hash = hashmap('a',1)
asthma_prob_dummy_hash = hashmap('a',1)
diabetes_prob_dummy_hash = hashmap('a',1)
smoker_prob_dummy_hash = hashmap('a',1)
heart_disease_prob_dummy_hash = hashmap('a',1)
hf_prob_dummy_hash = hashmap('a',1)
htn_prob_dummy_hash = hashmap('a',1)
bmi_dummy_hash = hashmap('a',1)

CF_dummy_hash = hashmap('a',1)
SNF_dummy_hash = hashmap('a',1)
homeless_dummy_hash = hashmap('a',1)
hcw_dummy_hash = hashmap('a', 1)
educ_dummy_hash = hashmap('a',1)


county_SNF = vector()
age_low = vector()
age_high = vector()
sex_SNF = vector()
race_overall = vector()
ethn_overall = vector()
eth_SNF = vector()
SNF_props = vector()
CF_props = vector()
educ_props = vector()
hcw_props = vector()
homeless_props = vector()
SNF_ages = c(0,4, 5,9, 10,14, 15,17,18,24,25,34,35,44,45,54,55,64,65,74,75,84,85,100)

# #Load in CHIS data tables
# chis_overall = read.csv(file = "20201115/CA_wide.csv")
# chis_ab40 = read.csv(file = "20201115/ab40_df.csv")
# chis_ab34 = read.csv(file = "20201115/ab34_df.csv")
# chis_ab29 = read.csv(file = "20201115/ab29_df.csv")
# chis_diabetes = read.csv(file = "20201115/diabetes_df.csv")
# chis_smkcur = read.csv(file = "20201115/smkcur_df.csv")
# chis_ab52 = read.csv(file = "20201115/ab52_df.csv")
# 
# chis_bmi_age = read.csv(file ="20201115/bmi_age_df.csv")
# chis_bmi_sex = read.csv(file = "20201115/bmi_sex_df.csv")
# chis_bmi_race = read.csv(file = "20201115/bmi_race_df.csv")
# chis_bmi_eth = read.csv(file = "20201115/bmi_eth_df.csv")
# 
# CA_bmi_mean = read.csv(file = "20201115/bmi_overall_df.csv")[1,2]
# CA_bmi_std = read.csv(file = "20201115/bmi_overall_df.csv")[1,4] * CA_bmi_mean
# 
# chis_age_county_ab34 = read.csv(file = "20201118/ab34_age_county_df.csv")
# chis_age_county_ab40 = read.csv(file = "20201118/ab40_age_county_df.csv")
# chis_age_county_ab29 = read.csv(file = "20201118/ab29_age_county_df.csv")
# chis_age_county_smkcur = read.csv(file = "20201118/smkcur_age_county_df.csv")
# chis_age_county_diabetes = read.csv(file = "20201118/diabetes_age_county_df.csv")
# chis_age_county_bmi = read.csv(file = "20201115/bmi_age_county_df.csv")
# chis_age_county_bmi_totals = read.csv(file = "20201115/bmi_age_county_totals.csv")
# 
# chis_age_county_ab34 = chis_age_county_ab34[order(chis_age_county_ab34$srage_code),]
# chis_age_county_ab40 = chis_age_county_ab40[order(chis_age_county_ab40$srage_code),]
# chis_age_county_ab29 = chis_age_county_ab29[order(chis_age_county_ab29$srage_code),]
# chis_age_county_diabetes = chis_age_county_diabetes[order(chis_age_county_diabetes$srage_code),]
# chis_age_county_smkcur = chis_age_county_smkcur[order(chis_age_county_smkcur$srage_code),]
# chis_age_county_bmi = chis_age_county_bmi[order(chis_age_county_bmi$srage_code),]
# chis_age_county_bmi_totals = chis_age_county_bmi_totals[order(chis_age_county_bmi_totals$srage_code),]
# 
# chis_age_county_ab34[chis_age_county_ab34 == "San Bernandino County, California"] = "San Bernardino County, California"
# chis_age_county_ab40[chis_age_county_ab40 == "San Bernandino County, California"] = "San Bernardino County, California"
# chis_age_county_ab29[chis_age_county_ab29 == "San Bernandino County, California"] = "San Bernardino County, California"
# chis_age_county_diabetes[chis_age_county_diabetes == "San Bernandino County, California"] = "San Bernardino County, California"
# chis_age_county_smkcur[chis_age_county_smkcur == "San Bernandino County, California"] = "San Bernardino County, California"
# chis_age_county_bmi[chis_age_county_bmi == "San Bernandino County, California"] = "San Bernardino County, California"
# chis_age_county_bmi_totals[chis_age_county_bmi_totals == "San Bernandino County, California"] = "San Bernardino County, California"
# 
# #overall_bmi_mean = (chis_bmi_sex[1,3]*chis_bmi_sex[1,7] + chis_bmi_sex[2,3]*chis_bmi_sex[2,7])/(chis_bmi_sex[1,7]+chis_bmi_sex[2,7])
# #overall_bmi_std = sqrt((chis_bmi_sex[1,3]*chis_bmi_sex[1,9])^2 + (chis_bmi_sex[2,3]*chis_bmi_sex[2,9])^2)
# ab40_overall = chis_overall[2,2]
# ab34_overall = chis_overall[12,2]
# ab29_overall = chis_overall[16,2]
# ab52_overall = chis_overall[14,2]
# diabetes_overall = chis_overall[6,2]
# smkcur_overall = chis_overall[9,2]

for (county in counties_in_CA) {
  occ_filtered = filter(census_occ, NAME == county)
  pop_filtered = filter(census, NAME == county)
  county_population = pop_filtered[3,4]
  num_male_educ = occ_filtered[1,4] * frac_without_archivists
  num_male_hcw = occ_filtered[2,4] + occ_filtered[5,4]
  num_female_educ = occ_filtered[6,4] * frac_without_archivists
  num_female_hcw = occ_filtered[7,4] + occ_filtered[10,4]
  
  prop_white_race = (pop_filtered[333,4] + pop_filtered[343,4])/county_population
  prop_african_american_race = (pop_filtered[334,4] + pop_filtered[344,4])/county_population
  prop_AIAN_race = (pop_filtered[335,4] + pop_filtered[345,4])/county_population
  prop_asian_race = (pop_filtered[336,4] + pop_filtered[346,4])/county_population
  prop_NHPI_race = (pop_filtered[337,4] + pop_filtered[347,4])/county_population
  prop_other_race = (pop_filtered[338,4] + pop_filtered[348,4])/county_population
  prop_multi_race = (pop_filtered[339,4] + pop_filtered[349,4])/county_population
  
  county_race_props = c(prop_white_race, prop_african_american_race, prop_AIAN_race,
                        prop_asian_race, prop_NHPI_race, prop_other_race, prop_multi_race)
  
  num_homeless_county = filter(homeless_data, CA_county == county)[1,18]
  frac_homeless_under_18_county = filter(homeless_data, CA_county == county)[1, 19]/num_homeless_county
  frac_homeless_18_24_county = filter(homeless_data, CA_county == county)[1, 20]/num_homeless_county
  frac_homeless_over_24_county = filter(homeless_data, CA_county == county)[1, 21]/num_homeless_county
  age_differential = (frac_homeless_under_18_county + frac_homeless_18_24_county) - prop_homeless_under_18 - prop_homeless_18_24
  frac_homeless_25_54 = prop_homeless_25_54 - age_differential/3
  frac_homeless_55_61 = prop_homeless_55_61 - age_differential/3
  frac_homeless_62_plus = prop_homeless_62_plus - age_differential/3
  
  frac_homeless_county_white = filter(homeless_data, CA_county == county)[1, 28]/num_homeless_county
  frac_homeless_county_african_american = filter(homeless_data, CA_county == county)[1, 29]/num_homeless_county
  frac_homeless_county_asian = filter(homeless_data, CA_county == county)[1, 30]/num_homeless_county
  frac_homeless_county_AIAN = filter(homeless_data, CA_county == county)[1, 31]/num_homeless_county
  frac_homeless_county_NHPI = filter(homeless_data, CA_county == county)[1, 32]/num_homeless_county
  #frac_homeless_county_other = 
  #  max((filter(homeless_data, CA_county == county)[1, 33] - prop_multi_race*num_homeless_county)/num_homeless_county, 0)
  #frac_homeless_county_multi = prop_multi_race
  frac_homeless_county_other = prop_other_race/(prop_other_race +prop_multi_race)*filter(homeless_data, CA_county == county)[1, 33]/num_homeless_county
  frac_homeless_county_multi = prop_multi_race/(prop_other_race +prop_multi_race)*filter(homeless_data, CA_county == county)[1, 33]/num_homeless_county
  
  county_homeless_age_props = c(0,0,0, frac_homeless_under_18_county, frac_homeless_18_24_county,
                                frac_homeless_25_54/3, frac_homeless_25_54/3, frac_homeless_25_54/3,
                                frac_homeless_55_61, frac_homeless_62_plus*2/3, frac_homeless_62_plus/3, 0)
  
  county_homeless_race_props = c(frac_homeless_county_white, frac_homeless_county_african_american,
                                 frac_homeless_county_AIAN, frac_homeless_county_asian, frac_homeless_county_NHPI,
                                 frac_homeless_county_other, frac_homeless_county_multi)
  #num_homeless_county = num_homeless_CA * pop_filtered[3,4] /total_CA_pop
  
  county_male_0_4 = pop_filtered[5,4]
  county_male_5_9 = pop_filtered[6,4]
  county_male_10_14 = pop_filtered[7,4]
  county_male_15_17 = pop_filtered[8,4]
  county_male_18_24 = pop_filtered[9,4] + pop_filtered[10,4] + pop_filtered[11,4] + pop_filtered[12,4]
  county_male_25_34 = pop_filtered[13,4] + pop_filtered[14,4]
  county_male_35_44 = pop_filtered[15,4] + pop_filtered[16,4]
  county_male_45_54 = pop_filtered[17,4] + pop_filtered[18,4]
  county_male_55_64 = pop_filtered[19,4] + pop_filtered[20,4] + pop_filtered[21,4]
  county_male_65_74 = pop_filtered[22,4] + pop_filtered[23,4] + pop_filtered[24,4]
  county_male_75_84 = pop_filtered[25,4] + pop_filtered[26,4]
  county_male_85_plus = pop_filtered[27,4]
  male_ages = c(county_male_0_4, county_male_5_9, county_male_10_14,
                county_male_15_17, county_male_18_24, county_male_25_34,
                county_male_35_44, county_male_45_54, county_male_55_64,
                county_male_65_74, county_male_75_84, county_male_85_plus)
  
  delta = 24
  county_female_0_4 = pop_filtered[5+delta,4]
  county_female_5_9 = pop_filtered[6+delta,4]
  county_female_10_14 = pop_filtered[7+delta,4]
  county_female_15_17 = pop_filtered[8+delta,4]
  county_female_18_24 = pop_filtered[9+delta,4] + pop_filtered[10+delta,4] + pop_filtered[11+delta,4] + pop_filtered[12+delta,4]
  county_female_25_34 = pop_filtered[13+delta,4] + pop_filtered[14+delta,4]
  county_female_35_44 = pop_filtered[15+delta,4] + pop_filtered[16+delta,4]
  county_female_45_54 = pop_filtered[17+delta,4] + pop_filtered[18+delta,4]
  county_female_55_64 = pop_filtered[19+delta,4] + pop_filtered[20+delta,4] + pop_filtered[21+delta,4]
  county_female_65_74 = pop_filtered[22+delta,4] + pop_filtered[23+delta,4] + pop_filtered[24+delta,4]
  county_female_75_84 = pop_filtered[25+delta,4] + pop_filtered[26+delta,4]
  county_female_85_plus = pop_filtered[27+delta,4]
  
  female_ages = c(county_female_0_4, county_female_5_9, county_female_10_14,
                county_female_15_17, county_female_18_24, county_female_25_34,
                county_female_35_44, county_female_45_54, county_female_55_64,
                county_female_65_74, county_female_75_84, county_female_85_plus)
  
  # # filter county data table from chis to only look at this county
  # #age_county_strat = filter(chis_age_county, county = "County")
  # filt_ab34_age_county = filter(chis_age_county_ab34, srcnty_code == county)
  # filt_ab40_age_county = filter(chis_age_county_ab40, srcnty_code == county)
  # filt_ab29_age_county = filter(chis_age_county_ab29, srcnty_code == county)
  # filt_diabetes_age_county = filter(chis_age_county_diabetes, srcnty_code == county)
  # filt_smkcur_age_county = filter(chis_age_county_smkcur, srcnty_code == county)
  # filt_bmi_age_county = filter(chis_age_county_bmi, srcnty_code == county)
  # filt_bmi_age_county_totals = filter(chis_age_county_bmi_totals, srcnty_code == county)
  # 
  # filt_ab34_age_county[filt_ab34_age_county == 'na'] = 0
  # filt_ab34_age_county[is.na(filt_ab34_age_county)] = 0
  # 
  # filt_ab40_age_county[filt_ab40_age_county == 'na'] = 0
  # filt_ab40_age_county[is.na(filt_ab40_age_county)] = 0
  # 
  # filt_ab29_age_county[filt_ab29_age_county == 'na'] = 0
  # filt_ab29_age_county[is.na(filt_ab29_age_county)] = 0
  # 
  # filt_diabetes_age_county[filt_diabetes_age_county == 'na'] = 0
  # filt_diabetes_age_county[is.na(filt_diabetes_age_county)] = 0
  # 
  # filt_smkcur_age_county[filt_smkcur_age_county == 'na'] = 0
  # filt_smkcur_age_county[is.na(filt_smkcur_age_county)] = 0
  # 
  # filt_bmi_age_county[filt_bmi_age_county == 'na'] = 0
  # filt_bmi_age_county[is.na(filt_bmi_age_county)] = 0
  # 
  # filt_bmi_age_county_totals[filt_bmi_age_county_totals == 'na'] = 0
  # filt_bmi_age_county_totals[is.na(filt_bmi_age_county_totals)] = 0
  
  for (j in c(1: length(race_groups))) {
    # #filter race table
    # race_char = race_groups[j]
    # race_index_chis = 0
    # if (race_char == "White") {
    #   race_index_chis = 19
    # } else if (race_char == "African American") {
    #   race_index_chis = 13
    # } else if (race_char == "AIAN") {
    #   race_index_chis = 14
    # } else if (race_char == "Asian Alone") {
    #   race_index_chis = 15
    # } else if (race_char == "Native Hawaiian And Other Pacific Islander Alone") {
    #   race_index_chis = 18
    # } else if (race_char == "Some other race alone") {
    #   race_index_chis = 17
    # } else {
    #   race_index_chis = 16
    # }
    # ab40_race = chis_ab40[race_index_chis, 4]
    # ab34_race = chis_ab34[race_index_chis, 4]
    # ab29_race = chis_ab29[race_index_chis, 4]
    # ab52_race = chis_ab52[race_index_chis, 4]
    # diabetes_race = chis_diabetes[race_index_chis, 4]
    # smkcur_race = chis_smkcur[race_index_chis, 3]
    # bmi_race_mean = chis_bmi_race[race_index_chis-12, 3]
    # bmi_race_std = chis_bmi_race[race_index_chis-12, 9] * bmi_race_mean
    
    for (sex in sex_list) {
      # #filter sex table
      # #filt_chis_sex = chis_sex[sex]
      # sex_index_chis = 1
      # if (sex == "Male") { 
      #   bmi_sex_mean = chis_bmi_sex[2,3]
      #   bmi_sex_std = bmi_sex_mean*chis_bmi_sex[2,9] 
      #   sex_index_chis = 2
      # } else {
      #   bmi_sex_mean = chis_bmi_sex[1,3]
      #   bmi_sex_std = bmi_sex_mean*chis_bmi_sex[1,9] 
      # }
      # ab40_sex = chis_ab40[sex_index_chis, 4]
      # ab34_sex = chis_ab34[sex_index_chis, 4]
      # ab29_sex = chis_ab29[sex_index_chis, 4]
      # ab52_sex = chis_ab52[sex_index_chis, 4]
      # diabetes_sex = chis_diabetes[sex_index_chis, 4]
      # smkcur_sex = chis_smkcur[sex_index_chis, 3]
      
      for (eth in ethnicity_list) {
        #filter eth table
        # if (eth == "Hispanic or Latino") {
        #   #filt_eth = chis_eth[,3]
        #   bmi_eth_mean = chis_bmi_eth[1,3]
        #   bmi_eth_std = chis_bmi_eth[1,9]*bmi_eth_mean
        # } else {
        #   #filt_eth = chis_eth[, 4]
        #   bmi_eth_mean = chis_bmi_eth[2,3]
        #   bmi_eth_std = chis_bmi_eth[2,9]*bmi_eth_mean
        # }
        for (i in c(1:12)) {
          # bmi_mean = 23
          # bmi_std = 2
          # ab40 = 0.071 #wtvr CA avg is
          # ab34 = 0 #assume no heart disease in <18
          # smkcur = 0.02 #5.8% of hs students smoke, roughly 1/3 of the 0-18 pop so assume 2%
          # diabetes = 0.0025 #Google search, assuming 0.25% diabetes in <18
          # ab29 = 0.03 #assume 3% htn for <18 
          # ab52 = 0 #assume no heart failure for < 18
          # if (i > 4) {
          #   #bmi_age_mean = filt_bmi_age_county[i-4,4]
          #   #bmi_age_std = filt_bmi_age_county[i-4,6]*bmi_age_mean
          #   
          #   bmi_age_mean = chis_bmi_age[i-4,3]
          #   bmi_age_std = bmi_age_mean *chis_bmi_age[i-4,9]
          #  
          #   ab52_age = chis_ab52[i,4]
          #   
          #   ab40_age = chis_ab40[i, 4]
          #   ab34_age = chis_ab34[i,4]
          #   #ab29_age = chis_ab29[i,4]
          #   ab29_age = chis_ab29[i, 4]
          #   diabetes_age = chis_diabetes[i,4]
          #   smkcur_age = chis_smkcur[i,3]
          #   #race, sex, eth adjustments
          #   
          #   #alternative age group mappings when county filtering:
          #   #1 65-74
          #   #2 75-84
          #   #3 55-64
          #   #if (filt_bmi_age_county_totals[i-4,5] == "na" | filt_bmi_age_county_totals[i-4,5] < 10 |
          #    #   is.na(filt_bmi_age_county_totals[i-4,5])) {
          #   #  bmi_age_mean = chis_bmi_age[i-4,3]
          #   #  bmi_age_std = bmi_age_mean *chis_bmi_age[i-4,9]
          #   #} 
          #   if (as.numeric(filt_ab29_age_county[i-4, 14]) + as.numeric(filt_ab29_age_county[i-4, 15]) > 10) {
          #     ab29_age = filt_ab29_age_county[i-4,4]
          #   }
          #   
          #   if (as.numeric(filt_ab34_age_county[i-4, 14]) + as.numeric(filt_ab34_age_county[i-4, 15]) > 10) {
          #     ab34_age = filt_ab34_age_county[i-4,5]
          #   }
          #   
          #   if (as.numeric(filt_ab40_age_county[i-4, 14]) + as.numeric(filt_ab40_age_county[i-4, 15]) > 10) {
          #     ab40_age = filt_ab40_age_county[i-4,5]
          #   }
          #   
          #   if (as.numeric(filt_diabetes_age_county[i-4, 14]) + as.numeric(filt_diabetes_age_county[i-4, 15]) > 10) {
          #     diabetes_age = filt_diabetes_age_county[i-4,5]
          #   }
          #   
          #   if (as.numeric(filt_smkcur_age_county[i-4, 14]) + as.numeric(filt_smkcur_age_county[i-4, 15]) > 10) {
          #     smkcur_age = filt_smkcur_age_county[i-4,5]
          #   }
          #   
          #   if (as.numeric(filt_bmi_age_county_totals[i-4, 5]) > 10) {
          #     bmi_age_mean = filt_bmi_age_county[i-4,4]
          #     bmi_age_std = filt_bmi_age_county[i-4,6]*bmi_age_mean
          #   }
          # 
          #   ## add some more conditions to deal with na's, NA's, and <10 vals for other sets
          #   ab40 = ab40_age * (ab40_sex/ab40_overall) * (ab40_race/ab40_overall)
          #   ab34 = ab34_age * (ab34_sex/ab34_overall) *(ab34_race/ab34_overall)
          #   ab29 = ab29_age * (ab29_sex/ab29_overall) *(ab29_race/ab29_overall)
          #   ab52 = ab52_age * (ab52_sex/ab52_overall) *(ab52_race/ab52_overall)
          #   diabetes = diabetes_age * (diabetes_sex/diabetes_overall) *(diabetes_race/diabetes_overall)
          #   smkcur = smkcur_age * (smkcur_sex/smkcur_overall) *(smkcur_race/smkcur_overall)
          #   bmi_mean = bmi_age_mean * (bmi_sex_mean/CA_bmi_mean) *(bmi_race_mean/CA_bmi_mean)
          #   bmi_std = bmi_age_std * (bmi_sex_std/CA_bmi_std) *(bmi_race_std/CA_bmi_std)
          # 
          #   
          #   ab40 = min(ab40, 1)
          #   diabetes = min(diabetes,1)
          #   smkcur = min(smkcur, 1)
          #   ab34 = min(ab34, 1)
          #   ab29 = min(ab29, 1)
          #   ab52 = min(ab52, 1)
          #   
          #   if (is.na(ab40)) {
          #     print(county)
          #     print(eth)
          #     print("ab40")
          #     print(i)
          #   }
          #   if (is.na(ab34)) {
          #     print(county)
          #     print(eth)
          #     print("ab34")
          #     print(i)
          #   }
          #   if (is.na(ab29)) {
          #     print(county)
          #     print(eth)
          #     print("ab29")
          #     print(i)
          #   }
          #   if (is.na(ab52)) {
          #     print(county)
          #     print(eth)
          #     print("ab52")
          #     print(i)
          #   }
          #   if (is.na(diabetes)) {
          #     print(county)
          #     print(eth)
          #     print("diabetes")
          #     print(i)
          #     print(x2)
          #   }
          #   if (is.na(smkcur)) {
          #     print(county)
          #     print(eth)
          #     print("smoker")
          #     print(i)
          #   }
          #   if (is.na(bmi_mean) | is.na(bmi_std)) {
          #     print(county)
          #     print(eth)
          #     print("bmi")
          #     print(i)
          #     print(bmi_age_mean)
          #     print(bmi_race_mean)
          #     print(bmi_sex_mean)
          #     print(bmi_mean)
          #     print(bmi_std)
          #   }
          # } 
          
          demog_pop = male_ages[i]
          plus_65_bin = male_ages[10] + male_ages[11]
          age_low = c(age_low, SNF_ages[i+i-1])
          age_high = c(age_high, SNF_ages[i+i])
          sex_SNF = c(sex_SNF, sex)
          race_overall = c(race_overall, race_groups[j])
          eth_SNF = c(eth_SNF, eth)
          county_SNF = c(county_SNF, county)
          prop_homeless_sex = prop_male_homeless_CA/CA_prop_male
          num_educ = num_male_educ
          num_hcw = num_male_hcw
          #What proportion of people with this age, sex, race, and ethnicity are in SNFs?
          offset = 0
          k = 0
          if (sex == "Female") {
            offset = 9
            demog_pop = female_ages[i]
            plus_65_bin = female_ages[10] + female_ages[11]
            k = 12
            num_hcw = num_female_hcw
            num_educ = num_female_educ
            prop_homeless_sex = prop_female_homeless_CA/CA_prop_female
          }
          #What proportion of individuals in this age group and county are
          # healthcare workers or educational occupations?
          hcw_prop = 0
          educ_prop = 0
          if (SNF_ages[i+i] > 17 & SNF_ages[i+i-1] <85 ) {
            if (SNF_ages[i+i] == 24) {
              frac_educ = educ_ages[1]/2 + educ_ages[2]
              frac_hcw = hcw_ages[1]/2 + hcw_ages[2]
              educ_prop = num_educ*frac_educ/demog_pop
              hcw_prop = num_hcw * frac_hcw/demog_pop
              
            } else if (SNF_ages[i+i] == 74) {
              frac_educ = educ_ages[7]*2/3
              frac_hcw = hcw_ages[7]*2/3
              educ_prop = num_educ*frac_educ/demog_pop
              hcw_prop = num_hcw * frac_hcw/demog_pop
              
            } else if (SNF_ages[i+i] == 84) {
              frac_educ = educ_ages[7]*1/3
              frac_hcw = hcw_ages[7]*1/3
              educ_prop = num_educ*frac_educ/demog_pop
              hcw_prop = num_hcw * frac_hcw/demog_pop
              
            } else {
              frac_educ = educ_ages[i-3]
              frac_hcw = hcw_ages[i-3]
              educ_prop = num_educ*frac_educ/demog_pop
              hcw_prop = num_hcw * frac_hcw/demog_pop
            }
          }
          
          educ_props = c(educ_props, educ_prop)
          hcw_props = c(hcw_props, hcw_prop)
          #prop_homeless_age_group = CA_homeless_props[i]*num_homeless_CA/(CA_age_totals[i] + CA_age_totals[i+12])
          #prop_homeless_age_group = CA_homeless_props[i]*num_homeless_county/(male_ages[i] + female_ages[i])
          prop_homeless_age_group = county_homeless_age_props[i]*num_homeless_county/(male_ages[i] + female_ages[i])
            
         #homeless_race_prop = CA_homeless_props_race[j]/(CA_race_totals[j]/total_CA_pop)
          rp = county_race_props[j]
          if (county_race_props[j] == 0) {
            rp = CA_race_totals[j]/total_CA_pop
          }
          homeless_race_prop = county_homeless_race_props[j]/rp
          homeless_prop = prop_homeless_age_group * homeless_race_prop *prop_homeless_sex
          if (homeless_race_prop < 0 | homeless_prop < 0 | prop_homeless_age_group < 0) {
            print(county)
            print(prop_homeless_age_group)
            print(homeless_race_prop)
          }
          if (homeless_prop > 0.2) {
            print(homeless_prop)
            print(county)
            print(i)
            print(prop_homeless_age_group)
            print(race_groups[j])
            print(homeless_race_prop)
            #print(SNF_ages[i+i-1])
          }
          
          homeless_props = c(homeless_props, homeless_prop)
          if (eth == "Hispanic or Latino") {
            eth_fraction = num_hispanic_SNF/total_SNF_CA
            eth_fraction = eth_fraction /(num_hispanic_CA/total_CA_pop)
            eth_frac_CF = (num_hispanic_CF/total_CA_prison_pop)/(num_hispanic_CA/total_CA_pop)
          } else {
            eth_fraction = ((total_SNF_CA - num_hispanic_SNF)/total_SNF_CA)/(num_not_hispanic_CA/total_CA_pop)
            eth_frac_CF = ((total_CA_prison_pop - num_hispanic_CF)/total_CA_prison_pop)/(num_not_hispanic_CA/total_CA_pop)
          }
          SNF_race_fraction = (SNF_race_totals[j]/total_SNF_CA)/(CA_race_totals[j]/total_CA_pop)
          CF_race_fraction = (CF_race_totals[j]/total_CA_prison_pop)/(CA_race_totals[j]/total_CA_pop)
          
          SNF_prop = (census_SNF[6+i+offset,4] * eth_fraction * SNF_race_fraction)/CA_age_totals[k+i]
          CF_prop = (census_CF[8+i+offset,4] * eth_frac_CF * CF_race_fraction)/CA_age_totals[k+i]
          #if ((SNF_prop > 1) | (CF_prop > 1)) {
          #  print("warning")
          #}
          if (SNF_ages[i+i-1] <= 15) {
            SNF_prop = 0
            CF_prop = 0
          }
          SNF_props = c(SNF_props, SNF_prop)
          CF_props = c(CF_props, CF_prop)
          
          age = SNF_ages[i+i-1]
          #This is where you have to insert the values
          #Replace simulation, sim1-7
          simulation = rnorm(1, mean = 0.25, sd = 0.01)
          sim2 = rnorm(1, mean = 0.12, sd = 0.01)
          sim3 = rnorm(1, mean = 0.11, sd = 0.01)
          sim4 = rnorm(1, mean = 0.11, sd = 0.01)
          sim5 = rnorm(1, mean = 0.39, sd = 0.03)
          sim6 = rnorm(1, mean = 23, sd = 3)
          sim7 = rnorm(1, mean = 0.09, sd = 0.003)
          
          string_mash = paste(county, age, sex, race_groups[j], eth)
          # stdev_hash$insert(string_mash, bmi_std)
          
          hcw_dummy_hash$insert(string_mash, hcw_prop)
          educ_dummy_hash$insert(string_mash, educ_prop)
          
          # asthma_prob_dummy_hash$insert(string_mash, ab40)
          # diabetes_prob_dummy_hash$insert(string_mash, diabetes)
          # smoker_prob_dummy_hash$insert(string_mash, smkcur)
          # #print(smkcur)
          # heart_disease_prob_dummy_hash$insert(string_mash, ab34)
          # hf_prob_dummy_hash$insert(string_mash, ab52)
          # htn_prob_dummy_hash$insert(string_mash, ab29)
          # bmi_dummy_hash$insert(string_mash, bmi_mean)
          # #stdev_hash$insert(string_mash, rnorm(1, mean = 4, sd =1))
          
          CF_dummy_hash$insert(string_mash, CF_prop)
          SNF_dummy_hash$insert(string_mash, SNF_prop)
          homeless_dummy_hash$insert(string_mash, homeless_prop)
          
          
          #print(diabetes_prob_dummy_hash$find(string_mash))
        }
      }
    }
  }
}

dummy_df = data.frame("Age_Lower" = age_low,
                      "Age_Higher" = age_high,
                      "County" = county_SNF,
                      "Sex" = sex_SNF,
                      "Race" = race_overall,
                      "Ethnicity" = eth_SNF,
                      "CF_prop" = CF_props,
                      "SNF_props" = SNF_props,
                      "Educ_props" = educ_props,
                      "HCW_props" = hcw_props,
                      "Homeless" = homeless_props)


#for loop to create fake data for comorbidity stuff 

ages = c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)


compute_special_pop <- function(county, age_bin, sex, race, ethnicity) {
  # 0 for no special population
  # 1 for healthcare workers
  # 2 for incarcerated
  # 3 for SNF
  row = filter(dummy_df, County == county, Age_Lower <= age, Age_Higher > age,
               Sex == sex, Race == race, Ethnicity == ethnicity)
  SNF = row[1,7]
  CF =  row[1,8]
  if (age < 18) {
    return (0)
  } 
  sample = rnorm(1, mean = 0.2, sd = 0.05)
  if (age > 65 & sample < 0.15) {
    return (3)
  } else if (age > 65 & sample > 0.29) {
    return (1)
  } else if (age < 65 & sample > 0.27) {
    return (1)
  } else if (age < 65 & sample < 0.09) {
    return (2)
  }
  return (1)
}

compute_special_pop2 <- function(county, age_bin, age, sex, race, ethnicity) {
  # 0 for no special population
  # 1 for healthcare workers
  # 2 for incarcerated
  # 3 for SNF
  # 4 for educational occupations
  # 5 for homeless
  # 6 for essential frontline
  # 7 for essential non-frontline
  # 8 for ALF resident
  string_mash = paste(county, age_bin, sex, race, ethnicity)
  cf = CF_dummy_hash$find(string_mash)
  snf = SNF_dummy_hash$find(string_mash)
  homeless = homeless_dummy_hash$find(string_mash)
  hcw = hcw_dummy_hash$find(string_mash)
  educ = educ_dummy_hash$find(string_mash)
  if (age < 18) {
    if (age >15) {
      #half the probability that the 16-19 group is frontline/non-frontline
      frontline = frontline_hash$find(15)
      non_frontline = non_frontline_hash$find(15)
      sample = runif(1)
      if (sample < frontline) {
        return (6)
      } else if (sample >= frontline & sample <non_frontline) {
        return (7)
      }
    }
    return (0)
  } 
  
  frontline = frontline_hash$find(age_bin)
  non_frontline = non_frontline_hash$find(age_bin)
  ALF = ALF_hash$find(age_bin)
  bin1 = hcw
  bin2 = cf + bin1
  bin3 = snf + bin2
  bin4 = educ + bin3
  bin5 = homeless + bin4
  bin6 = frontline + bin5
  bin7 = non_frontline + bin6
  bin8 = bin7 + ALF
  sample = runif(1)
  if (sample <= bin1) {
    return (1)
  } else if (sample > bin1 & sample <= bin2) {
    return (2)
  } else if (sample > bin2 & sample <= bin3) {
    return (3)
  } else if (sample > bin3 & sample <= bin4) {
    return (4)
  } else if (sample > bin4 & sample <= bin5) {
    return (5)
  } else if (sample > bin5 & sample <= bin6) {
    return (6)
  } else if (sample > bin6 & sample <= bin7) {
    return (7)
  } else if (sample > bin7 & sample <= bin8) {
    return (8)
  }
  return (0)
}

compute_comorbidity <- function(county, age, sex, race, ethnicity) {
  row = filter(dummy_data, Age_Lower<age, Age_Upper>=age,County == county,
               Ethnicity == ethnicity, Race == race, Sex == sex)
  prob = row[1,7]
  lower_bound =  row[1,8]
  upper_bound = row[1,9]
  n = row[1,10]
  std_dev = (lower_bound - prob)/(-1.960)
  #simulation = rnorm(1, mean = 0.25, sd = 0.1)
  simulation = 0.5
  if (simulation > 0.25) {
    return (0)
  }
  return (1)
}

compute_comorbidity2 <- function(county, age, sex, race, ethnicity) {
  row = filter(dummy_df, County == county, Age_Lower <= age, Age_Higher > age,
               Sex == sex, Race == race, Ethnicity == ethnicity)
  SNF = row[1,7]
  CF =  row[1,8]
  age = round((age - 2)/5)*5
  
  string_mash = paste(age, county, sex, race, ethnicity)
  #stdev = stdev_hash$find(string_mash)
  #prob = asthma_prob_dummy_hash$find(string_mash)
  #p2 = diabetes_prob_dummy_hash$find(string_mash)
  #p3 = smoker_prob_dummy_hash$find(string_mash)
  #p4 = heart_disease_prob_dummy_hash$find(string_mash)
  #p5 = htn_prob_dummy_hash$find(string_mash)
  #p6 = hf_prob_dummy_hash$find(string_mash)
  #bmi = bmi_dummy_hash$find(string_mash)
  co_morbs_list = c(0,0,0,0,0,0,0)
  stdev = 0.05
  prob = 0.13
  p2 = 0.12
  p3 = 0.11
  p4 = 0.12
  p5 = 0.01*age
  p6 = 0.01*age
  bmi = 24
  simulate = rnorm(1, mean = prob, sd = stdev)
  #s2 = rnorm(1, mean = p2, sd = stdev)
  #s3 = rnorm(1, mean = p3, sd = stdev)
  #s4 = rnorm(1, mean = p4, sd = stdev)
  #s5 = rnorm(1, mean = p5, sd = stdev)
  #s6 = rnorm(1, mean = p6, sd = stdev)
  s7 = rtruncnorm(1, a=5, b=50, mean = bmi, sd = 3)
  co_morbs_list[1] = rbinom(1,1, prob)
  co_morbs_list[2] = rbinom(1,1, p2)
  co_morbs_list[3] = rbinom(1,1, p3)
  co_morbs_list[4] = rbinom(1,1, p4)
  co_morbs_list[5] = rbinom(1,1, p5)
  co_morbs_list[6] = rbinom(1,1, p6)
  co_morbs_list[7] = s7
  #print(simulate)
  if (is.na(simulate)) {
    print(age)
    print(county)
    print(sex)
    print(race)
    print(ethnicity)
    print(prob)
    print(stdev)
  }
  #if (simulate < prob) {
  #  co_morbs_list[1] = 1
  #} 
  #if (s2 < p2) {
  #  co_morbs_list[2] = 1 
  #}
  #if (s3 < p3) {
  #  co_morbs_list[3] = 1 
  #} 
  #if (s4 < p4) {
  #  co_morbs_list[4] = 1 
  #} 
  #if (s5 < p5) {
  #  co_morbs_list[5] = 1 
  #} 
  #if (s6 < p6) {
  #  co_morbs_list[6] = 1 
  #} 
  return (co_morbs_list)
}

compute_comorbidity3 <- function(county, age_bin, sex, race, ethnicity) {
  #row = filter(dummy_df, County == county, Age_Lower <= age, Age_Higher > age,
  #             Sex == sex, Race == race, Ethnicity == ethnicity)
  #SNF = row[1,7]
  #CF =  row[1,8]

  string_mash = paste(county, age_bin, sex, race, ethnicity)
  stdev = stdev_hash$find(string_mash)
  prob = asthma_prob_dummy_hash$find(string_mash)
  p2 = diabetes_prob_dummy_hash$find(string_mash)
  p3 = smoker_prob_dummy_hash$find(string_mash)
  p4 = heart_disease_prob_dummy_hash$find(string_mash)
 
  p5 = hf_prob_dummy_hash$find(string_mash)
  p6 = htn_prob_dummy_hash$find(string_mash)
  bmi = bmi_dummy_hash$find(string_mash)
  co_morbs_list = c(0,0,0,0,0,0,0)
  #stdev = 0.05
  #prob = 0.13
  #p2 = 0.12
  #p3 = 0.11
  #p4 = 0.12
  #p5 = 0.01*age
  #p6 = 0.01*age
  #bmi = 24
  #print(string_mash)
  simulate = rnorm(1, mean = prob, sd = stdev)
  #s2 = rnorm(1, mean = p2, sd = stdev)
  #s3 = rnorm(1, mean = p3, sd = stdev)
  #s4 = rnorm(1, mean = p4, sd = stdev)
  #s5 = rnorm(1, mean = p5, sd = stdev)
  #s6 = rnorm(1, mean = p6, sd = stdev)
  s7 = rtruncnorm(1, a=5, b=50, mean = bmi, sd = stdev)
  if (is.na(s7) | s7 > 50 | s7 < 5) {
    print(county)
    print(age_bin)
    print(sex)
    print(race)
    print(ethnicity)
    print(s7)
  }
  if (is.na(rbinom(1,1, prob))) {
    print(prob)
  }
  co_morbs_list[1] = rbinom(1,1, prob)
  co_morbs_list[2] = rbinom(1,1, p2)
  co_morbs_list[3] = rbinom(1,1, p3)
  co_morbs_list[4] = rbinom(1,1, p4)
  co_morbs_list[5] = rbinom(1,1, p5)
  co_morbs_list[6] = rbinom(1,1, p6)
  co_morbs_list[7] = s7
  #print(p3)
  if (is.na(co_morbs_list[3])) {
    print(age)
    print(county)
    print(sex)
    print(race)
    print(ethnicity)
    print(p3)
    print(stdev)
  }
  #if (simulate < prob) {
  #  co_morbs_list[1] = 1
  #} 
  #if (s2 < p2) {
  #  co_morbs_list[2] = 1 
  #}
  #if (s3 < p3) {
  #  co_morbs_list[3] = 1 
  #} 
  #if (s4 < p4) {
  #  co_morbs_list[4] = 1 
  #} 
  #if (s5 < p5) {
  #  co_morbs_list[5] = 1 
  #} 
  #if (s6 < p6) {
  #  co_morbs_list[6] = 1 
  #} 
  return (co_morbs_list)
}

age_categories = c('Estimate!!Total!!Male!!Under 5 years',
                   'Estimate!!Total!!Male!!5 to 9 years',
                   'Estimate!!Total!!Male!!10 to 14 years',
                   'Estimate!!Total!!Male!!15 to 17 years',
                   'Estimate!!Total!!Male!!18 and 19 years',
                   # 'Estimate!!Total!!Male!!20 years',
                   # 'Estimate!!Total!!Male!!21 years',
                   'Estimate!!Total!!Male!!20 to 24 years',
                   'Estimate!!Total!!Male!!25 to 29 years',
                   'Estimate!!Total!!Male!!30 to 34 years',
                   'Estimate!!Total!!Male!!35 to 44 years',
                   # 'Estimate!!Total!!Male!!40 to 44 years',
                   'Estimate!!Total!!Male!!45 to 54 years',
                   # 'Estimate!!Total!!Male!!50 to 54 years',
                   'Estimate!!Total!!Male!!55 to 64 years',
                   # 'Estimate!!Total!!Male!!60 and 61 years',
                   #'Estimate!!Total!!Male!!62 to 64 years',
                   'Estimate!!Total!!Male!!65 to 74 years',
                   # 'Estimate!!Total!!Male!!67 to 69 years',
                   # 'Estimate!!Total!!Male!!70 to 74 years',
                   'Estimate!!Total!!Male!!75 to 84 years',
                   # 'Estimate!!Total!!Male!!80 to 84 years',
                   'Estimate!!Total!!Male!!85 years and over',
                   'Estimate!!Total!!Female!!Under 5 years',
                   'Estimate!!Total!!Female!!5 to 9 years',
                   'Estimate!!Total!!Female!!10 to 14 years',
                   'Estimate!!Total!!Female!!15 to 17 years',
                   'Estimate!!Total!!Female!!18 and 19 years',
                   #'Estimate!!Total!!Female!!20 years',
                   #'Estimate!!Total!!Female!!21 years',
                   'Estimate!!Total!!Female!!20 to 24 years',
                   'Estimate!!Total!!Female!!25 to 29 years',
                   'Estimate!!Total!!Female!!30 to 34 years',
                   'Estimate!!Total!!Female!!35 to 44 years',
                   #'Estimate!!Total!!Female!!40 to 44 years',
                   'Estimate!!Total!!Female!!45 to 54 years',
                   #'Estimate!!Total!!Female!!50 to 54 years',
                   'Estimate!!Total!!Female!!55 to 64 years',
                   #'Estimate!!Total!!Female!!60 and 61 years',
                   #'Estimate!!Total!!Female!!62 to 64 years',
                   'Estimate!!Total!!Female!!65 to 74 years',
                   #'Estimate!!Total!!Female!!67 to 69 years',
                   #'Estimate!!Total!!Female!!70 to 74 years',
                   'Estimate!!Total!!Female!!75 to 84 years',
                   #'Estimate!!Total!!Female!!80 to 84 years',
                   'Estimate!!Total!!Female!!85 years and over')

sample_co_morb <- function(fn, county, age_list, sex, race, ethnicity){
  output = vector()
  for (i in c(1:length(age_list))) {
    output = append(output, fn(county, age_list[i], sex, race, ethnicity[i]))
  }
  return (output)
  
}


simulate_county <-function(name) {
  
  county <- filter(census, NAME == name)
  # population estimate: B01001_001
  ID_vec <- c(1:county[3,4])
  alive_vec <- rep(1, county[3,4])
  
  #I commented out some of these age categories because while there is data availalbe
  # on the total number of ppl in that group, there isn't any race/ethnicity data
  
  
  #age_ranges = c(tuple(0,4), tuple(5,9), tuple(10,14), tuple(15,17), tuple(18,19),
  #               tuple(20,20), tuple(21,21), tuple(22,24), tuple(25,29), tuple(30,34),
  #               tuple(35,39), tuple(40,44), tuple(45,49), tuple(50,54), tuple(55,59),
  #               tuple(60,61), tuple(62,64), tuple(65,66), tuple(67,69), tuple(70,74),
  #               tuple(75,79), tuple(80,84), tuple(85,100))
  age_ranges = rep(c(c(0,4), c(5,9), c(10,14), c(15,17), c(18,19),
                     c(20,24), c(25,29), c(30,34), c(35,44), c(45,54), c(55,64),
                     c(65,74), c(75,84), c(85,100)), 2)
  
  #for the vector of each person's sex, first n entries are male, second m are female
  county_vec = rep(name, county[3,4])
  sex_vec = rep('Male', county[4,4])
  sex_vec = append(sex_vec, rep('Female', county[28,4]))
  race_vec = vector()
  age_vec = vector()
  ethnicity_vec = vector()
  co_morb_vec1 = vector()
  co_morb_vec2 = vector()
  co_morb_vec3 = vector()
  special_pop = vector()
  
  #pre compute probability of ethnicity based on race
  #if both the number of hispanic and non-hispanic individuals for a race is 0
  #then the proportion hispanic is 0
  white_prop_hispanic = 0
  african_american_prop_hispanic = 0
  AIAN_prop_hispanic = 0
  asian_prop_hispanic = 0
  native_hawaiian_prop_hispanic = 0
  native_hawaiian_prop_hispanic = 0
  other_prop_hispanic = 0
  two_plus_prop_hispanic = 0
  
  if (max(county[343,4], county[333,4]) > 0) {
    white_prop_hispanic = county[343,4]/(county[333,4] + county[343,4])
  } 
  if (max(county[344,4], county[334,4]) > 0) {
    african_american_prop_hispanic = county[344,4]/(county[334,4] + county[344,4])
  }
  if (max(county[345,4], county[335,4]) > 0) {
    AIAN_prop_hispanic = county[345,4]/(county[335,4] + county[345,4])
  }
  if (max(county[346,4], county[336,4]) > 0) {
    asian_prop_hispanic = county[346,4]/(county[336,4] + county[346,4])
  }
  if (max(county[347,4], county[337,4]) > 0) {
    native_hawaiian_prop_hispanic = county[347,4]/(county[337,4] + county[347,4])
  }
  if (max(county[348,4], county[338,4]) > 0) {
    other_prop_hispanic = county[348,4]/(county[338,4] + county[348,4])
  }
  if (max(county[349,4], county[339,4]) > 0) {
    two_plus_prop_hispanic = county[349,4]/(county[339,4] + county[349,4])
  }
  
  i = 1
  j = 2
  county2 = county %>% filter(concept != "SEX BY AGE")
  for (elem in age_categories) {
    curr_sex = "Male"
    if (i > 15) {
      curr_sex = "Female"
    }
    #first filter the table so that we only have racial categories for each age
    tr_county = filter(county2, label == elem)
    num_white = tr_county[1,4]
    num_african_american = tr_county[2,4]
    num_AIAN = tr_county[3,4]
    num_asian =  tr_county[4,4]
    num_native_hawaiian_pac_islander = tr_county[5,4]
    num_other = tr_county[6,4]
    num_two_plus =  tr_county[7,4]
  
    
    race_vec = append(race_vec, rep("White", num_white))
    race_vec = append(race_vec, rep("African American", num_african_american))
    race_vec = append(race_vec, rep("AIAN", num_AIAN))
    race_vec = append(race_vec, rep("Asian Alone", num_asian))
    race_vec = append(race_vec, rep("Native Hawaiian And Other Pacific Islander Alone", 
                                    num_native_hawaiian_pac_islander))
    race_vec = append(race_vec, rep("Some other race alone", num_other))
    race_vec = append(race_vec, rep("Two or more races", num_two_plus))
    
    age_group_total = sum(tr_county$estimate) - tr_county[8,4] - tr_county[9,4]
    age_vec = append(age_vec, sample(age_ranges[i]:age_ranges[j], 
                                     age_group_total, replace = T))
    
    prop_hisp = tr_county[9,4]/age_group_total
    estimated_hisp = (num_white*white_prop_hispanic + num_african_american*african_american_prop_hispanic
                    + num_AIAN*AIAN_prop_hispanic + num_asian*asian_prop_hispanic 
                    + num_native_hawaiian_pac_islander*native_hawaiian_prop_hispanic 
                    + num_other*other_prop_hispanic + num_two_plus*two_plus_prop_hispanic)/age_group_total
    

    age_specific_scaler = prop_hisp/(estimated_hisp)
    #print(age_specific_scaler)
    #age_specific_scaler = 1
  
    i = i + 2
    j = j + 2
    ethnicity_coin = c("Hispanic or Latino", "Not Hispanic or Latino")
    
    white_prop_rounded = round(min(1, age_specific_scaler*white_prop_hispanic) * num_white)
    aa_prop_rounded = round(min(1, age_specific_scaler *african_american_prop_hispanic) * num_african_american)
    aian_prop_rounded = round(min(1, age_specific_scaler*AIAN_prop_hispanic) * num_AIAN)
    asian_prop_rounded = round(min(1, age_specific_scaler*asian_prop_hispanic)* num_asian)
    nh_prop_rounded = round(min(1,age_specific_scaler*native_hawaiian_prop_hispanic) * num_native_hawaiian_pac_islander)
    other_prop_rounded = round(min(1,age_specific_scaler*other_prop_hispanic) * num_other)
    multi_prop_rounded = round(min(1,age_specific_scaler*two_plus_prop_hispanic) * num_two_plus)
    
    
    sampling_vec_white =  c(rep(c("Hispanic or Latino"), white_prop_rounded), rep(c("Not Hispanic or Latino"), num_white-white_prop_rounded))
    
    sampling_vec_aa =  c(rep(c("Hispanic or Latino"), aa_prop_rounded), rep(c("Not Hispanic or Latino"), num_african_american-aa_prop_rounded))
    
    sampling_vec_aian =  c(rep(c("Hispanic or Latino"), aian_prop_rounded), rep(c("Not Hispanic or Latino"), num_AIAN-aian_prop_rounded))
    
    sampling_vec_asian =  c(rep(c("Hispanic or Latino"), asian_prop_rounded), rep(c("Not Hispanic or Latino"), num_asian-asian_prop_rounded))
    
    sampling_vec_nh =  c(rep(c("Hispanic or Latino"), nh_prop_rounded), rep(c("Not Hispanic or Latino"), num_native_hawaiian_pac_islander-nh_prop_rounded))
    
    sampling_vec_other =  c(rep(c("Hispanic or Latino"), other_prop_rounded), rep(c("Not Hispanic or Latino"), num_other-other_prop_rounded))
    
    sampling_vec_multi =  c(rep(c("Hispanic or Latino"), multi_prop_rounded), rep(c("Not Hispanic or Latino"), num_two_plus-multi_prop_rounded))
    
    #ethnicity_vec = append(ethnicity_vec, sample(ethnicity_coin, size = num_white, 
    #                                             replace = TRUE, prob = c(white_prop_hispanic, 1- white_prop_hispanic)))
    #ethnicity_vec = append(ethnicity_vec, sample(ethnicity_coin, size = num_african_american, 
    #                                             replace = TRUE, prob = c(african_american_prop_hispanic, 1 - african_american_prop_hispanic)))
    #ethnicity_vec = append(ethnicity_vec, sample(ethnicity_coin, size = num_AIAN, 
    #                                             replace = TRUE, prob = c(AIAN_prop_hispanic, 1- AIAN_prop_hispanic)))
    #ethnicity_vec = append(ethnicity_vec, sample(ethnicity_coin, size = num_asian, 
    #                                             replace = TRUE, prob = c(asian_prop_hispanic, 1- asian_prop_hispanic)))
    #ethnicity_vec = append(ethnicity_vec, sample(ethnicity_coin, size = num_native_hawaiian_pac_islander, 
    #                                             replace = TRUE, prob = c(native_hawaiian_prop_hispanic, 1- native_hawaiian_prop_hispanic)))
    #ethnicity_vec = append(ethnicity_vec, sample(ethnicity_coin, size = num_other, 
    #                                             replace = TRUE, prob = c(other_prop_hispanic, 1- other_prop_hispanic)))
    #ethnicity_vec = append(ethnicity_vec, sample(ethnicity_coin, size = num_two_plus, 
    #                                             replace = TRUE, prob = c(two_plus_prop_hispanic, 1- two_plus_prop_hispanic)))
    
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_white, num_white, 
                                                 replace = FALSE))
    
    age_vec_len  = length(age_vec) - age_group_total
    age_vec_curr = age_vec_len + num_white
    eth_vec_len = length(ethnicity_vec)
    last_few = eth_vec_len - num_white
    ages = age_vec[age_vec_len:age_vec_curr]
    ethnicities = ethnicity_vec[last_few:eth_vec_len]
    #compute_comorbidity2(name, ages, curr_sex, "White", ehtnicities)
    #co_morb_vec1 = append(co_morb_vec1, sample_co_morb(compute_comorbidity, name, 
    #                                                   ages, curr_sex, "White", ethnicities))
    
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_aa, num_african_american, 
                                                 replace = FALSE))
    
    age_vec_len  = age_vec_curr
    age_vec_curr = age_vec_len + num_african_american
    eth_vec_len = length(ethnicity_vec)
    last_few = eth_vec_len - num_african_american
    ages = age_vec[age_vec_len:age_vec_curr]
    ethnicities = ethnicity_vec[last_few:eth_vec_len]
    #co_morb_vec1 = append(co_morb_vec1, sample_co_morb(compute_comorbidity, name,
     #                                                  ages, curr_sex, "African American", ethnicities))
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_aian, num_AIAN, 
                                                 replace = FALSE))
    
    age_vec_len  = age_vec_curr
    age_vec_curr = age_vec_len + num_AIAN
    eth_vec_len = length(ethnicity_vec)
    last_few = eth_vec_len - num_AIAN
    ages = age_vec[age_vec_len:age_vec_curr]
    ethnicities = ethnicity_vec[last_few:eth_vec_len]
    #co_morb_vec1 = append(co_morb_vec1, sample_co_morb(compute_comorbidity,name,
     #                                                  ages, curr_sex, "AIAN", ethnicities))
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_asian, num_asian, 
                                                 replace = FALSE))
    
    age_vec_len  = age_vec_curr
    age_vec_curr = age_vec_len + num_asian
    eth_vec_len = length(ethnicity_vec)
    last_few = eth_vec_len - num_asian
    ages = age_vec[age_vec_len:age_vec_curr]
    ethnicities = ethnicity_vec[last_few:eth_vec_len]
    #co_morb_vec1 = append(co_morb_vec1, sample_co_morb(compute_comorbidity, name,
    #                                                   ages, curr_sex, "Asian", ethnicities))
    
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_nh, num_native_hawaiian_pac_islander, 
                                                 replace = FALSE))
    
    age_vec_len  = age_vec_curr
    age_vec_curr = age_vec_len + num_native_hawaiian_pac_islander
    eth_vec_len = length(ethnicity_vec)
    last_few = eth_vec_len - num_native_hawaiian_pac_islander
    ages = age_vec[age_vec_len:age_vec_curr]
    ethnicities = ethnicity_vec[last_few:eth_vec_len]
    #co_morb_vec1 = append(co_morb_vec1, sample_co_morb(compute_comorbidity, name, ages, 
    #                curr_sex, "Native Hawaiian or Pacific Islander", ethnicities))
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_other, num_other, 
                                                 replace = FALSE))
    
    age_vec_len  = age_vec_curr
    age_vec_curr = age_vec_len + num_other
    eth_vec_len = length(ethnicity_vec)
    last_few = eth_vec_len - num_other
    ages = age_vec[age_vec_len:age_vec_curr]
    ethnicities = ethnicity_vec[last_few:eth_vec_len]
    #co_morb_vec1 = append(co_morb_vec1, sample_co_morb(compute_comorbidity, name, ages, 
    #                                                   curr_sex, "Other", ethnicities))
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_multi, num_two_plus, 
                                                 replace = FALSE))
    
    age_vec_len  = age_vec_curr
    age_vec_curr = age_vec_len + num_two_plus
    eth_vec_len = length(ethnicity_vec)
    last_few = eth_vec_len - num_two_plus
    ages = age_vec[age_vec_len:age_vec_curr]
    ethnicities = ethnicity_vec[last_few:eth_vec_len]
    #co_morb_vec1 = append(co_morb_vec1, sample_co_morb(compute_comorbidity, name, ages, 
    #                                                   curr_sex, "Multi", ethnicities))
  }
  asthma_vec = vector()
  diabetes_vec = vector()
  special_pop = vector()
  smoker_vec = vector()
  heart_disease = vector()
  heart_failure = vector()
  htn = vector()
  bmi = vector()
  status_df = 0
  stauts = vector()
  status_df <- foreach (i = 1:length(ID_vec), .combine = rbind) %dopar%{
    age = age_vec[i]
    if (age < 5) {
      age_bin = 0
    } else if (age < 10) {
      age_bin = 5
    } else if (age < 15) {
      age_bin = 10
    } else if (age < 18) {
      age_bin = 15
    } else if (age < 25) {
      age_bin = 18
    } else if (age < 35) {
      age_bin = 25
    } else if (age < 45) {
      age_bin = 35
    } else if (age < 55) {
      age_bin = 45
    } else if (age < 65) {
      age_bin = 55
    } else if (age < 75) {
      age_bin = 65
    } else if (age < 85) {
      age_bin = 75
    } else {
      age_bin = 85
    }
    sex = sex_vec[i]
    race = race_vec[i]
    county = county_vec[i]
    ethnicity = ethnicity_vec[i]
    status  = c(compute_comorbidity3(county, age_bin, sex, race, ethnicity), 
                compute_special_pop(county, age_bin, age, sex, race, ethnicity))
    #asthma_vec = c(asthma_vec, status[1])
    #diabetes_vec = c(diabetes_vec, status[2])
    #special_pop = c(special_pop, compute_special_pop(county, age, sex, race, ethnicity))
    #smoker_vec = c(smoker_vec, status[3])
    #heart_disease = c(heart_disease, status[4])
    #heart_failure = c(heart_failure, status[5])
    #htn = c(htn, status[6])
    #bmi = c(bmi, status[7])
  }
  special_pop = status_df[,8]
  asthma_vec = status_df[,1]
  diabetes_vec = status_df[,2]
  smoker_vec = status_df[,3]
  heart_disease = status_df[,4]
  heart_failure =  status_df[,5]
  htn = status_df[,6]
  bmi = status_df[,7]
  status_df = 0
  county_df = data.frame("ID" = ID_vec, 
                         "Alive/Dead_Stats" = alive_vec, 
                         "County" = county_vec,
                         "Age" = age_vec, 
                         "Sex" = sex_vec, 
                         "Race" = race_vec, 
                         "Ethnicity" = ethnicity_vec,
                         "BMI" = bmi,
                         "Special Population" = special_pop,
                         "Asthma" = asthma_vec,
                         "Diabetes" = diabetes_vec,
                         "Smoker" = smoker_vec,
                         "Heart Disease" = heart_disease,
                         "Heart Failure" = heart_failure,
                         "Hypertension" = htn)


  return (county_df)
}

simulate_county2 <-function(name) {
  
  county <- filter(census, NAME == name)
  # population estimate: B01001_001
  ID_vec <- c(1:county[3,4])
  alive_vec <- rep(1, county[3,4])
  
  #I commented out some of these age categories because while there is data availalbe
  # on the total number of ppl in that group, there isn't any race/ethnicity data
  
  
  #age_ranges = c(tuple(0,4), tuple(5,9), tuple(10,14), tuple(15,17), tuple(18,19),
  #               tuple(20,20), tuple(21,21), tuple(22,24), tuple(25,29), tuple(30,34),
  #               tuple(35,39), tuple(40,44), tuple(45,49), tuple(50,54), tuple(55,59),
  #               tuple(60,61), tuple(62,64), tuple(65,66), tuple(67,69), tuple(70,74),
  #               tuple(75,79), tuple(80,84), tuple(85,100))
  age_ranges = rep(c(c(0,4), c(5,9), c(10,14), c(15,17), c(18,19),
                     c(20,24), c(25,29), c(30,34), c(35,44), c(45,54), c(55,64),
                     c(65,74), c(75,84), c(85,100)), 2)
  
  #for the vector of each person's sex, first n entries are male, second m are female
  county_vec = rep(name, county[3,4])
  sex_vec = rep('Male', county[4,4])
  sex_vec = append(sex_vec, rep('Female', county[28,4]))
  race_vec = vector()
  age_vec = vector()
  ethnicity_vec = vector()
  co_morb_vec1 = vector()
  co_morb_vec2 = vector()
  co_morb_vec3 = vector()
  special_pop = vector()
  
  #pre compute probability of ethnicity based on race
  #if both the number of hispanic and non-hispanic individuals for a race is 0
  #then the proportion hispanic is 0
  white_prop_hispanic = 0
  african_american_prop_hispanic = 0
  AIAN_prop_hispanic = 0
  asian_prop_hispanic = 0
  native_hawaiian_prop_hispanic = 0
  native_hawaiian_prop_hispanic = 0
  other_prop_hispanic = 0
  two_plus_prop_hispanic = 0
  
  if (max(county[343,4], county[333,4]) > 0) {
    white_prop_hispanic = county[343,4]/(county[333,4] + county[343,4])
  } 
  if (max(county[344,4], county[334,4]) > 0) {
    african_american_prop_hispanic = county[344,4]/(county[334,4] + county[344,4])
  }
  if (max(county[345,4], county[335,4]) > 0) {
    AIAN_prop_hispanic = county[345,4]/(county[335,4] + county[345,4])
  }
  if (max(county[346,4], county[336,4]) > 0) {
    asian_prop_hispanic = county[346,4]/(county[336,4] + county[346,4])
  }
  if (max(county[347,4], county[337,4]) > 0) {
    native_hawaiian_prop_hispanic = county[347,4]/(county[337,4] + county[347,4])
  }
  if (max(county[348,4], county[338,4]) > 0) {
    other_prop_hispanic = county[348,4]/(county[338,4] + county[348,4])
  }
  if (max(county[349,4], county[339,4]) > 0) {
    two_plus_prop_hispanic = county[349,4]/(county[339,4] + county[349,4])
  }
  
  i = 1
  j = 2
  county2 = county %>% filter(concept != "SEX BY AGE")
  pop_so_far = 1
  county_df = 0
  for (elem in age_categories) {
    curr_sex = "Male"
    if (i > 15) {
      curr_sex = "Female"
    }
    #first filter the table so that we only have racial categories for each age
    tr_county = filter(county2, label == elem)
    num_white = tr_county[1,4]
    num_african_american = tr_county[2,4]
    num_AIAN = tr_county[3,4]
    num_asian =  tr_county[4,4]
    num_native_hawaiian_pac_islander = tr_county[5,4]
    num_other = tr_county[6,4]
    num_two_plus =  tr_county[7,4]
    
    
    race_vec = append(race_vec, rep("White", num_white))
    race_vec = append(race_vec, rep("African American", num_african_american))
    race_vec = append(race_vec, rep("AIAN", num_AIAN))
    race_vec = append(race_vec, rep("Asian Alone", num_asian))
    race_vec = append(race_vec, rep("Native Hawaiian And Other Pacific Islander Alone", 
                                    num_native_hawaiian_pac_islander))
    race_vec = append(race_vec, rep("Some other race alone", num_other))
    race_vec = append(race_vec, rep("Two or more races", num_two_plus))
    
    age_group_total = sum(tr_county$estimate) - tr_county[8,4] - tr_county[9,4]
    age_vec = append(age_vec, sample(age_ranges[i]:age_ranges[j], 
                                     age_group_total, replace = T))
    prop_hisp = tr_county[9,4]/age_group_total
    estimated_hisp = (num_white*white_prop_hispanic + num_african_american*african_american_prop_hispanic
                      + num_AIAN*AIAN_prop_hispanic + num_asian*asian_prop_hispanic 
                      + num_native_hawaiian_pac_islander*native_hawaiian_prop_hispanic 
                      + num_other*other_prop_hispanic + num_two_plus*two_plus_prop_hispanic)/age_group_total
    
    
    age_specific_scaler = prop_hisp/(estimated_hisp)
    #print(age_specific_scaler)
    #age_specific_scaler = 1
    
    i = i + 2
    j = j + 2
    ethnicity_coin = c("Hispanic or Latino", "Not Hispanic or Latino")
    
    white_prop_rounded = round(min(1, age_specific_scaler*white_prop_hispanic) * num_white)
    aa_prop_rounded = round(min(1, age_specific_scaler *african_american_prop_hispanic) * num_african_american)
    aian_prop_rounded = round(min(1, age_specific_scaler*AIAN_prop_hispanic) * num_AIAN)
    asian_prop_rounded = round(min(1, age_specific_scaler*asian_prop_hispanic)* num_asian)
    nh_prop_rounded = round(min(1,age_specific_scaler*native_hawaiian_prop_hispanic) * num_native_hawaiian_pac_islander)
    other_prop_rounded = round(min(1,age_specific_scaler*other_prop_hispanic) * num_other)
    multi_prop_rounded = round(min(1,age_specific_scaler*two_plus_prop_hispanic) * num_two_plus)
    
    
    sampling_vec_white =  c(rep(c("Hispanic or Latino"), white_prop_rounded), rep(c("Not Hispanic or Latino"), num_white-white_prop_rounded))
    
    sampling_vec_aa =  c(rep(c("Hispanic or Latino"), aa_prop_rounded), rep(c("Not Hispanic or Latino"), num_african_american-aa_prop_rounded))
    
    sampling_vec_aian =  c(rep(c("Hispanic or Latino"), aian_prop_rounded), rep(c("Not Hispanic or Latino"), num_AIAN-aian_prop_rounded))
    
    sampling_vec_asian =  c(rep(c("Hispanic or Latino"), asian_prop_rounded), rep(c("Not Hispanic or Latino"), num_asian-asian_prop_rounded))
    
    sampling_vec_nh =  c(rep(c("Hispanic or Latino"), nh_prop_rounded), rep(c("Not Hispanic or Latino"), num_native_hawaiian_pac_islander-nh_prop_rounded))
    
    sampling_vec_other =  c(rep(c("Hispanic or Latino"), other_prop_rounded), rep(c("Not Hispanic or Latino"), num_other-other_prop_rounded))
    
    sampling_vec_multi =  c(rep(c("Hispanic or Latino"), multi_prop_rounded), rep(c("Not Hispanic or Latino"), num_two_plus-multi_prop_rounded))
    
    
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_white, num_white, 
                                                 replace = FALSE))
    
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_aa, num_african_american, 
                                                 replace = FALSE))
    
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_aian, num_AIAN, 
                                                 replace = FALSE))
    
    
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_asian, num_asian, 
                                                 replace = FALSE))
    
    
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_nh, num_native_hawaiian_pac_islander, 
                                                 replace = FALSE))
    
    
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_other, num_other, 
                                                 replace = FALSE))
    
    
    ethnicity_vec = append(ethnicity_vec, sample(sampling_vec_multi, num_two_plus, 
                                                 replace = FALSE))
    asthma_vec = vector()
    diabetes_vec = vector()
    special_pop = vector()
    smoker_vec = vector()
    heart_disease = vector()
    heart_failure = vector()
    htn = vector()
    bmi = vector()
    status_df = 0
    statuts = vector()
    next_block = pop_so_far+age_group_total - 1
    if (next_block >= pop_so_far) {
      status_df <- foreach (i = pop_so_far:next_block, .combine = rbind) %dopar%{
        age = age_vec[i]
        if (age < 5) {
          age_bin = 0
        } else if (age < 10) {
          age_bin = 5
        } else if (age < 15) {
          age_bin = 10
        } else if (age < 18) {
          age_bin = 15
        } else if (age < 25) {
          age_bin = 18
        } else if (age < 35) {
          age_bin = 25
        } else if (age < 45) {
          age_bin = 35
        } else if (age < 55) {
          age_bin = 45
        } else if (age < 65) {
          age_bin = 55
        } else if (age < 75) {
          age_bin = 65
        } else if (age < 85) {
          age_bin = 75
        } else {
          age_bin = 85
        }
        sex = sex_vec[i]
        race = race_vec[i]
        county = county_vec[i]
        ethnicity = ethnicity_vec[i]
        # status  = c(compute_comorbidity3(county, age_bin, sex, race, ethnicity), 
        #             compute_special_pop2(county, age_bin, age, sex, race, ethnicity))
        status  = c(compute_comorbidity2(county, age_bin, sex, race, ethnicity), 
                    compute_special_pop2(county, age_bin, age, sex, race, ethnicity))
      }
      #print(status_df)
      if (next_block == pop_so_far) {
        special_pop = status_df[8]
        asthma_vec = status_df[1]
        diabetes_vec = status_df[2]
        smoker_vec = status_df[3]
        heart_disease = status_df[4]
        heart_failure =  status_df[5]
        htn = status_df[6]
        bmi = status_df[7]
      } else {
        special_pop = status_df[,8]
        asthma_vec = status_df[,1]
        diabetes_vec = status_df[,2]
        smoker_vec = status_df[,3]
        heart_disease = status_df[,4]
        heart_failure =  status_df[,5]
        htn = status_df[,6]
        bmi = status_df[,7]
      }
      status_df = 0
      next_block = pop_so_far+age_group_total -1
      this_df = data.frame("ID" = ID_vec[pop_so_far:next_block], 
                           "Alive/Dead_Stats" = alive_vec[pop_so_far:next_block], 
                           "County" = county_vec[pop_so_far:next_block],
                           "Age" = age_vec[pop_so_far:next_block], 
                           "Sex" = sex_vec[pop_so_far:next_block], 
                           "Race" = race_vec[pop_so_far:next_block], 
                           "Ethnicity" = ethnicity_vec[pop_so_far:next_block],
                           "BMI" = bmi,
                           "Special Population" = special_pop,
                           "Asthma" = asthma_vec,
                           "Diabetes" = diabetes_vec,
                           "Smoker" = smoker_vec,
                           "Heart Disease" = heart_disease,
                           "Heart Failure" = heart_failure,
                           "Hypertension" = htn)
      if (pop_so_far == 1) {
        county_df = this_df
      } else{
        county_df = rbind(county_df, this_df)
      }
    }
    
    pop_so_far = pop_so_far + age_group_total
  }
  
  return (county_df)
}

#all_CA_counties = unique(census$NAME)
all_CA_counties = c("Los Angeles County, California", "Santa Clara County, California",
                    "San Francisco County, California")
final_df = data.frame("ID" = c(), 
                      "Alive/Dead_Stats" = c(), 
                      "County" = c(),
                      "Age" = c(), 
                      "Sex" = c(), 
                      "Race" = c(), 
                      "Ethnicity" = c(),
                      "BMI" = c(),
                      "Special Population" = c(),
                      "Asthma" = c(),
                      "Diabetes" = c(),
                      "Smoker" = c(),
                      "Heart Disease" = c(),
                      "Heart Failure" = c(),
                      "Hypertension" = c())
dir.create("../Data/Simulation1226/",recursive = T)
for (county in counties_in_CA) {
  df = mclapply(county, simulate_county2, mc.cores = numCores)
  print(county)
  write.csv(df, paste0("../Data/Simulation1226/",county," df.csv"))
  #final_df = rbind(final_df, df)
}

for (county in counties_in_CA) {
  df = read.csv(paste0("../Data/Simulation1226/",county," df.csv"))
  final_df = rbind(final_df, df)
}
write.csv(final_df, "../Data/Simulation1226/All_CA_Counties.csv")

# start = proc.time()
# Sonoma_df = simulate_county2("Sonoma County, California")
# stop = proc.time() - start
# print(stop)
# 
# start = proc.time()
# Sonoma_df = mclapply("Sonoma County, California", simulate_county2, mc.cores = numCores)
# stop = proc.time() - start
# print(stop)
# 
# SC_county <- census %>% filter(NAME == 'Santa Clara County, California')
# SC_df = simulate_county("Santa Clara County, California")
# LA_df = simulate_county("Los Angeles County, California")
# Al_df = simulate_county("Alameda County, California")
# Merced_df = simulate_county("Merced County, California")
# Alpine_df = simulate_county("Alpine County, California")
# write.csv(Alpine_df, "/Users/poojanshukla/Downloads/Alpine_county.csv", row.names = TRUE)
# 
# SF_df = simulate_county("San Francisco County, California")
# write.csv(SF_df, "/Users/poojanshukla/Downloads/SF_county.csv", row.names = TRUE)
# 
# start = proc.time()
# Fresno_df = simulate_county("Fresno County, California")
# stop = proc.time() - start
# print(stop)
# #write.csv(Fresno_df, "/Users/poojanshukla/Downloads/Fresno_county.csv", row.names = TRUE)
# Sonoma_df = 0
# #more parameters that might be useful:
# #healthcare workers total: 	B24114_198 to B24114_257, indexing: c(15133:15192, 15699: 15758, 16265:16324)
# #healthcare workers male: B24115_198 to B24115_257
# #healthcare workers female: B24116_198 to B24116_257
# 
# #Nursing facility parameters
# # B26101_130 to B26101_156 for 3 types, indexing: c(22367:22393, 22464:22538, 23254:23280, 23392:23482)
# # B26103_005 to B26103I_005
# # B26201_130 to B26201_156 for 5 types
# # B26203_005 to B26203I_005
# 
# #Number of teachers
# # B24114_155 to B24114_161, indexing: c(15090: 15096, 15656: 15662, 16222: 16228 )
# # B24115_155 to B24115_161
# # B24116_155 to B24116_161
# 
# #total population in occupied housing units: B25008_001, indexing: c(20106)
# 
# #correctional facilities populations: B26101_100 to B26101_129
# #B26103_004 to B26103I_004
# #B26201_100 to B26201_129
# #B26203_004 to B26203I_004
# # indexing: c(22337:22366, 22463:22536, 23224:23253, 23391:23481)
# housing_var <- test$name[c(20106)]
# housing_label1 <- rep(test$label[c(20106)], 1) # 58 counties in CA 
# housing_label2 <- rep(test$concept[c(20106)],1) # 58 counties in CA 
# 
# census_housing <- get_acs(geography = "county",
#                      variables =  housing_var,
#                      state = "CA",
#                      year = 2018) # obtain data 
# census_housing <- data.frame(census_housing, label=housing_label1, concept=housing_label2) # label variables 
# 
# county_validation <- function(county_name, df, census_data){
#   county_filtered = census_data %>% filter(NAME == county_name)
#   hist(df$Age, xlab = "Age", main = paste(county_name,"Age Distribution in Simulation"), 
#        #breaks = c(0,5,10,15,18,20,25,30,35,45,55,65,75,85,100))
#        breaks = c(0,4,9,14,17,19,24,29,34,44,54,64,74,84,100))
#   
#   
#   age_mids = rep(c(c(4), c(7), c(12), c(16), c(19),
#                      c(22), c(27), c(32), c(40), c(50), c(60),
#                      c(70), c(80), c(90)), 2)
#   i = 1
#   age_data = c()
#   for (elem in age_categories) {
#     tr_county = filter(county_filtered, label == elem, concept != "SEX BY AGE")
#     age_group_total = sum(tr_county$estimate) - tr_county[8,4] - tr_county[9,4]
#     age_data = c(age_data, rep(age_mids[i], age_group_total))
#     i = i+1
#   }
#   hist(age_data, xlab = "Age", main = paste(county_name,"Age Distribution in Real Data"), 
#        #breaks = c(0,5,10,15,18,20,25,30,35,45,55,65,75,85,100))
#        breaks = c(0,4,9,14,17,19,24,29,34,44,54,64,74,84,100))
#   ### Now to make bar plots for race/ethnicity
#   ### Need to calculate number of hispanic/not hispanic for each race.
#   white_hisp_vec <- c(nrow(filter(df, Race == "White", Ethnicity == "Hispanic or Latino")),
#                       nrow(filter(df, Race == "White", Ethnicity == "Not Hispanic or Latino")))
#   aa_hisp_vec <- c(nrow(filter(df, Race == "African American", Ethnicity == "Hispanic or Latino")),
#                       nrow(filter(df, Race == "African American", Ethnicity == "Not Hispanic or Latino")))
#   
#   asian_hisp_vec <- c(nrow(filter(df, Race == "Asian Alone", Ethnicity == "Hispanic or Latino")),
#                       nrow(filter(df, Race == "Asian Alone", Ethnicity == "Not Hispanic or Latino")))
#   
#   AIAN_hisp_vec <- c(nrow(filter(df, Race == "AIAN", Ethnicity == "Hispanic or Latino")),
#                       nrow(filter(df, Race == "AIAN", Ethnicity == "Not Hispanic or Latino")))
#   
#   nhpi_hisp_vec <- c(nrow(filter(df, Race == "Native Hawaiian And Other Pacific Islander Alone", Ethnicity == "Hispanic or Latino")),
#                       nrow(filter(df, Race == "Native Hawaiian And Other Pacific Islander Alone", Ethnicity == "Not Hispanic or Latino")))
#   
#   other_hisp_vec <- c(nrow(filter(df, Race == "Some other race alone", Ethnicity == "Hispanic or Latino")),
#                       nrow(filter(df, Race == "Some other race alone", Ethnicity == "Not Hispanic or Latino")))
#   
#   multi_hisp_vec <- c(nrow(filter(df, Race == "Two or more races", Ethnicity == "Hispanic or Latino")),
#                       nrow(filter(df, Race == "Two or more races", Ethnicity == "Not Hispanic or Latino")))
#   
#   ###Now compute this from the census data
#   true_white_hisp_vec =c(county_filtered[343,4], county_filtered[333,4])
#   true_aa_hisp_vec = c(county_filtered[344,4], county_filtered[334,4])
#   true_AIAN_hisp_vec = c(county_filtered[345,4], county_filtered[335,4])
#   true_asian_hisp_vec = c(county_filtered[346,4], county_filtered[336,4])
#   true_nhpi_hisp_vec = c(county_filtered[347,4], county_filtered[337,4])
#   true_other_hisp_vec = c(county_filtered[348,4], county_filtered[338,4])
#   true_multi_hisp_vec = c(county_filtered[349,4], county_filtered[339,4])
#   
#   white_hisp_final <- c(white_hisp_vec[1], true_white_hisp_vec[1])
#   aa_hisp_final <- c(aa_hisp_vec[1], true_aa_hisp_vec[1])
#   AIAN_hisp_final <- c(AIAN_hisp_vec[1], true_AIAN_hisp_vec[1])
#   asian_hisp_final <- c(asian_hisp_vec[1], true_asian_hisp_vec[1])
#   nhpi_hisp_final <- c(nhpi_hisp_vec[1], true_nhpi_hisp_vec[1])
#   other_hisp_final <- c(other_hisp_vec[1], true_other_hisp_vec[1])
#   multi_hisp_final <- c(multi_hisp_vec[1], true_multi_hisp_vec[1])
#   
# 
#   white_not_hisp_final <- c(white_hisp_vec[2], true_white_hisp_vec[2])
#   aa_not_hisp_final <- c(aa_hisp_vec[2], true_aa_hisp_vec[2])
#   AIAN_not_hisp_final <- c(AIAN_hisp_vec[2], true_AIAN_hisp_vec[2])
#   asian_not_hisp_final <- c(asian_hisp_vec[2], true_asian_hisp_vec[2])
#   nhpi_not_hisp_final <- c(nhpi_hisp_vec[2], true_nhpi_hisp_vec[2])
#   other_not_hisp_final <- c(other_hisp_vec[2], true_other_hisp_vec[2])
#   multi_not_hisp_final <- c(multi_hisp_vec[2], true_multi_hisp_vec[2])
#   
#   hisp_df <- cbind(white_hisp_final, aa_hisp_final, AIAN_hisp_final,
#                    asian_hisp_final, nhpi_hisp_final, other_hisp_final, multi_hisp_final)
#   colnames(hisp_df) <- c("WH", "AA", "AIAN", "AS",
#                          "NHPI", "Other","Multi")
#   rownames(hisp_df) <- c("Hispanic or Latino", "Not Hispanic or Latino")
#   
#   not_hisp_df <- cbind(white_not_hisp_final, aa_not_hisp_final, AIAN_not_hisp_final,
#                    asian_not_hisp_final, nhpi_not_hisp_final, other_not_hisp_final, multi_not_hisp_final)
#   colnames(not_hisp_df) <- c("WH", "AA", "AIAN", "AS",
#                          "NHPI", "Other","Multi")
#   rownames(not_hisp_df) <- c("Hispanic or Latino", "Not Hispanic or Latino")
#   
#   
#   barplot(hisp_df,
#           main = paste(county_name, "Number Hispanic by Race, Simulation vs Real Data"),
#           xlab = "Race",
#           col = c("blue","yellow"), beside = TRUE
#   )
#   legend("top",
#          c("Simulated","Census Data"),
#          fill = c("blue","yellow")
#   )
#   barplot(not_hisp_df,
#           main = paste(county_name,"Number Not Hispanic by Race, Simulation vs Real Data"),
#           xlab = "Race",
#           col = c("blue","yellow"), beside = TRUE
#   )
#   legend("top",
#          c("Simulated","Census Data"),
#          fill = c("blue","yellow")
#   )
#   
#   
#   # First distribution
#   hist(df$Age, #breaks= c(0,5,10,15,18,20,25,30,35,45,55,65,75,85,100), 
#        breaks = c(0,4,9,14,17,19,24,29,34,44,54,64,74,84,100),
#        xlim=c(0,110), col=rgb(1,0,0,0.5), xlab="Age", 
#        ylab="Number of People", main=paste(county_name,"Age Distributions of Simulated vs Census Data") )
#   
#   # Second with add=T to plot on top
#   hist(age_data, #breaks= c(0,5,10,15,18,20,25,30,35,45,55,65,75,85,100), 
#        breaks = c(0,4,9,14,17,19,24,29,34,44,54,64,74,84,100),
#        xlim=c(0,110), col=rgb(0,0,1,0.5), add=T)
#   
#   # Add legend
#   legend("topright", legend=c("Simulated","Census Data"), col=c(rgb(1,0,0,0.5), 
#                                                         rgb(0,0,1,0.5)), pt.cex=2, pch=15 )
# }
# county_validation("Los Angeles County, California", LA_df, census)
# for (county_name in all_CA_counties) {
#   x <- simulate_county(county_name)
#   county_validation(county_name, x, census)
#   final_df = rbind(final_df, x)
# }
# 
# 
# special_populations_validation <- function(county_name, county_df) {
#   # 0 for no special population
#   # 1 for healthcare workers
#   # 2 for incarcerated
#   # 3 for SNF
#   # 4 for educational occupations
#   # 5 for homeless
#   num_hcw = nrow(filter(county_df, colnames(county_df)[10] == 1))
#   num_prisoners = nrow(filter(county_df, colnames(county_df)[10] == 2))
#   num_SNF = nrow(filter(county_df, colnames(county_df)[10] == 3))
#   num_educ = nrow(filter(county_df, colnames(county_df)[10] == 4))
#   num_homeless = nrow(filter(county_df, colnames(county_df)[10] == 5))
#   
#   occ_filtered = filter(census_occ, NAME == county_name)
#   pop_filtered = filter(census, NAME == county_name)
#   num_male_educ = occ_filtered[1,4] * frac_without_archivists
#   num_male_hcw = occ_filtered[2,4] + occ_filtered[5,4]
#   num_female_educ = occ_filtered[6,4] * frac_without_archivists
#   num_female_hcw = occ_filtered[7,4] + occ_filtered[10,4]
#   
#   print(num_hcw -  num_male_hcw -  num_female_hcw)
#   print(num_educ - num_male_educ - num_female_educ)
#   print(num_SNF)
#   print(num_homeless)
#   
# }
# special_populations_validation("Alpine County, California", df)
# 
# SC_county <- census %>% filter(NAME == 'Santa Clara County, California')
# SC_df = simulate_county("Santa Clara County, California")
# 
# 
# age_categories_abbrev = c('Male <5',
#                    'Male 5 to 9',
#                    'Male 10 to 14',
#                    'Male 15 to 17',
#                    'Male 18 and 19',
#                    'Male!!20 to 24',
#                    'Male 25 to 29',
#                    'Male 30 to 34',
#                    'Male 35 to 44',
#                    'Male 45 to 54',
#                    'Male 55 to 64',
#                    'Male 65 to 74',
#                    'Male 75 to 84',
#                    'Male 85+',
#                    'Female <5',
#                    'Female 5 to 9',
#                    'Female 10 to 14',
#                    'Female 15 to 17',
#                    'Female 18 and 19',
#                    'Female!!20 to 24',
#                    'Female 25 to 29',
#                    'Female 30 to 34',
#                    'Female 35 to 44',
#                    'Female 45 to 54',
#                    'Female 55 to 64',
#                    'Female 65 to 74',
#                    'Female 75 to 84',
#                    'Female 85+')
# #test ethnicity assumptions
# age_and_eth_validation <- function(county_name) {
#   county_filtered = filter(census, NAME == county_name)
#   county2 = filter(county_filtered, concept != "SEX BY AGE")
#   age_mids = rep(c(c(4), c(7), c(12), c(16), c(19),
#                    c(22), c(27), c(32), c(40), c(50), c(60),
#                    c(70), c(80), c(90)), 2)
#   hispanic_props = c()
#   for (age in age_categories) {
#     tr_county = filter(county2, label == age)
#     next_prop = tr_county[9,4]/(tr_county[1,4] +tr_county[2,4]+tr_county[3,4]+tr_county[4,4]+tr_county[5,4]+tr_county[6,4]+tr_county[7,4])
#     hispanic_props = c(hispanic_props, next_prop)
#     
#   }
#   barplot(hispanic_props[1:5],
#           main = paste(county_name, "Proportion Hispanic by Age, Real Data"),
#           xlab = "Age Group",
#           names.arg = age_categories_abbrev[1:5],
#           col = c("blue")
#   )
#   barplot(hispanic_props[6:10],
#           main = paste(county_name, "Proportion Hispanic by Age, Real Data"),
#           xlab = "Age Group",
#           names.arg = age_categories_abbrev[6:10],
#           col = c("blue")
#   )
#   barplot(hispanic_props[11:15],
#           main = paste(county_name, "Proportion Hispanic by Age, Real Data"),
#           xlab = "Age Group",
#           names.arg = age_categories_abbrev[11:15],
#           col = c("blue")
#   )
#   barplot(hispanic_props[16:20],
#           main = paste(county_name, "Proportion Hispanic by Age, Real Data"),
#           xlab = "Age Group",
#           names.arg = age_categories_abbrev[16:20],
#           col = c("blue")
#   )
#   barplot(hispanic_props[21:25],
#           main = paste(county_name, "Proportion Hispanic by Age, Real Data"),
#           xlab = "Age Group",
#           names.arg = age_categories_abbrev[21:25],
#           col = c("blue")
#   )
#   barplot(hispanic_props[26:30],
#           main = paste(county_name, "Proportion Hispanic by Age, Real Data"),
#           xlab = "Age Group",
#           names.arg = age_categories_abbrev[26:30],
#           col = c("blue")
#   )
#   barplot(hispanic_props,
#           main = paste(county_name, "Proportion Hispanic by Age, Real Data"),
#           xlab = "Age Group",
#           names.arg = age_categories_abbrev,
#           col = c("blue")
#   )
#   
#   return (hispanic_props)
# }
# SF_df =  simulate_county("San Francisco County, California")
# county_validation("Alameda County, California", Al_df, census)
# hispanic_props = age_and_eth_validation("Alameda County, California")
# 
# age_ranges = rep(c(c(0,4), c(5,9), c(10,14), c(15,17), c(18,19),
#                    c(20,24), c(25,29), c(30,34), c(35,44), c(45,54), c(55,64),
#                    c(65,74), c(75,84), c(85,100)), 2)
# i = 1
# j = 2
# curr_sex = "Male"
# df_prop_hispanic = c()
# for (elem in age_categories) {
#   if (age_ranges[i] == 0 & i >1) {
#     curr_sex =  "Female"
#   }
#   low = age_ranges[i]
#   high = age_ranges[j]
#   prop_hispanic = nrow(filter(Al_df, Age >= low, Age <= high, Sex == curr_sex, Ethnicity == "Hispanic or Latino"))/nrow(filter(Al_df, Age >= low, Age <= high, Sex == curr_sex))
#   df_prop_hispanic = c(df_prop_hispanic, prop_hispanic)
#   i = i+2
#   j = j+2
# }
# print(hispanic_props)
# print(df_prop_hispanic)
# eth_by_age_df = rbind(hispanic_props,df_prop_hispanic)
# colnames(eth_by_age_df) <- age_categories_abbrev
# rownames(eth_by_age_df) <- c("Census", "Simulated")
# 
# barplot(eth_by_age_df,
#         main = paste("Alameda County","Proportion Hispanic by Age, Simulation vs Real Data"),
#         xlab = "Age",
#         col = c("blue","yellow"), beside = TRUE
# )
# legend("top",
#        c("Census","Simulated"),
#        fill = c("blue","yellow"))
# 
# essential_workers_validation("Yuba County, California", Yuba_df)
# ALF_validation("Yuba County, California", Yuba_df)
# 
# essential_workers_validation("Sonoma County, California", Sonoma_df)
# ALF_validation("Sonoma County, California", Sonoma_df)
# 
# Alpine_comp_df = chis_data_validation("Alpine County, California", Alpine_df)
# Yuba_comp_df = chis_data_validation("Yuba County, California", Yuba_df)
# chis_data_validation("Merced County, California", Merced_df)
# 
# special_populations_validation("Yuba County, California", Yuba_df)
# special_populations_validation("Sonoma County, California", Sonoma_df)
# 
# 
# 
# special_populations_validation("Yuba County, California", Yuba_df)
# 
# Yuba_df = simulate_county2("Yuba County, California")
# Merced_df = simulate_county2("Merced County, California")
# Alpine_df = simulate_county2("Alpine County, California")
# Sonoma_df = simulate_county2("Sonoma County, California")
# #SF_df = simulate_county2("San Francisco County, California")
# # SF_df = read.csv("sf_df.csv")
# San_Bern_df = simulate_county2("San Bernardino County, California")
# 
# 
# # SF_homeless = filter(SF_df, Special.Population == 5)
# # SF_homeless_expected_total = 8035
# # SF_homeless_simulated = nrow(SF_homeless)
# # frac_african_american = nrow(filter(SF_homeless, Race == "African American"))/SF_homeless_simulated
# # frac_white = nrow(filter(SF_homeless, Race == "White"))/SF_homeless_simulated
# # frac_asian = nrow(filter(SF_homeless, Race == "Asian Alone"))/SF_homeless_simulated
# # frac_AIAN = nrow(filter(SF_homeless, Race == "AIAN"))/SF_homeless_simulated
# # frac_other = nrow(filter(SF_homeless, Race == "Some other race alone"))/SF_homeless_simulated
# 
# for (county in counties_in_CA) {
#   df = mclapply(county, simulate_county2, mc.cores = numCores)
#   print(county)
#   write.csv(df, paste("/Users/poojanshukla/Downloads/counties/", paste(county, "df.csv")))
#   #final_df = rbind(final_df, df)
# }
# 
# for (county in counties_in_CA) {
#   df = read.csv(paste("/Users/poojanshukla/Downloads/counties/", paste(county, "df.csv")))
#   final_df = rbind(final_df, df)
# }
# write.csv(final_df, "/Users/poojanshukla/Downloads/counties/All_CA_Counties.csv")
# #pops = c()
# #for (county in counties_in_CA) {
# #  pop = filter(census, NAME == county)[3,4]
# #  pops = c(pops, pop)
# #}
# #
# #cnty_totals = data.frame("County" = counties_in_CA, "Pop Total" = pops)
# #write.csv(cnty_totals,"county_totals.csv")
# 
# 
# ###############################################################################
# #variable names for occupations
# #282 to 
# #C24010_018 1302
# #C24010_017 1301
# #C24010_016 1300
# #C24010_054 1338
# #C24010_053 1337
# #C24010_052 1336
# 
# #some add ons
# #C24010_020 Male healthcare support occupations total 26013
# #C24010_014 education, males 26007
# #C24010_050 female education 26043
# #C24010_056 female support healtchare 26049
# 
# load_variables2018 <- load_variables(year=2018,dataset="acs1")
# 
# var_occ <- test$name[c(26007,26009,26010,26011,26013,26043,26045,26046,26047, 26049)] # Extract section of possibly relevant variables 
# var_occ_label1 <- rep(test$label[c(26007,26009,26010,26011,26013,26043,26045,26046,26047, 26049)], 58) # 58 counties in CA 
# var_occ_label2 <- rep(test$concept[c(26007,26009,26010,26011,26013,26043,26045,26046,26047, 26049)],58) # 58 counties in CA 
# 
# census_occ <- get_acs(geography = "county",
#                   variables =  var_occ,
#                   state = "CA",
#                   year = 2018, survey = "acs5") # obtain data 
# census_occ <- data.frame(census_occ, label=var_occ_label1, concept=var_occ_label2) # label variables 
# 
# 
# #B24010_043 to B24010_051 for teachers, librarians Males, 16119 to 16127      
# #B24010_056 to B24010_063 for healthcare Males, 16132 to 16138
# #B24010_065 to B24010_068 healthcare support Males 16141 to 16144
# #B24010_194 to B24010_202 teachers, librarians, female 16270 to 16278
# #B24010_207 to B24010_214 for healthcare, female 16283 to 16290
# #B24010_216 to B24010_219 for healthcare support, female 16292 to 16295
# 
# #B24010A_014 white alone, total education and librarian occupation, male 16393 
# #B24010A_017, B24010A_018 healthcare diagnosing and tech, white alone, male 16396,16397 16398
# #B24010A_020 healtchare support, white alone, male 16399
# #B24010A_050 education, white alone, female 16329
# #B24010A_052 to B24010A_054 white alone, healthcare, female 16431 to 16433
# #B24010A_056 healthcare support, white alone, female 16435
# 
# #apply offset of 73 for all subsequent groups
# 
# #same as above, but B for african americans
# #C for AIAN
# #D for asian alone
# #E for native hawaiian
# #F for other
# #G for multi
# var_occ_sex_alone = load_variables2018$name[c(16119:16127,16132:16138,16141:16144,16270:16278,16283:16290,16292:16295)]
# var_occ_sex_and_race_1 = c(16393,16396,16397,16398,16399,16329,16431:16433,16435)
# 
# var_occ_sex_and_race = c(var_occ_sex_and_race_1, var_occ_sex_and_race_1+73) #white, african american
# var_occ_sex_and_race = c(var_occ_sex_and_race, var_occ_sex_and_race_1+2*73) #AIAN
# var_occ_sex_and_race = c(var_occ_sex_and_race, var_occ_sex_and_race_1+3*73) #asian
# var_occ_sex_and_race = c(var_occ_sex_and_race, var_occ_sex_and_race_1+4*73) #native hawaiian
# var_occ_sex_and_race = c(var_occ_sex_and_race, var_occ_sex_and_race_1+5*73) #other
# var_occ_sex_and_race = c(var_occ_sex_and_race, var_occ_sex_and_race_1+6*73) #multi
# var_occ_sex_and_race = c(var_occ_sex_and_race, var_occ_sex_and_race_1+8*73) #hispanic or latino
# 
# var_occ_sex_and_race_names = load_variables2018$name[var_occ_sex_and_race]
# 
# var_occ_sex_alone_l1 <- rep(load_variables2018$label[c(16119:16127,16132:16138,16141:16144,16270:16278,16283:16290,16292:16295)], 40) # 58 counties in CA 
# var_occ_sex_alone_l2 <- rep(load_variables2018$concept[c(16119:16127,16132:16138,16141:16144,16270:16278,16283:16290,16292:16295)],40) # 58 counties in CA 
# 
# var_occ_sex_and_race_l1 <- rep(load_variables2018$label[var_occ_sex_and_race], 40) # 58 counties in CA 
# var_occ_sex_and_race_l2 <- rep(load_variables2018$concept[var_occ_sex_and_race],40) # 58 counties in CA 
# 
# 
# census_occ_sex_alone <- get_acs(geography = "state",
#                       variables =  var_occ_sex_alone,
#                       state = "CA",
#                       year = 2018, survey = "acs1") # obtain data 
# census_occ_sex_alone <- data.frame(census_occ_sex_alone, label=var_occ_sex_alone_l1, concept=var_occ_sex_alone_l2) # label variables 
# 
# census_occ_sex_and_race <- get_acs(geography = "state",
#                                 variables =  var_occ_sex_and_race_names,
#                                 state = "CA",
#                                 year = 2018, survey = "acs1") # obtain data 
# census_occ_sex_and_race <- data.frame(census_occ_sex_and_race, label=var_occ_sex_and_race_l1, concept=var_occ_sex_and_race_l2) # label variables 
# 
# 
