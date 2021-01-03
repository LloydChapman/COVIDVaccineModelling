rm(list=ls())

source("~/Dropbox/COVIDVaccineModelling/Code/prediction_functions.R")
source("~/Dropbox/COVIDVaccineModelling/Code/calc_hosp_ICU_death_risk.R")

dir <- "~/Dropbox/COVIDVaccineModelling/Data/"
setwd(dir)
fnms <- list.files(dir,pattern = "^regression_output_.*_1.RDS")

# Relative risks for HCWs, inmates, SNF residents, teachers and homeless individuals
RR <- c(HCW = 3,inmate = 5.5,SNF = 9,teacher = 1.8,homeless = 1.65)

# Set prediction parameters
t_sim <- 28 # prediction time horizon in days
r <- 1.1^(1/t_sim) # per day multiplier for force of infection = 10% increase per month, currently

# Get probabilities of death by age and sex
res <- calc_hosp_ICU_death_risk()
p_death <- res$p_death_given_clin_by_age_and_sex

# Read in names of aggregated counties in CDPH data
county_names <- read.csv("county_names.csv",stringsAsFactors = F)

### COVIDForecastHub predicted cases ###
# CFH <- read.csv("COVIDForecastHubData/2020-09-21-all-forecasted-cases-model-data.csv",stringsAsFactors = F)
# CFH <- read.csv("COVIDForecastHubData/2020-10-19-all-forecasted-cases-model-data.csv",stringsAsFactors = F)
CFH <- read.csv("COVIDForecastHubData/2020-11-23-all-forecasted-cases-model-data.csv",stringsAsFactors = F)
CFH <- CFH[CFH$model=="Ensemble" & (CFH$State=="California" & CFH$location_name!="California"),]
CFH$location_name <- gsub(" County, California","",CFH$location_name)

# Overwrite county names for aggregated counties
CFH$county_res <- county_names$cdph_county_res[match(CFH$location_name,county_names$county_res)]

# Aggregate cases by county
x0 <- aggregate(cbind(point,quantile_0.025,quantile_0.25,quantile_0.75,quantile_0.975) ~ county_res,CFH,sum) 

### CalCAT predicted deaths ###
CalCAT <- read.csv("CalCATdata/CalCAT_Custom_Data_2020-12-11.csv",stringsAsFactors = F)
View(CalCAT)
names(CalCAT) <- CalCAT[2,]
names(CalCAT)[names(CalCAT)=="Modeler/Model"] <- "Model"
names(CalCAT)[names(CalCAT)=="Region Name"] <- "region_name"
names(CalCAT)[names(CalCAT)=="Type of Estimate"] <- "type"
CalCAT <- CalCAT[3:(nrow(CalCAT)-5),]

CalCAT$Estimate <- as.numeric(CalCAT$Estimate)
CalCAT$Date <- as.Date(CalCAT$Date)

# Subtract cumulative number of deaths at end of period from cumulative number at start to get deaths during period
# idx <- (CalCAT$type=="deaths" & CalCAT$Date==max(CalCAT$Date))
idx <- (CalCAT$type=="deaths" & CalCAT$Date==min(CalCAT$Date) + 28)
CalCAT$Estimate[idx] <- CalCAT$Estimate[idx]-CalCAT$Estimate[CalCAT$type=="deaths" & CalCAT$Date==min(CalCAT$Date)]
CalCAT <- CalCAT[CalCAT$type=="hosp"|idx,]

# Overwrite county names for aggregated counties
CalCAT$county_res <- county_names$cdph_county_res[match(CalCAT$region_name,county_names$county_res)]

# Aggregate to get cumulative hospitalisations and deaths
x1 <- aggregate(Estimate ~ county_res + type, CalCAT,sum)
# Cast to wide format
x1_wide <- dcast(x1,county_res ~ type,value.var = "Estimate")
x2 <- merge(x0,x1_wide)

# # Calculate case and death multipliers to match state-level forecast
# cases_ratio <- numeric(length(fnms))
# deaths_ratio <- numeric(length(fnms))
# fdir <- "../Figures/PredictionsVsForecastHubDecisionTree/"
# for (i in 1:length(fnms)){
#   datestr <- gsub("regression_output|\\.RData","",fnms[i])
#   # Work out case multiplier to match total cases at state level
#   x <- predict_cases_hosps_deaths(fnms[i],t_sim,r,1,1,p_death)
#   fnm <- paste0("pred",datestr,"_",t_sim,"days.csv")
#   write.csv(x,paste0(fnm),row.names = F)
#   cases_ratio[i] <- compare_predictions(fnm,x2,fdir)$cases_ratio
#   # Correct total cases with case multiplier to work out death multiplier to match total deaths
#   x <- predict_cases_hosps_deaths(fnms[i],t_sim,r,cases_ratio[i],1,p_death)
#   fnm <- paste0("pred",datestr,"_",t_sim,"days_corrected.csv")
#   write.csv(x,paste0(fnm),row.names = F)
#   deaths_ratio[i] <- compare_predictions(fnm,x2,fdir)$deaths_ratio
#   # Check total cases and deaths at state and county level
#   x <- predict_cases_hosps_deaths(fnms[i],t_sim,r,cases_ratio[i],deaths_ratio[i],p_death)
#   fnm <- paste0("pred",datestr,"_",t_sim,"days_corrected_deaths.csv")
#   write.csv(x,paste0(fnm),row.names = F)
#   compare_predictions(fnm,x2,fdir)
# }
# saveRDS(cases_ratio,file="cases_multiplier.RDS")
# saveRDS(deaths_ratio,file="deaths_multiplier.RDS")
# 
# # Plot 4-week-ahead predictions
# x <- read.csv(paste0("pred",gsub("regression_output|\\.RData","",fnms[3]),"_",t_sim,"days_corrected_deaths.csv"),stringsAsFactors = F)
# plot_four_week_ahead_predictions(x,RR,"../Figures/FourWeekAheadPredictionsDecisionTree/")

# Calculate case multiplier to match state-level forecast
cases_ratio <- numeric(length(fnms))
deaths_ratio <- numeric(length(fnms))
dir1 <- "Calibration2/"
dir.create(dir1,recursive = T)
fdir <- "../Figures/PredictionsVsForecastHub2/"
for (i in 1:length(fnms)){
  datestr <- gsub("regression_output|\\.RDS","",fnms[i])
  p_death <- readRDS(paste0("prob_death_given_clin_by_age_and_sex",gsub("regression_output","",fnms[i])))
  # Work out case multiplier to match total cases at state level
  x <- predict_cases_hosps_deaths(fnms[i],t_sim,r,1,1,p_death)
  fnm <- paste0("pred",datestr,"_",t_sim,"days.csv")
  write.csv(x,paste0(dir1,fnm),row.names = F)
  cases_ratio[i] <- compare_predictions(dir1,fnm,x2,fdir)$cases_ratio
  # Correct total cases with case multiplier
  x <- predict_cases_hosps_deaths(fnms[i],t_sim,r,cases_ratio[i],1,p_death)
  fnm <- paste0("pred",datestr,"_",t_sim,"days_estd_death_risk.csv")
  write.csv(x,paste0(dir1,fnm),row.names = F)
  deaths_ratio[i] <- compare_predictions(dir1,fnm,x2,fdir)$deaths_ratio
  # Correct total deaths with death multiplier and check deaths match state-level and county-level forecasts
  x <- predict_cases_hosps_deaths(fnms[i],t_sim,r,cases_ratio[i],deaths_ratio[i],p_death)
  fnm <- paste0("pred",datestr,"_",t_sim,"days_corrected_death_risk.csv")
  write.csv(x,paste0(dir1,fnm),row.names = F)
  compare_predictions(dir1,fnm,x2,fdir)
}
saveRDS(cases_ratio,file=paste0(dir1,"cases_multiplier1.RDS"))
saveRDS(deaths_ratio,file=paste0(dir1,"deaths_multiplier1.RDS"))

# Plot 4-week-ahead predictions
# x <- read.csv(paste0(dir1,"pred",gsub("regression_output|\\.RDS","",fnms[3]),"_",t_sim,"days_estd_death_risk.csv"),stringsAsFactors = F)
x <- read.csv(paste0(dir1,"pred",gsub("regression_output|\\.RDS","",fnms[3]),"_",t_sim,"days_corrected_death_risk.csv"),stringsAsFactors = F)
plot_four_week_ahead_predictions(x,RR,"../Figures/FourWeekAheadPredictions2/")
