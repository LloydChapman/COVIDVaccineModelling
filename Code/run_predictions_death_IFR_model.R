rm(list=ls())

library(reshape2)
library(ggplot2)

source("prediction_functions.R")
source("calc_hosp_ICU_death_risk.R")

fnms <- list.files("../Data/",pattern = "death_regression_output_.*_2.RDS")

# Set prediction parameters
t_sim <- 28 # prediction time horizon in days
r <- 1.1^(1/round(365/2)) # per day multiplier for death risk = 10% increase over 6 months

# Read in names of aggregated counties in CDPH data
county_names <- read.csv("../Data/county_names.csv",stringsAsFactors = F)

### COVIDForecastHub predicted cases ###
CFH <- read.csv("../Data/COVIDForecastHubData/2020-12-21-all-forecasted-cases-model-data.csv",stringsAsFactors = F)
CFH <- CFH[CFH$model=="Ensemble" & (CFH$State=="California" & CFH$location_name!="California"),]
CFH$location_name <- gsub(" County, California","",CFH$location_name)

# Overwrite county names for aggregated counties
CFH$county_res <- county_names$cdph_county_res[match(CFH$location_name,county_names$county_res)]

# Aggregate cases by county
x0 <- aggregate(cbind(point,quantile_0.025,quantile_0.25,quantile_0.75,quantile_0.975) ~ county_res,CFH,sum) 

### CalCAT predicted deaths ###
CalCAT <- read.csv("../Data/CalCATdata/CalCAT_Custom_Data_2021-01-04.csv",stringsAsFactors = F)
# View(CalCAT)
names(CalCAT) <- CalCAT[2,]
names(CalCAT)[names(CalCAT)=="Modeler/Model"] <- "Model"
names(CalCAT)[names(CalCAT)=="Region Name"] <- "region_name"
names(CalCAT)[names(CalCAT)=="Type of Estimate"] <- "type"
CalCAT <- CalCAT[3:(nrow(CalCAT)-5),]

CalCAT$Estimate <- as.numeric(CalCAT$Estimate)
CalCAT$Date <- as.Date(CalCAT$Date)

# Subtract cumulative number of deaths at end of period from cumulative number at start to get deaths during period
idx <- (CalCAT$type=="deaths" & CalCAT$Date==min(CalCAT$Date) + t_sim)
CalCAT$Estimate[idx] <- CalCAT$Estimate[idx]-CalCAT$Estimate[CalCAT$type=="deaths" & CalCAT$Date==min(CalCAT$Date)]
CalCAT <- CalCAT[idx,]

# Overwrite county names for aggregated counties
CalCAT$county_res <- county_names$cdph_county_res[match(CalCAT$region_name,county_names$county_res)]

# Aggregate to get cumulative hospitalisations and deaths
x1 <- aggregate(Estimate ~ county_res + type, CalCAT,sum)
# Cast to wide format
x1_wide <- dcast(x1,county_res ~ type,value.var = "Estimate")
x2 <- merge(x0,x1_wide)

# Read in IFR data
IFR_long <- readRDS("../Data/IFR_by_age_ODriscoll.RDS")
IFR_ratio <- readRDS("../Data/IFR_ratio.RDS") # estimated ratio of CA IFR to O'Driscoll IFR 
IFR_long$median_perc <- IFR_long$median_perc*IFR_ratio
IFR_long$CI_95_LB <- IFR_long$CI_95_LB*IFR_ratio
IFR_long$CI_95_UB <- IFR_long$CI_95_UB*IFR_ratio

# Calculate death multiplier to match state-level forecast
deaths_ratio <- numeric(length(fnms))
dir1 <- "../Data/CalibrationDeathIFRModel1/"
dir.create(dir1,recursive = T)
fdir <- "../Figures/PredictionsVsCalCATDeathIFRModel1/"
for (i in 1:length(fnms)){
  datestr <- gsub("death_regression_output|\\.RDS","",fnms[i])
  # Work out death multiplier to match total deaths at state level
  glm_fit <- readRDS(paste0("../Data/",fnms[i]))
  x <- predict_deaths(glm_fit,t_sim,r,1,IFR_long)
  fnm <- paste0("pred",datestr,"_",t_sim,"days.csv")
  write.csv(x,paste0(dir1,fnm),row.names = F)
  deaths_ratio[i] <- compare_predictions(dir1,fnm,x2,fdir,T)$deaths_ratio
  # Correct total deaths with death multiplier and check deaths match state-level and county-level forecasts
  x <- predict_deaths(glm_fit,t_sim,r,deaths_ratio[i],IFR_long)
  fnm <- paste0("pred",datestr,"_",t_sim,"days_corrected_death_risk.csv")
  write.csv(x,paste0(dir1,fnm),row.names = F)
  compare_predictions(dir1,fnm,x2,fdir,T)
}
saveRDS(deaths_ratio,file=paste0(dir1,"deaths_multiplier.RDS"))

# # Plot 4-week-ahead predictions
# x <- read.csv(paste0(dir1,"pred",gsub("regression_output|\\.RDS","",fnms[1]),"_",t_sim,"days_corrected_death_risk.csv"),stringsAsFactors = F)
# plot_four_week_ahead_predictions(x,RR,"Figures/FourWeekAheadPredictionsDeathIFRModel/")
