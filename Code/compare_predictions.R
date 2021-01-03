rm(list=ls())

library(reshape2)
library(ggplot2)

source("~/Dropbox/COVIDVaccineModelling/Code/prediction_functions.R")

dir <- "~/Dropbox/COVIDVaccineModelling/Data/"
setwd(dir)

county_names <- read.csv("county_names.csv",stringsAsFactors = F)

### COVIDForecastHub predicted cases ###
CFH <- read.csv("COVIDForecastHubData/2020-09-21-all-forecasted-cases-model-data.csv",stringsAsFactors = F)
CFH <- CFH[CFH$model=="Ensemble" & (CFH$State=="California" & CFH$location_name!="California"),]
CFH$location_name <- gsub(" County, California","",CFH$location_name)

# Overwrite county names for aggregated counties
CFH$county_res <- county_names$cdph_county_res[match(CFH$location_name,county_names$county_res)]

# Aggregate cases by county
x0 <- aggregate(cbind(point,quantile_0.025,quantile_0.25,quantile_0.75,quantile_0.975) ~ county_res,CFH,sum) 
  
### CalCAT predicted deaths ###
CalCAT <- read.csv("CalCATdata/CalCAT_Custom_Data_2020-11-02.csv",stringsAsFactors = F)
View(CalCAT)
names(CalCAT) <- CalCAT[2,]
names(CalCAT)[names(CalCAT)=="Modeler/Model"] <- "Model"
names(CalCAT)[names(CalCAT)=="Region Name"] <- "region_name"
names(CalCAT)[names(CalCAT)=="Type of Estimate"] <- "type"
CalCAT <- CalCAT[3:(nrow(CalCAT)-5),]

CalCAT$Estimate <- as.numeric(CalCAT$Estimate)
CalCAT$Date <- as.Date(CalCAT$Date)

# Subtract cumulative number of deaths at end of period from cumulative number at start to get deaths during period
idx <- (CalCAT$type=="deaths" & CalCAT$Date==max(CalCAT$Date))
CalCAT$Estimate[idx] <- CalCAT$Estimate[idx]-CalCAT$Estimate[CalCAT$type=="deaths" & CalCAT$Date==min(CalCAT$Date)]
CalCAT <- CalCAT[CalCAT$type=="hosp"|idx,]

# Overwrite county names for aggregated counties
CalCAT$county_res <- county_names$cdph_county_res[match(CalCAT$region_name,county_names$county_res)]

# Aggregate to get cumulative hospitalisations and deaths
x1 <- aggregate(Estimate ~ county_res + type, CalCAT,sum)
# Cast to wide format
x1_wide <- dcast(x1,county_res ~ type)
x2 <- merge(x0,x1_wide)

# fnms <- list.files(dir,pattern = "pred_.*31days.csv")
# fnms <- list.files(dir,pattern = "pred_.*28days.csv")
# fnms <- list.files(dir,pattern = "pred_.*28days_corrected.csv")
fnms <- list.files(dir,pattern = "pred_.*28days_corrected_deaths.csv")

cum_cases_ratio <- numeric(length(fnms))
cum_deaths_ratio <- numeric(length(fnms))
for (i in 1:length(fnms)){
  x <- read.csv(paste0(dir,fnms[i]),stringsAsFactors = F)
  agg_x <- aggregate(cbind(cum_cases_pred,cum_hosps_pred,cum_deaths_pred) ~ county_res,x,sum)
  agg_x <- merge(agg_x,x2)
  # names(agg_x)[names(agg_x)=="cum_deaths_pred"] <- "risk model"
  # names(agg_x)[names(agg_x)=="deaths"] <- "CalCAT"
  cum_cases <- melt(agg_x[,c("county_res","cum_cases_pred","point")],id.vars = "county_res",value.name = "cases")
  pdf(paste0("../Figures/compare_cases_",sub(".csv","",fnms[i]),"_county.pdf"),width = 7,height = 4)
  print(ggplot(cum_cases,aes(fill=variable,y=cases,x=county_res)) + geom_bar(position = "dodge", stat="identity") + xlab("County") + ylab("Predicted cases in next 4 weeks") + scale_fill_discrete(name = "",labels = c("Model prediction","Forecast Hub")) + theme(axis.text.x = element_text(angle = 45,hjust = 1)))
  dev.off()
  # cum_cases_state <- data.frame(model=names(agg_x[,c("cum_cases_pred","point")]),cases=apply(agg_x[,c("cum_cases_pred","point")],2,sum))
  cum_cases_state <- data.frame(model=factor(c("Model prediction","Forecast Hub"),level=c("Model prediction","Forecast Hub")),cases=apply(agg_x[,c("cum_cases_pred","point")],2,sum))
  cum_cases_ratio[i] <- cum_cases_state$cases[2]/cum_cases_state$cases[1]
  pdf(paste0("../Figures/compare_cases_",sub(".csv","",fnms[i]),"_state.pdf"),width = 4,height = 4)
  print(ggplot(cum_cases_state,aes(x=model,y=cases)) + geom_bar(stat = "identity") + theme(axis.title.x = element_blank()))
  dev.off()
  cum_deaths <- melt(agg_x[,c("county_res","cum_deaths_pred","deaths")],id.vars = "county_res",value.name = "deaths")
  pdf(paste0("../Figures/compare_deaths_",sub(".csv","",fnms[i]),"_county.pdf"),width = 7,height = 4)
  print(ggplot(cum_deaths,aes(fill=variable,y=deaths,x=county_res)) + geom_bar(position = "dodge", stat="identity") + xlab("County") + ylab("Predicted deaths in next 4 weeks") + scale_fill_discrete(name = "",labels = c("Model prediction","CalCAT")) + theme(axis.text.x = element_text(angle = 45,hjust = 1)))
  dev.off()
  cum_deaths_state <- data.frame(model=names(agg_x[,c("cum_deaths_pred","deaths")]),deaths=apply(agg_x[,c("cum_deaths_pred","deaths")],2,sum))
  cum_deaths_state <- data.frame(model=factor(c("Model prediction","CalCAT"),level=c("Model prediction","CalCAT")),deaths=apply(agg_x[,c("cum_deaths_pred","deaths")],2,sum))
  cum_deaths_ratio[i] <- cum_deaths_state$deaths[2]/cum_deaths_state$deaths[1]
  pdf(paste0("../Figures/compare_deaths_",sub(".csv","",fnms[i]),"_state.pdf"),width = 4,height = 4)
  print(ggplot(cum_deaths_state,aes(x=model,y=deaths)) + geom_bar(stat = "identity") + theme(axis.title.x = element_blank()))
  dev.off()
}
# saveRDS(cum_cases_ratio,file="cases_multiplier.RDS")
# saveRDS(cum_deaths_ratio,file="deaths_multiplier.RDS")

# Relative risks for HCWs, inmates, SNF residents, teachers and homeless individuals
RR <- c(HCW = 3,inmate = 5.5,SNF = 9,teacher = 1.8,homeless = 1.65)

# Plot 4-week-ahead predictions
plot_four_week_ahead_predictions(x,RR,"../Figures/FourWeekAheadPredictions/")

