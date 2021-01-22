rm(list=ls())

library(ggplot2)

## CDPH data
dir <- "../Data/"
fnms <- list.files(dir,pattern = "^regression_output_.*_2.RDS")

for (i in 1:length(fnms)){
  glm_fit <- readRDS(paste0(dir,fnms[i]))
  agg_x <- aggregate(cum_deaths ~ age_cat,glm_fit$data,sum)
  pdf(paste0("../Figures/DeathAgeDistns/","death_age_distn",gsub("regression_output|\\.RDS","",fnms[i]),".pdf"),width = 5, height = 4)
  # ggplot(agg_x,aes(x=age_cat,y=cum_deaths)) + geom_bar(stat="identity")
  print(ggplot(agg_x,aes(x=age_cat,y=cum_deaths/sum(cum_deaths))) + geom_bar(stat="identity") + xlab("Age") + ylab("Proportion of deaths"))
  dev.off()
}

## CDC US data
# Read in data
deaths <- read.csv("../Data/CDCDeathsData/Provisional_COVID-19_Death_Counts_by_Sex__Age__and_State_2021-01-21.csv",stringsAsFactors = F)

# Subset to national data
# deaths <- deaths[deaths$State=="United States",]
deaths$Age.group <- sub(" year.*","",deaths$Age.group)
deaths$Age.group[deaths$Age.group=="85"] <- "85+"

plot_deaths <- function(deaths,fnm){
  pdf(paste0(paste0("../Figures/DeathAgeDistns/",fnm,".pdf")),width = 5, height = 4)
  print(ggplot(deaths,aes(x=Age.group,y=COVID.19.Deaths/sum(COVID.19.Deaths,na.rm=T))) + geom_bar(stat="identity") + xlab("Age") + ylab("Proportion of deaths"))
  dev.off()
  ggplot
}

# Plot US death age distn
USdeaths <- deaths[deaths$State=="United States" & deaths$Sex=="All Sexes" & !(deaths$Age.group %in% c("All Ages","Under 1","0-17","18-29","30-49","50-64")),]
USdeaths <- USdeaths[order(as.numeric(gsub("-.*|\\+","",USdeaths$Age.group))),]
USdeaths$Age.group <- factor(USdeaths$Age.group,levels = USdeaths$Age.group)
plot_deaths(USdeaths,"US_death_age_distn2")
USdeaths1 <- deaths[deaths$Sex=="All Sexes" & (deaths$Age.group %in% c("0-17","18-29","30-49","50-64","65-74","75-84","85+")),]
plot_deaths(USdeaths1,"US_death_age_distn3")

# Plot CA death age distn
CAdeaths <- deaths[deaths$State=="California" & !(deaths$Age.group %in% c("All Ages","Under 1","0-17","18-29","30-49","50-64")),]
CAdeaths$COVID.19.Deaths[is.na(CAdeaths$COVID.19.Deaths)] <- sample.int(9,sum(is.na(CAdeaths$COVID.19.Deaths)),replace = T)
CAdeaths <- aggregate(COVID.19.Deaths ~ Age.group,CAdeaths,sum)
CAdeaths <- CAdeaths[order(as.numeric(sub("-.*|\\+","",CAdeaths$Age.group))),]
CAdeaths$Age.group <- factor(CAdeaths$Age.group,levels = CAdeaths$Age.group)
plot_deaths(CAdeaths,"CA_death_age_distn2")

CAdeaths1 <- deaths[deaths$State=="California" & (deaths$Age.group %in% c("0-17","18-29","30-49","50-64","65-74","75-84","85+")),]
CAdeaths1$COVID.19.Deaths[is.na(CAdeaths1$COVID.19.Deaths)] <- sample.int(9,sum(is.na(CAdeaths1$COVID.19.Deaths)),replace = T)
CAdeaths1 <- aggregate(COVID.19.Deaths ~ Age.group,CAdeaths1,sum)
plot_deaths(CAdeaths1,"CA_death_age_distn3")
