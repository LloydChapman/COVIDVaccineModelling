rm(list=ls())

source("~/Dropbox/COVIDVaccineModelling/Code/seroprevalence_functions.R")

setwd("~/Dropbox/COVIDVaccineModelling/Data")

# Read in seroprevalence data
seroprev_long <- read.csv("CA_seroprevalence.csv",stringsAsFactors = F)
seroprev_long$k <- round(seroprev_long$seroprev * seroprev_long$n)
seroprev <- aggregate(cbind(n,k) ~ group,seroprev_long,sum)
seroprev <- seroprev[!(seroprev$group %in% c("Cumulative","Female","Male")),]
names(seroprev)[names(seroprev)=="group"] <- "age"

seroprev$age[seroprev$age=="65+"] <- "65-100"
seroprev$age_low <- as.numeric(sub("-.*","",seroprev$age))
seroprev$age_upp <- as.numeric(sub(".*-","",seroprev$age))

ggplot(seroprev,aes(x=age,y=k/n)) + geom_bar(stat="identity") + ylab("seroprevalence")

# Minimise negative log-likelihood with optim # [ ] - Beeds updating as should be constrained optimization, which I think requires specification of the gradient of the likelihood
res <- optim(c(0.06,4e-4),seroprev_negative_log_likelihood,x=seroprev) #method="L-BFGS-B",lower=c(0,0),upper=c(1,100),
saveRDS(res,"seroprev_pars.RDS")
