rm(list=ls())

county_data <- read.csv("../county_data.csv",stringsAsFactors = F)

la_times_county_data <- read.csv("../california-coronavirus-data/latimes-county-totals.csv")
la_times_county_data$date <- as.Date(la_times_county_data$date)
la_times_county_data <- la_times_county_data[la_times_county_data$date==as.Date("2020-10-19"),]

x <- merge(county_data,la_times_county_data[,c("county","confirmed_cases")],by.x="county_res",by.y="county",all.x = T)

y <- aggregate(x$confirmed_cases,by=list(cdph_county_res=x$cdph_county_res),sum)
x <- merge(x,y,all.x=T)
