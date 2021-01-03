


SNFdata <- read.delim("../SNFdata/COUNTY_DETAILS_data.csv",sep = "\t",fileEncoding = "UTF-16LE")
names(SNFdata)[names(SNFdata)=="AS_OF_DATE"] <- "Date"

