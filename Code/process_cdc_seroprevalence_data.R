rm(list=ls())

# Read in CDC seroprevalence survey data
seroprev <- read.csv("../Data/CDCSeroprevalenceData/Nationwide_Commercial_Laboratory_Seroprevalence_Survey_2021-01-04.csv",stringsAsFactors = F)

# Subset to California
seroprev <- seroprev[seroprev$Site=="CA",] #grep("Rate",names(seroprev))

# Change variable names and dates
names(seroprev)[names(seroprev)=="Date.Range.of.Specimen.Collection"] <- "Date"
seroprev[,c("Start_date","End_date")] <- do.call(rbind,strsplit(seroprev$Date," - "))
seroprev$Start_date <- trimws(seroprev$Start_date)
seroprev$Start_date <- as.Date(paste0(seroprev$Start_date,", 2020"),format = "%b %d, %Y")
seroprev$End_date <- trimws(seroprev$End_date)
seroprev$End_date <- as.Date(paste0(seroprev$End_date,", 2020"),format = "%b %d, %Y")
# names(seroprev)[names(seroprev)=="Catchment.population"] <- "Population"

# Melt to long format
names(seroprev) <- sub("([0-9])\\.([0-9])","\\1-\\2",sub("\\.$","+",gsub("\\.{2}","",sub("\\.$","",gsub("Years\\.|Prevalence\\.","",names(seroprev))))))
group <- c("0-17","18-49","50-64","65+","Male","Female","Cumulative")
seroprev_long <- reshape(seroprev,varying = list(paste0("n",group),paste0("Rate",group),paste0("Lower.CI",group),paste0("Upper.CI",group)),v.names = c("n","Rate","Lower.CI","Upper.CI"),idvar = c("Site","Start_date","End_date","Round"),direction = "long",drop = c("Date",names(seroprev)[grep("Catchment|CI.Flag|Estimated|All",names(seroprev))]))
seroprev_long$group <- group[seroprev_long$time]
seroprev_long$time <- NULL
names(seroprev_long)[names(seroprev_long)=="Rate"] <- "seroprev"
row.names(seroprev_long) <- 1:nrow(seroprev_long)
# seroprev_long <- melt(seroprev,id.vars = c("Site","Start_date","End_date","Round","Population"),measure.vars = grep("n\\.\\.|Rate",names(seroprev)),variable.name = "group",value.name = "seroprev")
# seroprev_long$group <- sub("-$","+",sub("\\.","-",sub("\\.$","",gsub("Rate[.]+|Years\\.|Prevalence\\.","",seroprev_long$group))))

# Change seroprevalence to proportion
seroprev_long$seroprev <- seroprev_long$seroprev/100
seroprev_long$Lower.CI <- seroprev_long$Lower.CI/100
seroprev_long$Upper.CI <- seroprev_long$Upper.CI/100

# Reorder
seroprev_long <- seroprev_long[order(seroprev_long$Round,seroprev_long$group),]

# Save
write.csv(seroprev_long,"../Data/CA_seroprevalence1.csv",row.names = F)
