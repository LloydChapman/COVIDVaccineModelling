load_synth_pop <- function(dir,county_names){
  fnms <- list.files(dir)
  ncounties <- length(fnms)
  
  # Read in files
  ldf <- vector("list",ncounties)
  for (i in 1:ncounties){
    ldf[[i]] <- read.csv(paste0(dir,fnms[i]), stringsAsFactors = F)[,-1]
  }
  # rbind into one data frame
  df <- do.call(rbind,ldf)
  
  # Recode variables
  df$County <- gsub(" County, California","",df$County)
  # Change county names to those in CDPH data
  df$County <- county_names$cdph_county_res[match(df$County,county_names$county_res)]
  names(df)[names(df)=="County"] <- "county_res"
  
  # Assign age groups
  age_cat_lbs <- c(seq(0,80,by=10),max(df$Age)+1)
  df$age_cat <- cut(df$Age,age_cat_lbs,right = F,labels = c("<10","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"))
  
  # Recode variables
  df$Sex <- ifelse(df$Sex=="Male",1,0)
  # df$Ethnicity <- ifelse(df$Ethnicity=="Hispanic or Latino",1,0)
  df$Obesity <- ifelse(df$BMI>30,1,0)
  
  # Recode race/ethnicity categories
  df$race_ethnicity <- ""
  df$race_ethnicity[df$Race=="African American"] <- "non-Hispanic Black"
  df$race_ethnicity[df$Race=="White"] <- "non-Hispanic White"
  df$race_ethnicity[df$Race %in% c("AIAN","Asian Alone","Native Hawaiian And Other Pacific Islander Alone","Some other race alone","Two or more races")] <- "Other"
  # Overwrite race/ethnicity for people with Ethnicity="Hispanic or Latino"
  # df$race_ethnicity[df$Ethnicity==1] <- "Hispanic/Latino"
  df$race_ethnicity[df$Ethnicity=="Hispanic or Latino"] <- "Hispanic/Latino"
  
  names(df) <- tolower(names(df))
  
  return(df)
}

add_essential_workers <- function(df){
  age_16to49 <- (df$special.population==0 & df$age>=16 & df$age<50)
  pop_EW <- 87/328*nrow(df)-sum(df$special.population %in% c(1,4))
  prop_16to49_essential <- 2/3*pop_EW/sum(age_16to49) # number of essential workers in CA who are 16-49 from CEPR frontline workers analysis
  df$special.population[age_16to49 & (runif(nrow(df)) < prop_16to49_essential)] <- 6
  age_50plus <- (df$special.population==0 & df$age>=50)
  prop_50plus_essential <- 1/3*pop_EW/sum(age_50plus) # number of essential workers in CA who are 50+ from CEPR frontline workers analysis
  df$special.population[age_50plus & (runif(nrow(df)) < prop_50plus_essential)] <- 6
  
  return(df)
}

calc_analysis_vars <- function(agg_cases,start_date,agg_pop){
  # Subset data
  z <- agg_cases[agg_cases$first_report_date >= start_date,]
  
  # Create data frame of all combinations of dates and demographic variables
  x <- expand.grid(first_report_date=seq.Date(start_date,max(z$first_report_date),by=1),county_res=unique(z$county_res),age_cat=unique(z$age_cat),sex=c(0,1),race_ethnicity=unique(z$race_ethnicity),stringsAsFactors = F)
  
  # Merge with aggregated cases data frame
  y <- merge(z,x,all=T)
  # Fill in dates without recorded cases and deaths with zeros
  y$n[is.na(y$n)] <- 0
  y$n_deaths[is.na(y$n_deaths)] <- 0
  
  # Create exposure time variable
  y$time <- as.numeric(y$first_report_date - start_date)
  
  # Calculate cumulative cases by date for each risk factor group
  y$cum_cases <- ave(y$n,list(y$county_res,y$age_cat,y$sex,y$race_ethnicity),FUN=cumsum)
  y$cum_deaths <- ave(y$n_deaths,list(y$county_res,y$age_cat,y$sex,y$race_ethnicity),FUN=cumsum)
  # # Check
  # y[order(y$race_ethnicity,y$sex,y$age_cat,y$county_res,y$first_report_date),][1:40,]
  
  ### Merge with population data ###
  # # Read in population data
  # agg_pop <- read.csv("/mnt/nlo_shared/data/agg_pop1.csv",stringsAsFactors=F)
  # agg_pop <- read.csv("../Data/agg_pop1.csv",stringsAsFactors=F)
  
  y <- merge(y,agg_pop,all.x=T)
  y <- y[,c("first_report_date",names(y)[names(y)!="first_report_date"])]
  y <- y[do.call(order,y),]
  
  # Exclude counts for individuals with unknown county or of unknown race-ethnicity
  y <- y[!(y$county_res=="UNASSIGNED"|y$race_ethnicity=="Unknown"),]
  
  return(y)
}

calc_death_analysis_vars <- function(agg_deaths,start_date,agg_pop){
  # Subset data
  z <- agg_deaths[agg_deaths$date_of_death >= start_date,]
  
  # Create data frame of all combinations of dates and demographic variables
  x <- expand.grid(date_of_death=seq.Date(start_date,max(z$date_of_death),by=1),county_res=unique(z$county_res),age_cat=unique(z$age_cat),sex=c(0,1),race_ethnicity=unique(z$race_ethnicity),stringsAsFactors = F)
  
  # Merge with aggregated cases data frame
  y <- merge(z,x,all=T)
  # Fill in dates without recorded cases and deaths with zeros
  y$n_deaths[is.na(y$n_deaths)] <- 0
  
  # Create survival time variable
  y$time_death <- as.numeric(y$date_of_death - start_date)
  
  # Calculate cumulative cases by date for each risk factor group
  y$cum_deaths <- ave(y$n_deaths,list(y$county_res,y$age_cat,y$sex,y$race_ethnicity),FUN=cumsum)
  # # Check
  # y[order(y$race_ethnicity,y$sex,y$age_cat,y$county_res,y$date_of_death),][1:40,]
  
  ### Merge with population data ###
  # # Read in population data
  # agg_pop <- read.csv("/mnt/nlo_shared/data/agg_pop1.csv",stringsAsFactors=F)
  # agg_pop <- read.csv("../Data/agg_pop1.csv",stringsAsFactors=F)
  
  y <- merge(y,agg_pop,all.x=T)
  y <- y[,c("date_of_death",names(y)[names(y)!="date_of_death"])]
  y <- y[do.call(order,y),]
  
  # Exclude counts for individuals with unknown county or of unknown race-ethnicity
  y <- y[!(y$county_res=="UNASSIGNED"|y$race_ethnicity=="Unknown"),]
  
  return(y)
}

