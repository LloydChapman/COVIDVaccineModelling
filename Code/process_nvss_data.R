library(readxl)

dir <- "~/Dropbox/COVIDVaccineModelling/Data/"

setwd(dir)
fnms <- list.files("USLifeTables/")

llt <- vector("list",length(fnms))
colnms <- c("all_male","all_female","Hispanic/Latino_male","Hispanic/Latino_female","non-Hispanic White_male","non-Hispanic White_female","non-Hispanic Black_male","non-Hispanic Black_female")
for (i in 1:length(fnms)){
  tmp <- read_excel(paste0(dir,"USLifeTables/",fnms[i]),range = "A3:G104")
  tmp <- tmp[,c(1,ncol(tmp))]
  names(tmp) <- c("age",colnms[i])
  tmp$age <- as.numeric(sub("–.*| .*","",tmp$age))
  llt[[i]] <- tmp[,2]
}

lt <- cbind(age = tmp$age,do.call(cbind,llt))

lt_long <- melt(lt,id.vars = "age",value.name = "life_expectancy")
lt_long[,c("race_ethnicity","sex")] <- do.call(rbind,strsplit(as.character(lt_long$variable),"_"))
lt_long <- lt_long[,c("age","sex","race_ethnicity","life_expectancy")]

# Use overall average life expectancies for "Other" race/ethnicity group
lt_long$race_ethnicity[lt_long$race_ethnicity=="all"] <- "Other"

# Save compiled life table
write.csv(lt_long,"US_life_table.csv",row.names = F)
