rm(list=ls())

library(sf)
library(ggplot2)
library(viridis)

# Load aggregated CDPH case and death data
load("../Data/case_and_death_demogrphcs.RData")

# Read in CA county shape files and county names
shape <- st_read("../Data/CA_counties/")
county_names <- read.csv("../Data/county_names.csv",stringsAsFactors = F)

# Function for plotting numbers and cumulative incidence
plot_num_and_inc_by_demogrphcs <- function(agg_mrg,demogrphc,vrble,fdir,shape,county_names){
  # for (i in 1:length(demogrphcs)){#
    # Cumulative numbers
    p <- ggplot(agg_mrg,aes(x=agg_mrg[,demogrphc],y=x)) + geom_bar(stat="identity") + xlab(lbls[i]) + ylab(substr(vrble,5,nchar(vrble)))
    if (demogrphc == "county_res"){
      p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1))
    }
    pdf(paste0("../Figures/",fdir,substr(vrble,5,nchar(vrble)),"_by_",demogrphc,".pdf"),width=4.375,height=3.5)
    print(p)
    dev.off()
    # Cumulative incidence
    p <- ggplot(agg_mrg,aes(x=agg_mrg[,demogrphc],y=cum_inc)) + geom_bar(stat="identity") + xlab(lbls[i]) + ylab(paste0("Cumulative incidence of ",substr(vrble,5,nchar(vrble))," (%)"))
    if (demogrphc %in% c("county_res","race_ethnicity")){
      p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1))
    }
    pdf(paste0("../Figures/",fdir,substr(vrble,1,nchar(vrble)-1),"_inc_by_",demogrphc,".pdf"),width=4.375,height=3.5)
    print(p)
    dev.off()
    # Plot map of numbers and incidence
    if (demogrphc == "county_res"){
      agg_mrg <- merge(county_names,agg_mrg,by.x="cdph_county_res",by.y="county_res",all.x=T)
      shape <- merge(shape,agg_mrg,by.x="NAME",by.y="county_res",all.x=T)
      pdf(paste0("../Figures/",fdir,substr(vrble,1,nchar(vrble)),"_map.pdf"),width=4.7,height=4)
      print(ggplot() + geom_sf(data = shape, aes(fill = x)) + scale_fill_viridis(name = substr(vrble,5,nchar(vrble)-1),discrete=FALSE) + theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank()))
      dev.off()
      pdf(paste0("../Figures/",fdir,substr(vrble,1,nchar(vrble)-1),"_inc_map.pdf"),width=4.7,height=4)
      print(ggplot() + geom_sf(data = shape, aes(fill = cum_inc)) + scale_fill_viridis(name = paste0("Cum. inc. (%)"),discrete=FALSE) + theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank()))
      dev.off()
    }
  # }
}

# Plot cases and deaths and cumulative case and death incidence by demographic factors
demogrphcs <- c("county_res","age_cat","sex","race_ethnicity")
lbls <- c("County","Age","Sex","Race/ethnicity")
fdir <- "CaseAndDeathDemographics/"
dir.create(paste0("../Figures/",fdir))

for (i in 1:length(demogrphcs)){
  if (demogrphcs[i]=="sex"){
    l$lagg_mrg[[i]]$sex <- ifelse(l$lagg_mrg[[i]]$sex==1,"male","female")
    l_deaths$lagg_mrg[[i]]$sex <- ifelse(l_deaths$lagg_mrg[[i]]$sex==1,"male","female")
  }
  # Cases
  plot_num_and_inc_by_demogrphcs(l$lagg_mrg[[i]],demogrphcs[i],"cum_cases",fdir,shape,county_names)
  # Deaths
  plot_num_and_inc_by_demogrphcs(l_deaths$lagg_mrg[[i]],demogrphcs[i],"cum_deaths",fdir,shape,county_names)
  # Print ranges of cumulative case and death incidences over demographic factor
  print(range(l$lagg_mrg[[i]]$cum_inc))
  print(range(l_deaths$lagg_mrg[[i]]$cum_inc))
}




