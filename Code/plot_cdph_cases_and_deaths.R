rm(list=ls())

library(ggplot2)

# Load aggregated CDPH case and death data
load("../Data/case_and_death_demogrphcs.RData")
  
# Function for plotting numbers and cumulative incidence
plot_num_and_inc_by_demogrphcs <- function(agg_y,agg_mrg,vrble,fdir){
  # for (i in 1:length(demogrphcs)){#
    # Cumulative numbers
    p <- ggplot(agg_y,aes(x=agg_y[,demogrphcs[i]],y=x)) + geom_bar(stat="identity") + xlab(lbls[i]) + ylab(substr(vrble,5,nchar(vrble)))
    if (demogrphcs[i] == "county_res"){
      p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1))
    }
    pdf(paste0("../Figures/",fdir,substr(vrble,5,nchar(vrble)),"_by_",demogrphcs[i],".pdf"),width=6,height=4)
    print(p)
    dev.off()
    # Cumulative incidence
    p <- ggplot(agg_mrg,aes(x=agg_mrg[,demogrphcs[i]],y=cum_inc)) + geom_bar(stat="identity") + xlab(lbls[i]) + ylab(paste0("Cumulative ",substr(vrble,5,nchar(vrble)-1)," incidence (%)"))
    if (demogrphcs[i] == "county_res"){
      p <- p + theme(axis.text.x = element_text(angle = 45,hjust = 1))
    }
    pdf(paste0("../Figures/",fdir,substr(vrble,1,nchar(vrble)-1),"_inc_by_",demogrphcs[i],".pdf"),width=6,height=4)
    print(p)
    dev.off()
  # }
}

# Plot cases and deaths and cumulative case and death incidence by demographic factors
demogrphcs <- c("county_res","age_cat","sex","race_ethnicity")
lbls <- c("County","Age","Sex","Race/ethnicity")
fdir <- "CaseAndDeathDemographics/"
dir.create(paste0("../Figures/",fdir))

for (i in 1:length(demogrphcs)){
  if (demogrphcs[i]=="sex"){
    l$lagg_y[[i]]$sex <- ifelse(l$lagg_y[[i]]$sex==1,"male","female")
    l$lagg_mrg[[i]]$sex <- ifelse(l$lagg_mrg[[i]]$sex==1,"male","female")
    l_deaths$lagg_y[[i]]$sex <- ifelse(l_deaths$lagg_y[[i]]$sex==1,"male","female")
    l_deaths$lagg_mrg[[i]]$sex <- ifelse(l_deaths$lagg_mrg[[i]]$sex==1,"male","female")
  }
  # Cases
  plot_num_and_inc_by_demogrphcs(l$lagg_y[[i]],l$lagg_mrg[[i]],"cum_cases",fdir)
  # Deaths
  plot_num_and_inc_by_demogrphcs(l_deaths$lagg_y[[i]],l_deaths$lagg_mrg[[i]],"cum_deaths",fdir)
  # Print ranges of cumulative case and death incidences over demographic factor
  print(range(l$lagg_mrg[[i]]$cum_inc))
  print(range(l_deaths$lagg_mrg[[i]]$cum_inc))
}




