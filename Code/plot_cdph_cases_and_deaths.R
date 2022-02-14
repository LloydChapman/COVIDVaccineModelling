rm(list=ls())

library(sf)
library(ggplot2)
library(viridis)
library(cowplot)

# Set plot theme
theme_set(theme_cowplot(font_size = 13) + theme(
  strip.background = element_blank(),
  plot.background = element_rect(fill="white"),
  legend.background = element_rect(fill="white"),
  panel.background = element_rect(fill="white")))

# Load aggregated CDPH case and death data
load("../Data/case_and_death_demogrphcs.RData")

# Read in CA county shape files and county names
shape <- st_read("../Data/CA_counties/")
county_names <- read.csv("../Data/county_names.csv",stringsAsFactors = F)

# Function for plotting numbers and cumulative incidence
plot_num_and_inc_by_demogrphcs <- function(agg_mrg,demogrphc,vrble,fdir,shape,county_names){
  if (demogrphc != "county_res"){
    # Cumulative numbers
    p1 <- ggplot(agg_mrg,aes(x=agg_mrg[,demogrphc],y=x)) + geom_bar(stat="identity") + xlab(lbls[i]) + ylab(substr(vrble,5,nchar(vrble)))
    if (demogrphc == "county_res"){
      p1 <- p1 + theme(axis.text.x = element_text(angle = 45,hjust = 1))
    }
    # Cumulative incidence
    p2 <- ggplot(agg_mrg,aes(x=agg_mrg[,demogrphc],y=cum_inc)) + geom_bar(stat="identity") + xlab(lbls[i]) + ylab(paste0("Cumulative incidence of ",substr(vrble,5,nchar(vrble))," (%)"))
    if (demogrphc %in% c("county_res","age_cat","race_ethnicity")){
      p2 <- p2 + theme(axis.text.x = element_text(angle = 45,hjust = 1))
    }
  } else { # Plot map of numbers and incidence
    agg_mrg <- merge(county_names,agg_mrg,by.x="cdph_county_res",by.y="county_res",all.x=T)
    shape <- merge(shape,agg_mrg,by.x="NAME",by.y="county_res",all.x=T)
    p1 <- ggplot() + 
      geom_sf(data = shape, aes(fill = x)) + 
      scale_fill_viridis(name = substr(vrble,5,nchar(vrble)-1),discrete=FALSE) + 
      theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank())
    p2 <- ggplot() + 
      geom_sf(data = shape, aes(fill = cum_inc)) + 
      scale_fill_viridis(name = paste0("Cum. inc.\n",substr(vrble,5,nchar(vrble))," (%)"),discrete=FALSE) + 
      theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank())
  }
  return(list(p1,p2))
}

# Plot cases and deaths and cumulative case and death incidence by demographic factors
demogrphcs <- c("county_res","age_cat","sex","race_ethnicity")
lbls <- c("County","Age","Sex","Race/ethnicity")
fdir <- "CaseAndDeathDemographics/"
dir.create(paste0("../Figures/",fdir))

lp <- vector("list",2*length(demogrphcs)) 
for (i in 1:length(demogrphcs)){
  if (demogrphcs[i]=="sex"){
    l$lagg_mrg[[i]]$sex <- ifelse(l$lagg_mrg[[i]]$sex==1,"male","female")
    l_deaths$lagg_mrg[[i]]$sex <- ifelse(l_deaths$lagg_mrg[[i]]$sex==1,"male","female")
  }
  # Cases
  lp[[2*(i-1)+1]] <- plot_num_and_inc_by_demogrphcs(l$lagg_mrg[[i]],demogrphcs[i],"cum_cases",fdir,shape,county_names)
  # Deaths
  lp[[2*i]] <- plot_num_and_inc_by_demogrphcs(l_deaths$lagg_mrg[[i]],demogrphcs[i],"cum_deaths",fdir,shape,county_names)
  # Print ranges of cumulative case and death incidences over demographic factor
  print(range(l$lagg_mrg[[i]]$cum_inc))
  print(range(l_deaths$lagg_mrg[[i]]$cum_inc))
}

lp1 <- lapply(lp,"[[",2)
for (i in 1:2){
  lp1[[i]] <- lp1[[i]] + theme(axis.line=element_blank(),axis.ticks=element_blank())
}
p <- plot_grid(plotlist=lp1,align="v",axis="l",nrow=4,ncol=2,rel_heights = c(1.5,1.1,1,1.3),labels="AUTO",hjust=rep(c(-0.5,0.6),4))
ggsave(paste0("../Figures/",fdir,"case_and_death_inc_by_demogrphcs.pdf"),p,width=10,height=15)



