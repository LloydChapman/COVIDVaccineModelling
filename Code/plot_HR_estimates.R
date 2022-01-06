rm(list=ls())

library(sf)
library(ggplot2)
library(viridis)
library(cowplot)

# Set plot theme
theme_set(theme_cowplot(font_size = 11) + theme(
  strip.background = element_blank(),
  plot.background = element_rect(fill="white"),
  legend.background = element_rect(fill="white"),
  panel.background = element_rect(fill="white")))

# Read in CA county shape files
shape <- st_read("../Data/CA_Counties/")

# Read in hazard ratio estimates
HRs <- read.csv("../Data/coeffs_death_2020-02-05_3.csv",stringsAsFactors = F)

## County HRs
# Subset to county HRs
HRs_county <- HRs[grep("county_res",HRs$term),]
# Remove "county_res" from county names
HRs_county$term <- sub("county_res","",HRs_county$term)

# Read in and merge county names
county_names <- read.csv("../Data/county_names.csv",stringsAsFactors = F)
HRs_county <- merge(county_names,HRs_county,by.x="cdph_county_res",by.y="term",all.x=T)
# Fill in data for Alameda
HRs_county[HRs_county$cdph_county_res=="Alameda",3:ncol(HRs_county)] <- list(0,NA,NA,NA,1,NA,NA,"1 (Ref)")

# Merge with shape file data frame
shape <- merge(shape,HRs_county,by.x="NAME",by.y="county_res",all.x=T)

# Plot
fdir <- "../Figures/HRplots/"
dir.create(fdir)
p1 <- ggplot() + geom_sf(data = shape, aes(fill = exp_estimate)) + 
  scale_fill_viridis(name = "HR",discrete=FALSE) + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.line=element_blank(),axis.ticks=element_blank())

## Age HRs
HRs_age <- HRs[grep("age_cat",HRs$term),]
HRs_age$term <- sub("age_cat","",HRs_age$term)
# Fill in data for <10 year-olds
HRs_age <- rbind(list("<10",0,NA,NA,NA,1,NA,NA,"1 (Ref)"),HRs_age)
p2 <- ggplot(HRs_age) + geom_point(aes(x=term,y=exp_estimate)) + 
  geom_errorbar(aes(x=term,ymin=X2.5..,ymax=X97.5..),width=0.5,position=position_dodge(0.9)) + 
  xlab("Age (years)") + ylab("HR") +
  theme(axis.text.x=element_text(angle=45,hjust=1))

## Sex HRs
HRs_sex <- HRs[HRs$term=="sex",]
HRs_sex$term <- "male"
# Fill in data for <10 year-olds
HRs_sex <- rbind(list("female",0,NA,NA,NA,1,NA,NA,"1 (Ref)"),HRs_sex)
p3 <- ggplot(HRs_sex,aes(x=term)) + geom_point(aes(y=exp_estimate)) + 
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..),width=0.5,position=position_dodge(0.9)) + 
  xlab("Sex") + ylab("HR") + ylim(c(0,max(HRs_sex$X97.5..)))

## Race-ethnicity HRs
HRs_race_ethnicity <- HRs[grep("race_ethnicity",HRs$term),]
HRs_race_ethnicity$term <- sub("race_ethnicity","",HRs_race_ethnicity$term)
# Fill in data for <10 year-olds
HRs_race_ethnicity <- rbind(list("Hispanic/Latino",0,NA,NA,NA,1,NA,NA,"1 (Ref)"),HRs_race_ethnicity)
p4 <- ggplot(HRs_race_ethnicity,aes(x=term)) + geom_point(aes(y=exp_estimate)) + 
  geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..),width=0.5,position=position_dodge(0.9)) + 
  xlab("Race/ethnicity") + ylab("HR") + ylim(0,1) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggsave(paste0(fdir,"HR_plots.pdf"),plot_grid(p1,p2,p3,p4,nrow=2,ncol=2,align="h",axis="b",rel_widths=c(1,1.1),labels="AUTO"),width=7.5,height=8)
