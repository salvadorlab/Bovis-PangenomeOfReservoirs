#Noah A Legall
#Salvador Lab
#Description: An Rscript that will be used to generate necessary figures for this project.
#Last Update: May 6th

###0. libraries
install.packages("ggpubr")
library("ggpubr")
library(dplyr)

###1. Description of the data we have.
accession_meta <- read.csv("ProjectBovisReservoir_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
accession_bar <- accession_meta %>% group_by(Species, Host) %>%summarise(no_in_data = length(Species))
accession_bar <- accession_bar[order(accession_bar$Host,accession_bar$no_in_data),]
#A bar plot separated by species, colored by reservoir status.
ggbarplot(accession_bar,"Species","no_in_data", fill = "Host", color = "Host",palette = c("steelblue","red"), xlab = "Host Species", ylab = "No. of isolates", lab.size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size = 15), axis.title = element_text(size = 18)) + scale_y_log10()