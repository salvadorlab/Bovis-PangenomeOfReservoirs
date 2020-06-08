#Noah A Legall
#Salvador Lab
#Description: An Rscript that will be used to generate necessary figures for this project.
#Last Update: May 6th

###0. libraries
install.packages("ggpubr")
install.packages("patchwork")
library("ggpubr")
library(dplyr)
library(ggplot2)
library(patchwork)
library(magrittr)

###1. Description of the data we have.
accession_meta <- read.csv("ProjectBovisReservoir_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
accession_bar <- accession_meta %>% group_by(Species, Host) %>%summarise(no_in_data = length(Species))
accession_bar <- accession_bar[order(accession_bar$Host,accession_bar$no_in_data),]
#A bar plot separated by species, colored by reservoir status.
ggbarplot(accession_bar,"Species","no_in_data", fill = "Host", color = "Host",palette = c("steelblue","red"), xlab = "Host Species", ylab = "No. of isolates", lab.size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size = 15), axis.title = element_text(size = 18)) + scale_y_log10()

###2. Distribution of N50 for the assemblies
setwd("/Users/noah_/Documents/Bovis-PangenomeOfReservoirs/Resources/")
assembly_stats <- read.csv("mbovis_reservoir_genome_stat.tsv.csv",header = TRUE)

N50_hist <- assembly_stats %>% ggplot(aes(x=N50)) +
                   geom_histogram(binwidth = 5000, fill = "steelblue", color = "steelblue") +
                   theme_light()

genome_length_hist <- assembly_stats  %>% ggplot(aes(x=Total.length)) +
  geom_histogram(binwidth = 1000, fill = "steelblue", color = "steelblue") +
  theme_light()

GC_hist <- assembly_stats %>% ggplot(aes(x=GC....)) +
  geom_histogram(binwidth = 0.005, fill = "steelblue", color = "steelblue") +
  theme_light()

genome_frac_hist <- assembly_stats  %>% ggplot(aes(x=Genome.fraction....)) + 
  geom_histogram(binwidth = 0.01, fill = "steelblue", color = "steelblue") +
  theme_light()

patch_assemblies <- (N50_hist | genome_length_hist)/(GC_hist | genome_frac_hist)

assembly_stats %<>%
  filter(N50 > 45000) %>%
  filter(Total.length < 5000000) %>%
  filter(GC.... > 65.0 & GC.... < 65.99) %>%
  filter(Genome.fraction.... > 97.0)

###3. Impact of Filtering on Data
#goal: simply, just include the isolates whose names make it past the filtering scheme I devised.
#a. apply the filters on the assembly_stats table, extract the names, and filter the metadata table based on whats available
#b. after that, I'll just recreate the Data Description figure.

#a. 
names <- gsub('.scaffold','',assembly_stats$Assembly)
filtered_accession_metadata <- accession_meta[accession_meta$??..Sample %in% names,]

#b.
accession_bar <- filtered_accession_metadata %>% group_by(Species, Host) %>%summarise(no_in_data = length(Species))
accession_bar <- accession_bar[order(accession_bar$Host,accession_bar$no_in_data),]
#A bar plot separated by species, colored by reservoir status.
ggbarplot(accession_bar,"Species","no_in_data", fill = "Host", color = "Host",palette = c("steelblue","red"), xlab = "Host Species", ylab = "No. of isolates", lab.size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size = 15), axis.title = element_text(size = 18)) + scale_y_log10()
write.csv(assembly_stats,file = "filtered_mbovis_stats.csv")
