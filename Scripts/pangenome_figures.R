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

###4. Heatmap of data 
gene_pres_abs <- read.csv("mbovis_prab.csv", header = TRUE, stringsAsFactors = FALSE, row.names = "Gene")
accessory_genome <- gene_pres_abs[!(is.na(gene_pres_abs$Accessory.Fragment)),]
core_genome <- gene_pres_abs[is.na(gene_pres_abs$Accessory.Fragment),]

accessory_genome_non_unique <- accessory_genome[accessory_genome$No..isolates > 1 & accessory_genome$Annotation != "hypothetical protein",]
auxil <- gene_pres_abs %>% select(2:14)

accessory_pa <- accessory_genome_non_unique %>% select(14:(ncol(accessory_genome_non_unique)))
##edit the gene presence absence to be numeric.
accessory_pa[!(accessory_pa=="")] <- 1
accessory_pa[accessory_pa==""] <- 0

##let's transpose the dataframe by turning it into a matrix first. 
pa_transpose <- t(data.matrix(accessory_pa))
heatmap(pa_transpose, scale = "none", Rowv = NA, Colv = NA, col = c("white","blue"), main = "Accessory Genome Composition")
##Would be nice to merge a core genome phylogeny to this metadata 

pa_dist_jaccard <- proxy::dist(pa_transpose,method = "Jaccard")
pa_dist_euclid <- dist(pa_transpose)

mthds <- c( "average", "single", "complete", "ward")
names(mthds) <- c( "average", "single", "complete", "ward")

jaccard_ac <- function(x) {
  cluster::agnes(pa_dist_jaccard, method = x)$ac
}

euclidean_ac <- function(x) {
  cluster::agnes(pa_dist_euclid, method = x)$ac
}


best_linkage_euclidean <- map_dbl(mthds, euclidean_ac)
best_linkage_jaccard <- map_dbl(mthds, jaccard_ac)

prab_hc <- hclust(pa_dist_jaccard, method = "ward.D" )
plot(prab_hc, cex = 0.4)
rect.hclust(prab_hc, k = 3, border = 2:5)
prab_hc_viz <- plot(prab_hc, cex = 0.4);rect.hclust(prab_hc, k = 3, border = 2:5)

mbov_meta <- read.csv("ProjectBovisReservoir_metadata.csv",header = TRUE)

prab_dend <- as.dendrogram(prab_hc) %>% color_branches(k=3) %>% hang.dendrogram(hang_height=0.1)
prab_dend_labels <- as.data.frame(labels(prab_dend)) %>% left_join(mbov_meta,by=c("labels(prab_dend)" = "??..Sample"))
prab_species <- prab_dend_labels$Host

labels_colors(prab_dend) <- rainbow_hcl(2)[sort_levels_values(as.numeric(prab_dend_labels$Host))]
labels(prab_dend) <- prab_dend_labels$Host
labels_cex(prab_dend) <- 0.4
reservoir_viz <- plot(prab_dend, horiz = TRUE)

prab_dend <- as.dendrogram(prab_hc) %>% color_branches(k=3) %>% hang.dendrogram(hang_height=0.1)


mbov_meta <- read.csv("filtered_isolate_list.csv",header = TRUE)
logsvd_model = logisticSVD(pa_transpose, k = 3)

hostgg <- plot(logsvd_model, type = "scores") + 
            geom_point(aes(colour = mbov_meta$Host)) + 
            ggtitle("M. bovis Accessory Genome Clustering") +
            theme_bw()
hostgg$labels$colour <- "Reservoir Status"

countrygg <- plot(logsvd_model, type = "scores") + geom_point(aes(colour = mbov_meta$Country)) + ggtitle("M. bovis Accessory Genome Clustering")
countrygg$labels$colour <- "Country of Origin"

wildlifegg <- plot(logsvd_model, type = "scores") + geom_point(aes(colour = mbov_meta$Species)) + ggtitle("M. bovis Accessory Genome Clustering")
wildlifegg$labels$colour <- "Species"

pdf("logisticPCA_AccessoryGenome.pdf")
hostgg
countrygg
wildlifegg
dev.off()


mbov.pca <- prcomp(pa_transpose,center = TRUE,scale. = TRUE)
plot(mbov.pca, type = "scores") + geom_point(aes(colour = mbov_meta$Species))
fviz_pca_ind(mbov.pca,geom = "point",habillage = mbov_meta$Species) + xlim(-10,20)

mbov.pcoa <- cmdscale(pa_dist_jaccard,k=3)
mbov.pcoa <- as.data.frame(mbov.pcoa)
mbov.pcoa$Host <- mbov_meta$Host
mbov.pcoa$Country <- mbov_meta$Country
mbov.pcoa$Species <- mbov_meta$Species
names(mbov.pcoa)[1:3] <- c('PC1', 'PC2','PC3')

Tr_PcoA <- ggplot(mbov.pcoa, aes(x = PC2, y = PC3, colour = Country, 
                                 label = row.names(mbov.pcoa)))
Tr_PcoA + 
  geom_point() +
  ggtitle("PCoA of M. bovis Accessory Genome Distance Matrix") +
  theme_bw()

plot_ly(x=mbov.pcoa$PC1, y=mbov.pcoa$PC2, z=mbov.pcoa$PC3, type="scatter3d", mode="markers", color=mbov.pcoa$Host)