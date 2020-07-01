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
assembly_stats <- read.csv("mbovis_transposed_report.csv",header = TRUE)

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

accessory_genome_non_unique <- accessory_genome[accessory_genome$No..isolates > 1,]
auxil <- gene_pres_abs %>% select(2:14)

accessory_pa <- accessory_genome_non_unique %>% select(14:(ncol(accessory_genome_non_unique)))
##edit the gene presence absence to be numeric.
accessory_pa[!(accessory_pa=="")] <- 1
accessory_pa[accessory_pa==""] <- 0

##let's transpose the dataframe by turning it into a matrix first. 
pa_transpose <- t(data.matrix(accessory_pa))
heatmap(pa_transpose, scale = "none", Rowv = NA, Colv = NA, col = c("white","blue"), main = "Accessory Genome Composition")
##Would be nice to merge a core genome phylogeny to this metadata 


### 5. Hierarchical Clustering w/ Dendextend
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
prab_dend_labels <- as.data.frame(labels(prab_dend)) %>% left_join(mbov_meta,by=c("labels(prab_dend)" = "Sample"))
prab_species <- prab_dend_labels$Host
prab_country <- prab_dend_labels$Country
countrydend <- rainbow_hcl(3)[prab_dend_labels$Country]
hostdend <- rainbow_hcl(2)[prab_dend_labels$Host]
bars <- cbind(countrydend,hostdend)
prab_dend <- prab_dend %>% set("labels","")



plot(prab_dend)
colored_bars(bars,prab_dend)

labels_colors(prab_dend) <- rainbow_hcl(2)[sort_levels_values(as.numeric(prab_dend_labels$Host))]
labels(prab_dend) <- prab_dend_labels$Host
labels_cex(prab_dend) <- 0.4
reservoir_viz <- plot(prab_dend, horiz = TRUE)

prab_dend <- as.dendrogram(prab_hc) %>% color_branches(k=3) %>% hang.dendrogram(hang_height=0.1)


mbov_meta <- read.csv("filtered_isolate_list.csv",header = TRUE)

#### 6. PCoA
#compute percent explained variance
var_ex <- cmdscale(pa_dist_jaccard,k = 699)
vars <- apply(var_ex, 2, var)
pcoa_var_explained <- (vars/sum(vars)) * 100
pcoa_var_explained <- pcoa_var_explained[1:10]

prin_comp <- paste0("PC",1:10)
data <- data.frame(
  PC=prin_comp,  
  value=pcoa_var_explained
)

#create dataframe with the first 3 PCs
mbov.pcoa <- cmdscale(pa_dist_jaccard,k=4)
mbov.pcoa <- as.data.frame(mbov.pcoa)
mbov.pcoa$Host <- mbov_meta$Host
mbov.pcoa$Country <- mbov_meta$Country
mbov.pcoa$Species <- mbov_meta$Species
mbov.pcoa$Cattle <- mbov_meta$Domesticated
names(mbov.pcoa)[1:4] <- c('PC1', 'PC2','PC3','PC4')

#Scree plot
pcoa_scree <- data %>% ggplot(aes(x=reorder(PC, -value), y=value)) + 
  geom_bar(stat = "identity",fill = "steelblue") +
  geom_text(aes(label = round(value,1)), color = "black", position=position_dodge(width=0.9), vjust=-0.25) +
  theme_bw() +
  ggtitle("Principal Coordinate Analysis - Scree Plot") +
  xlab("Principal Component") + 
  ylab("Percent of Explained Variance")

pdf("PCoA_Scree.pdf")
pcoa_scree
dev.off()

#PCoA with different groupings explained. 
hostgg <- ggplot(mbov.pcoa, aes(x = PC3, y = PC4, colour = Host, label = row.names(mbov.pcoa))) + 
  geom_point() + 
  theme_bw()

speciesgg <- ggplot(mbov.pcoa, aes(x = PC3, y = PC4, colour = Species, label = row.names(mbov.pcoa))) + 
  geom_point() + 
  theme_bw()

countrygg <- ggplot(mbov.pcoa, aes(x = PC3, y = PC4, colour = Country, label = row.names(mbov.pcoa))) + 
  geom_point() + 
  theme_bw()

wildlifegg <- ggplot(mbov.pcoa, aes(x = PC3, y = PC4, colour = Cattle, label = row.names(mbov.pcoa))) + 
  geom_point() + 
  theme_bw()

pdf("PCoA_AccessoryGenome.pdf")
hostgg
countrygg
speciesgg
wildlifegg
dev.off()

mbov.pcoa_patch <- (hostgg | countrygg)/(wildlifegg | speciesgg)
#These are 3D versions of the plots, good for further visualization but not necessary

#PC's 1,2,3
plot_ly(x=mbov.pcoa$PC1, y=mbov.pcoa$PC2, z=mbov.pcoa$PC3, type="scatter3d", mode="markers", color = mbov.pcoa$PC4)

#PC's 1,3,4
plot_ly(x=mbov.pcoa$PC1, y=mbov.pcoa$PC3, z=mbov.pcoa$PC4, type="scatter3d", mode="markers", color = mbov.pcoa$Species)

# 4 Dimensions + Country of Origin 
mbov.pcoa_5D <- ggplot(mbov.pcoa, aes(x = PC3, y = PC4, color = PC1, shape = Country, label = row.names(mbov.pcoa))) + 
  geom_point(aes(size = PC2)) + 
  theme_bw() +
  ggtitle("All 4 Dimensions of PCoA")

pdf("PCoA_AllDimension.pdf")
mbov.pcoa_5D
dev.off()

pcoa_patch <- mbov.pcoa_5D | pcoa_scree
#### 7.
