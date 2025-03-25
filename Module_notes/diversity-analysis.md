### working in terminal using QIIME2 tools 
generate core metrics
```
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table spaceflight-filtered-table.qza \
  --p-sampling-depth 4008 \
  --m-metadata-file /mnt/datasets/project_2/nasa/ulcers_metadata.tsv \
  --output-dir core-metrics-results-spaceflight-filtered
```
generate Simpson Diversity 
```
qiime diversity alpha \
  --i-table spaceflight-filtered-table.qza \
  --p-metric simpson \
  --o-alpha-diversity core-metrics-results-spaceflight-filtered/simpson_diversity.qza
```

### to transfer files from server to computer for use in R
```
### WORKING ON SERVER 

qiime tools export \
  --input-path /data/ulcers/rooted-tree.qza \
  --output-path ulcers_export_div
  
qiime tools export \
--input-path /data/ulcers/spaceflight-filtered-table.qza \
--output-path table_export

biom convert \
-i feature-table.biom \
--to-tsv \
-o feature-table.txt

qiime tools export \
  --input-path /data/ulcers/taxonomy.qza \
  --output-path ulcers_export_div

### to move feature-table.txt to export file = working within table_export
mv feature-table.txt /data/ulcers/ulcers_export_div

### working in /mnt/datasets/project_2/nasa
cp ulcers_metadata.tsv /data/ulcers/ulcers_export_div

### NOW ALL 4 FILES ARE IN ulcers_export_div

### WORKING ON LOCAL COMPUTER
### exporting files into computer so we can work on them in R

scp -r root@10.19.139.161:~/data/ulcers/ulcers_export_div .

```
### Creating Phyloseq Object in R
```
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(picante)
library(ggplot2)
library(ggsignif)


#### Load data ####
# Change file paths as necessary
metafp <- "ulcers_export_div/ulcers_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "ulcers_export_div/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "ulcers_export_div/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "ulcers_export_div/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# OTU tables should be a matrix
# with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU) # to check if we have a phyloseq object 

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP) # to check if we have a phyloseq object 

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>% # removes confidence column
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>% # to separate into taxonomic ranks
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
ulcers_phylo <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(ulcers_phylo)
sample_data(ulcers_phylo)
tax_table(ulcers_phylo)
phy_tree(ulcers_phylo)

#### Subset out only the samples where Spaceflight is not Not Applicable 
ulcers_phylo_final <- subset_samples(ulcers_phylo, Spaceflight !="Not Applicable")

sample_data(ulcers_phylo_final)

#### saving data after filtering
save(ulcers_phylo_final, file="ulcers_phylo_final.RData")
```

### Generate Alpha Diversity Metrics and Plots
```
#### rarefaction

rarefaction_curve <- rarecurve(t(as.data.frame(otu_table(ulcers_phylo_final))), cex=0.1) # gives us a rarefaction curve 

ulcers_rare <- rarefy_even_depth(ulcers_phylo_final, rngseed = 1, sample.size = 4008) # sampling depth was decided during data wrangling by Byron

#### save rarefaction data
save(ulcers_rare, file="ulcers_rare.RData")

#### Alpha diversity - Shannon and Simpson #####
plot_richness(ulcers_rare) 

# only want to look at Shannon and Simpson 
plot_richness(ulcers_rare, measures = c("Shannon", "Simpson")) 

# make it into a box plot
gg_richness <- plot_richness(ulcers_rare, x = "Spaceflight", measures = c("Shannon","Simpson")) +
  xlab("Spaceflight") +
  geom_boxplot()
gg_richness

# save image
ggsave(filename = "shannon_simpson.png"
       , gg_richness
       , height=4, width=6)

# calculate Faith's phylogenetic diversity as PD
# measures phylogenetic distance 
phylo_dist <- pd(t(otu_table(ulcers_rare)), phy_tree(ulcers_rare),
                 include.root=F) 

# add PD to metadata table
sample_data(ulcers_rare)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(ulcers_rare), aes(Spaceflight, PD)) + 
  geom_boxplot() +
  xlab("Spaceflight") +
  ylab("Faith's PD")

# view plot
plot.pd

# save plot
ggsave(filename = "faith's_pd.png"
       , plot.pd
       , height=4, width=6)
```

### Generating Beta Diversity Metrics
```
### bray-curtis

bc_dm <- distance(ulcers_rare, method="bray")

# generate coordinate system 
pcoa_bc <- ordinate(ulcers_rare, method="PCoA", distance=bc_dm)

# take our rarefied phyloseq objects and combine our coordinates that we just generated 
# for different points, the body size will be identified based on colour 
# for different points, the subject will be identified based on shape 
plot_ordination(ulcers_rare, pcoa_bc, color = "Spaceflight")

gg_pcoa_bc <- plot_ordination(ulcers_rare, pcoa_bc, color = "Spaceflight") + # saving it under a variable
  labs(col = "Spaceflight") # rename legend names;col = colour channel
gg_pcoa_bc

ggsave("bray_curtis_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)

### weighted unifrac
wuf_dm <- distance(ulcers_rare, method = "unifrac", weighted = TRUE)

pcoa_wuf <- ordinate(ulcers_rare, method="PCoA", distance=wuf_dm)

plot_ordination(ulcers_rare, pcoa_wuf, color = "Spaceflight")

gg_pcoa_wuf <- plot_ordination(ulcers_rare, pcoa_wuf, color = "Spaceflight") + # saving it under a variable
  labs(col = "Spaceflight") # rename legend names;col = colour channel
gg_pcoa_wuf

ggsave("weightedunifrac_pcoa.png"
       , gg_pcoa_wuf
       , height=4, width=5)
```
### Wilcoxon Rank Sums for Alpha Diversity Metrics
```
#### WILCOXON RANK SUMS 

# generate table of alpha diversity metrics 
alphadiv <- estimate_richness(ulcers_rare)

# look only at metadata from phyloseq object
samp_dat <- sample_data(ulcers_rare)

# generate table combining metadata with alpha diversity metrics = now we have our predictor + response in one table! 
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

# Wilcoxon test for Shannon
wilcox.test(Shannon ~ Spaceflight, data=samp_dat_wdiv, exact = FALSE) ## gives us p-value = 0.02819

# Graph with Wilcoxon for Shannon
shan_spaceflight <- ggplot(samp_dat_wdiv, aes(x=`Spaceflight`, y=Shannon)) +
  geom_boxplot() +
  geom_jitter(width = 0, height = 0) +  # Add jittered dots
  geom_signif(comparisons = list(c("Space Flight","Ground Control")),
              y_position = c(2.8),
              annotations = c("p=0.02819"))+
  ylim(1.6,max(2.9))+
                       ggtitle("Shannon Diversity of Diabetic Foot Ulcers Microbiome Samples\nAcross Space and Ground Control Samples")

shan_spaceflight

ggsave("Shannon Diversity with Wilcoxon.png", plot = shan_spaceflight, width = 6, height = 4, dpi = 300)

# Wilcoxon test for Simpson
wilcox.test(Simpson ~ Spaceflight, data=samp_dat_wdiv, exact = FALSE) ## gives us p-value = 0.02819

# Graph with Wilcoxon for Simpson
simp_spaceflight <- ggplot(samp_dat_wdiv, aes(x=`Spaceflight`, y=Simpson)) +
  geom_boxplot() +
  geom_jitter(width = 0, height = 0) +  # Add jittered dots
  geom_signif(comparisons = list(c("Space Flight","Ground Control")),
              y_position = c(0.93),
              annotations = c("p=0.02819"))+
  ylim(0.75,max(0.95))+
  ggtitle("Simpson Diversity of Diabetic Foot Ulcers Microbiome Samples\nAcross Space and Ground Control Samples")

simp_spaceflight

ggsave("Simpson Diversity with Wilcoxon.png", plot = simp_spaceflight, width = 6, height = 4, dpi = 300)

# Wilcoxon test for Faith's PD
wilcox.test(PD ~ Spaceflight, data=samp_dat_wdiv, exact = FALSE) ## gives us p-value = 0.05182

# Graph with Wilcoxon for Faith's PD
PD_spaceflight <- ggplot(samp_dat_wdiv, aes(x=`Spaceflight`, y=`PD`)) +
  geom_boxplot() +
  geom_jitter(width = 0, height = 0) +  # Add jittered dots
  geom_signif(comparisons = list(c("Space Flight","Ground Control")),
              y_position = c(3.6),
              annotations = c("p=0.05182"))+
  ylim(0.2,max(4))+
  ggtitle("Faith's PD of Diabetic Foot Ulcers Microbiome Samples\nAcross Space and Ground Control Samples")

PD_spaceflight

ggsave("Faith's PD with Wilcoxon.png", plot = PD_spaceflight, width = 6, height = 4, dpi = 300)
```
### Wilcoxon Rank Sums for Beta Diversity Metrics
```

# Wilcoxon test for bray-curtis
## convert to matrix format 
bray_matrix <- as.matrix(bc_dm)

# Get sample metadata from the phyloseq object
bray_metadata <- data.frame(sample_data(ulcers_rare))

# Make sure metadata and distance matrix rows match
bray_metadata <- bray_metadata[rownames(bray_matrix), ]

# Identify groups
bray_group_A_samples <- rownames(bray_metadata[bray_metadata$Spaceflight == "Ground Control", ])
bray_group_B_samples <- rownames(bray_metadata[bray_metadata$Spaceflight == "Space Flight", ])

# Compute pairwise **within-group** distances
bray_group_A_dist <- bray_matrix[bray_group_A_samples, bray_group_A_samples][lower.tri(bray_matrix[bray_group_A_samples, bray_group_A_samples])]
bray_group_B_dist <- bray_matrix[bray_group_B_samples, bray_group_B_samples][lower.tri(bray_matrix[bray_group_B_samples, bray_group_B_samples])]

wilcox.test(bray_group_A_dist, bray_group_B_dist, alternative = "two.sided") #p-value = 0.009116

# Wilcoxon test for weighted unifrac
## convert to matrix format 
uni_matrix <- as.matrix(wuf_dm)

# Get sample metadata from the phyloseq object
uni_metadata <- data.frame(sample_data(ulcers_rare))

# Make sure metadata and distance matrix rows match
uni_metadata <- uni_metadata[rownames(uni_matrix), ]

# Identify groups
uni_group_A_samples <- rownames(uni_metadata[uni_metadata$Spaceflight == "Ground Control", ])
uni_group_B_samples <- rownames(uni_metadata[uni_metadata$Spaceflight == "Space Flight", ])

# Compute pairwise **within-group** distances
uni_group_A_dist <- uni_matrix[uni_group_A_samples, uni_group_A_samples][lower.tri(uni_matrix[uni_group_A_samples, uni_group_A_samples])]
uni_group_B_dist <- uni_matrix[uni_group_B_samples, uni_group_B_samples][lower.tri(uni_matrix[uni_group_B_samples, uni_group_B_samples])]

wilcox.test(uni_group_A_dist, uni_group_B_dist, alternative = "two.sided") #p-value = 0.01716
```
