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

### working in R... is this needed? 
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

#### STRUGGLING WITH THIS SECTION OF CODE ####
# Wilcoxon test for bray-curtis
## convert to matrix format 
bray_matrix <- as.matrix(bc_dm)
# Ensure metadata contains a column 'Spaceflight' with values 'Ground Control' and 'ISS'

bray_metadata <- data.frame(
  SampleID = rownames(bray_matrix),
  Group = ifelse(bray_metadata$Spaceflight == "Ground Control", "A",
                 ifelse(bray_metadata$Spaceflight == "Space Flight", "B", NA)))  # Assign groups

# Extract pairwise distances between groups
group_A <- bray_matrix[bray_metadata$Spaceflight == "Ground Control", bray_metadata$Spaceflight == "Ground Control"]
group_B <- bray_matrix[bray_metadata$Spaceflight == "Space Flight", bray_metadata$Spaceflight == "Space Flight"]

# Remove diagonal elements (self-comparisons)

group_A <- group_A[lower.tri(bray_matrix[bray_metadata$Spaceflight == "Ground Control", bray_metadata$Spaceflight == "Ground Control"])]
group_B <- group_B[lower.tri(bray_matrix[bray_metadata$Spaceflight == "Space Flight", bray_metadata$Spaceflight == "Space Flight"])]

wilcox.test(group_A, group_B, alternative = "two.sided") #p-value = 0.7128
length(group_A)  # Should be > 0
length(group_B)  

#### STRUGGLED ABOVE #### 

#Wilcoxon test for weighted unifrac
unifrac_matrix <- phyloseq::distance(ulcers_rare, method = "wunifrac")

# Now, extract this distance matrix and the metadata from the phyloseq object

unifrac_matrix <- as.matrix(unifrac_matrix)  # Convert to matrix format
uni_metadata <- sample_data(ulcers_rare)  # Metadata from phyloseq object

# Making sure samples in uni_metadata match those in unifrac_matrix
rownames(uni_metadata) <- uni_metadata$SampleID  # Set row names of metadata to sample IDs
uni_metadata <- uni_metadata[rownames(unifrac_matrix), , drop = FALSE]  # Subset metadata based on UniFrac samples

# Extract pairwise distances between groups (Ground Control vs ISS)
group_A <- unifrac_matrix[uni_metadata$Group == "Ground Control", uni_metadata$Group == "Ground Control"]
group_B <- unifrac_matrix[uni_metadata$Group == "Space Flight", uni_metadata$Group == "Space Flight"]

# Remove diagonal elements (self-comparisons)
group_A <- group_A[lower.tri(group_A)]
group_B <- group_B[lower.tri(group_B)]

# Run the Wilcoxon test
wilcox_test_result <- wilcox.test(group_A, group_B, alternative = "two.sided")

# Output the result
wilcox_test_result
```
