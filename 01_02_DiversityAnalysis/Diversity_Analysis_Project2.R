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

#### rarefaction

rarefaction_curve <- rarecurve(t(as.data.frame(otu_table(ulcers_phylo_final))), cex=0.1) # gives us a rarefaction curve 

set.seed(1)

ulcers_rare <- rarefy_even_depth(ulcers_phylo_final, rngseed = 1, sample.size = 4008) # sampling depth was decided during data wrangling by Byron

#### save rarefaction data

save(ulcers_rare, file="ulcers_rare.RData")

#### Alpha diversity - Shannon and Simpson #####
plot_richness(ulcers_rare) 

# only want to look at Shannon, Simpson, and Observed
plot_richness(ulcers_rare, measures = c("Shannon", "Simpson", "Observed")) 

# make it into a box plot
gg_richness <- plot_richness(ulcers_rare, x = "Spaceflight", measures = c("Shannon","Simpson", "Observed")) +
  xlab("Spaceflight") +
  geom_boxplot()
gg_richness

# save image
ggsave(filename = "shannon_simpson_observed.png"
       , gg_richness
       , height=4, width=6)

# calculate Faiths phylogenetic diversity as PD
# measures phylogenetic distance 
phylo_dist <- pd(t(otu_table(ulcers_rare)), phy_tree(ulcers_rare),
                 include.root=F) 

# add PD to metadata table
sample_data(ulcers_rare)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(ulcers_rare), aes(Spaceflight, PD, fill = Spaceflight)) + 
  geom_boxplot(colour = "black") +
  scale_fill_manual(values = c("blue", "red")) +
  xlab("Spaceflight") +
  ylab("Faith's Phylogenetic Diversity")

# view plot
plot.pd

# save plot
ggsave(filename = "faith's_pd.png"
       , plot.pd
       , height=4, width=6)

##### GENERATE BETA DIVERSITY METRICS

### bray-curtis

bc_dm <- distance(ulcers_rare, method="bray")

# generate coordinate system 
pcoa_bc <- ordinate(ulcers_rare, method="PCoA", distance=bc_dm)

# take our rarefied phyloseq objects and combine our coordinates that we just generated 

ord.bray <- plot_ordination(ulcers_rare, pcoa_bc, color = "Spaceflight")

gg_pcoa_bc <- plot_ordination(ulcers_rare, pcoa_bc, color = "Spaceflight") + 
  scale_color_manual(values = c("blue", "red")) +
  labs(col = "Spaceflight") # rename legend names;col = colour channel
gg_pcoa_bc

ggsave("bray_curtis_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)

### weighted unifrac
wuf_dm <- distance(ulcers_rare, method = "unifrac", weighted = TRUE)

pcoa_wuf <- ordinate(ulcers_rare, method="PCoA", distance=wuf_dm)

plot_ordination(ulcers_rare, pcoa_wuf, color = "Spaceflight")

gg_pcoa_wuf <- plot_ordination(ulcers_rare, pcoa_wuf, color = "Spaceflight") +
  scale_color_manual(values = c("blue", "red")) +
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
shan_spaceflight <- ggplot(samp_dat_wdiv, aes(x=`Spaceflight`, y=Shannon, fill = `Spaceflight`)) +
  geom_boxplot(colour = "black") +
  geom_jitter(width = 0, height = 0) +  # Add jittered dots
  geom_signif(comparisons = list(c("Space Flight","Ground Control")),
              y_position = c(2.8),
              annotations = c("p=0.03"), color = "black")+
  ylim(1.6,max(2.9))+
                       ggtitle("Shannon Diversity of Diabetic Foot Ulcers Microbiome Samples\nAcross Space and Ground Control Samples") +
  scale_fill_manual(values = c("Space Flight" = "red", "Ground Control" = "blue")) +
  scale_color_manual(values = c("Space Flight" = "red", "Ground Control" = "blue")) 

shan_spaceflight

ggsave("Shannon_Wilcoxon_Final.png", plot = shan_spaceflight, width = 6, height = 4, dpi = 300)

# Wilcoxon test for Simpson
wilcox.test(Simpson ~ Spaceflight, data=samp_dat_wdiv, exact = FALSE) ## gives us p-value = 0.02819

# Graph with Wilcoxon for Simpson
simp_spaceflight <- ggplot(samp_dat_wdiv, aes(x=`Spaceflight`, y=Simpson, fill = `Spaceflight`)) +
  geom_boxplot(colour = "black") +
  geom_jitter(width = 0, height = 0) +  # Add jittered dots
  geom_signif(comparisons = list(c("Space Flight","Ground Control")),
              y_position = c(0.93),
              annotations = c("p=0.03"), colour = "black")+
  ylim(0.75,max(0.95))+
  ggtitle("Simpson Diversity of Diabetic Foot Ulcers Microbiome Samples\nAcross Space and Ground Control Samples") +
  scale_fill_manual(values = c("Space Flight" = "red", "Ground Control" = "blue")) +
  scale_color_manual(values = c("Space Flight" = "red", "Ground Control" = "blue")) 

simp_spaceflight

ggsave("Simpson_Wilcoxon_Final.png", plot = simp_spaceflight, width = 6, height = 4, dpi = 300)

# Wilcoxon test for Observed
wilcox.test(Observed ~ Spaceflight, data=samp_dat_wdiv, exact = FALSE) ## gives us p-value = 0.04986

# Graph with Wilcoxon for Observed
obs_spaceflight <- ggplot(samp_dat_wdiv, aes(x=`Spaceflight`, y=Observed, fill = `Spaceflight`)) +
  geom_boxplot(colour = "black") +
  geom_jitter(width = 0, height = 0) +  # Add jittered dots
  geom_signif(comparisons = list(c("Space Flight","Ground Control")),
              y_position = c(23.3),
              annotations = c("p=0.05"), colour = "black")+
  ylim(8,max(24))+
  ggtitle("Observed Diversity of Diabetic Foot Ulcers Microbiome Samples\nAcross Space and Ground Control Samples") +
  scale_fill_manual(values = c("Space Flight" = "red", "Ground Control" = "blue")) +
  scale_color_manual(values = c("Space Flight" = "red", "Ground Control" = "blue")) 

obs_spaceflight

ggsave("Observed_Wilcoxon_Final.png", plot = obs_spaceflight, width = 6, height = 4, dpi = 300)

# Wilcoxon test for Faith's PD
wilcox.test(PD ~ Spaceflight, data=samp_dat_wdiv, exact = FALSE) ## gives us p-value = 0.05182

# Graph with Wilcoxon for Faith's PD
PD_spaceflight <- ggplot(samp_dat_wdiv, aes(x=`Spaceflight`, y=`PD`, fill = `Spaceflight`)) +
  geom_boxplot(colour = "black") +
  geom_jitter(width = 0, height = 0) +  # Add jittered dots
  geom_signif(comparisons = list(c("Space Flight","Ground Control")),
              y_position = c(3.6),
              annotations = c("p=0.05"), colour = "black")+
  ylim(0.2,max(4))+
  ggtitle("Faith's PD of Diabetic Foot Ulcers Microbiome Samples\nAcross Space and Ground Control Samples") +
  scale_fill_manual(values = c("Space Flight" = "red", "Ground Control" = "blue")) +
  scale_color_manual(values = c("Space Flight" = "red", "Ground Control" = "blue")) +
  ylab("Faith's Phylogenetic Diversity") 

PD_spaceflight

ggsave("PD_Wilcoxon_Final.png", plot = PD_spaceflight, width = 6, height = 4, dpi = 300)


#### PERMANOVA FOR BETA DIVERSITY 

dm_braycurtis <- vegdist(t(otu_table(ulcers_rare)), method="bray") # Bray-curtis

set.seed(123)
adonis2(dm_braycurtis ~ `Spaceflight`, data=samp_dat_wdiv) # p-value = 0.707

# Replot Bray-Curtis
ord.bray <- ordinate(ulcers_rare, method = "PCoA", distance = "bray")

BC_spaceflight <- plot_ordination(ulcers_rare, ord.bray, color = "Spaceflight") +
  scale_color_manual(values = c("blue", "red")) +
  stat_ellipse(type = "norm") + 
  ggtitle("Bray-Curtis for Diabetic Foot Ulcers Microbiome Samples\nAcross Space and Ground Control Samples")

BC_spaceflight
ggsave("bray_curtis_permanova_final.png", plot = BC_spaceflight, width = 6, height = 4, dpi = 300)

#Weighted UniFrac
dm_unifrac <- UniFrac(ulcers_rare, weighted = TRUE, normalized = TRUE, fast = TRUE)

ord.unifrac <- ordinate(ulcers_rare, method="PCoA", distance="dm_unifrac")

plot_ordination(ulcers_rare, ord.unifrac, color="Spaceflight")

set.seed(123)
adonis2(dm_unifrac ~ `Spaceflight`, data=samp_dat_wdiv) # p-value = 0.724

# Replot Weighted UniFrac

uni_spaceflight <- plot_ordination(ulcers_rare, ord.unifrac, color = "Spaceflight") +
  scale_color_manual(values = c("blue", "red")) +
  stat_ellipse(type = "norm") +
  ggtitle("Weighted UniFrac for Diabetic Foot Ulcers Microbiome Samples\nAcross Space and Ground Control Samples")

uni_spaceflight

ggsave("weighted_unifrac_permanova_final.png", plot = uni_spaceflight, width = 6, height = 4, dpi = 300)

#Unweighted UniFrac
dm_ununifrac <- UniFrac(ulcers_rare, weighted = FALSE, fast = TRUE) 

ord.ununifrac <- ordinate(ulcers_rare, method="PCoA", distance="dm_ununifrac")

plot_ordination(ulcers_rare, ord.ununifrac, color="Spaceflight")

set.seed(123)
adonis2(dm_ununifrac ~ `Spaceflight`, data=samp_dat_wdiv) # p-value = 0.034

ununi_spaceflight <- plot_ordination(ulcers_rare, ord.ununifrac, color = "Spaceflight") +
  scale_color_manual(values = c("blue", "red")) +
  stat_ellipse(type = "norm") +
  ggtitle("UniFrac for Diabetic Foot Ulcers Microbiome Samples\nAcross Space and Ground Control Samples")

ununi_spaceflight

ggsave("unweighted_unifrac_permanova_final.png", plot = ununi_spaceflight, width = 6, height = 4, dpi = 300)

