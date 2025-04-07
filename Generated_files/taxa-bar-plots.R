#### Install packages ####
# Install and load the necessary packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(ape)
library(cowplot)



#### Load the data ####

# 1. Read in OTU table, sample meta data, taxonomy table, and phylogenetic tree
# Change file paths as necessary

metafp <- "ulcers_export_div/ulcers_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "ulcers_export_div/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "ulcers_export_div/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "ulcers_export_div/tree.nwk"
phylotree <- read.tree(phylotreefp)



#### Create phyloseq object ####

# 1. Format OTU table

# OTU tables should be a matrix
# with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 


# 2. Format sample metadata

# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)


# 3. Format taxonomy table

# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)

# 4. Merge all into a phyloseq object
ulcers <- phyloseq(OTU, SAMP, TAX, phylotree)



#### Filter the data ####

# 1. Remove samples where Spaceflight is "Not Applicable"
ulcers_filt <- subset_samples(ulcers, Spaceflight !="No Applicable")

# 2. Remove non-bacterial sequences, if any
ulcers_final <- subset_taxa(ulcers_filt,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")



#### Saving phyloseq object #####

# You could save the final phyloseq object and load it in a new R script
# or continue using the object in the same script
save(ulcers_final, file="ulcers_final.RData")



#### Creating taxonomy bar plots ####

# 1. Convert to relative abundance
ulcers_RA <- transform_sample_counts(ulcers_final, fun=function(x) x/sum(x))


# 2v1. "glom" at taxonomic levels Phylum, Family and Genus
# this version does not sorting by abundance
ulcers_phylum <- tax_glom(ulcers_RA, taxrank = "Phylum", NArm=FALSE)
ulcers_family <- tax_glom(ulcers_RA, taxrank = "Family", NArm=FALSE)
ulcers_genus <- tax_glom(ulcers_RA, taxrank = "Genus", NArm=FALSE)


# 2v2. Order taxa groups from highest to lowest relative abundance within each sample
# Repeat the following code for each taxonomic rank

# Convert to a dataframe
ulcers_df <- psmelt(ulcers_RA)

# Sort the data
# - group by sample
# - then arrange the samples from highest to lowest relative abundance (RA) within each sample group
# - then reassign Phylumn column as factor and set the level order based on the RA sorting
sortedRA_phylum <- ulcers_df %>%
  group_by(Sample) %>%   
  arrange(desc(Abundance), .by_group = TRUE) %>%     
  mutate(Phylum = factor(Phylum, levels = unique(Phylum)))

# Set the factor level to reflect the sorting (for the legend)
# - Compute total relative abundance by Taxa group across all samples
# - and order the taxa groups from highest to lowest relative abundance
phylum_order <- ulcers_df %>%
  group_by(Phylum) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  pull(Phylum)

# - Then factor level globally
sortedRA_phylum$Phylum <- factor(sortedRA_phylum$Phylum, levels = phylum_order)


# 3v1a. Create taxa bar plots (with unsorted phyloseq object)

gg_phylum <- plot_bar(ulcers_phylum, fill="Phylum") + 
  facet_wrap(.~Spaceflight, scales = "free_x") +
  geom_bar(stat = "identity", position = "stack", colour = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("Relative Abundance") +
  guides(fill = guide_legend(ncol = 1)) + 
  theme(legend.position = "right", legend.justification.right="left") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14) 
  ) +
  theme(plot.margin = unit(c(1,1,1,1), "lines"))

gg_family <- plot_bar(ulcers_family, fill="Family") + 
  facet_wrap(.~Spaceflight, scales = "free_x") +
  geom_bar(stat = "identity", position = "stack", colour = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("Relative Abundance") +
  guides(fill = guide_legend(ncol = 1)) + 
  theme(legend.position = "right", legend.justification.right="left") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14) 
  ) +
  theme(plot.margin = unit(c(1,1,1,1), "lines"))

gg_genus <- plot_bar(ulcers_genus, fill="Genus") + 
  facet_wrap(.~Spaceflight, scales = "free_x") +
  geom_bar(stat = "identity", position = "stack", colour = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("Relative Abundance") +
  guides(fill = guide_legend(ncol = 1)) + 
  theme(legend.position = "right", legend.justification.right ="left") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14) 
  ) +
  theme(plot.margin = unit(c(1,1,1,1), "lines"))

# 3v1b. Create taxa bar plots (with unsorted phyloseq object), diff format
gg_genus_2 <- gg_genus +
  guides(fill = guide_legend(ncol = 3)) + 
  theme(legend.position = "bottom", legend.justification.bottom ="left")


# 3v2. Create taxa bar plots (with sorted dataframe)
# Repeat the following code for each taxonomic rank

gg_phy <- ggplot(sortedRA_phylum, aes(x = Sample, y = Abundance, fill=Phylum)) + 
  facet_wrap(.~Spaceflight, scales = "free_x") +
  geom_bar(stat = "identity", position = "stack", colour = NA, width = 0.9) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Relative Abundance") +
  guides(fill = guide_legend(ncol = 1)) + 
  theme(legend.position = "right") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14) 
  ) +
  theme(plot.margin = unit(c(1,1,1,1), "lines"))


# 4. Combine plots
gg_all <- plot_grid(gg_phylum, gg_family, gg_genus,
                    labels = c("A", "B", "C"),
                    label_size = 40,
                    align = "v",
                    hjust = -1,
                    ncol = 3)

# 5. Save plots
ggsave("plot_phylum.png"
       , gg_phylum
       , height=10, width =10)

ggsave("plot_family.png"
       , gg_family
       , height=10, width =10)

ggsave("plot_genus.png"
       , gg_genus
       , height=10, width =10)

ggsave("plot_phylum-family-genus.png"
       , gg_all
       , height=9, width =34)
