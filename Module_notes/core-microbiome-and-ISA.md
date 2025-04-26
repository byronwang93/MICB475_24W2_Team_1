Note: Ended up scrapping core microbiome and indicator species analysis due to insufficient results

# Creating phyloseq object and filter dataset to exclude samples where Spaceflight is Not Applicable
```
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

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
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)

#### Formatting taxonomy ####
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

#### Create phyloseq object ####
# Merge all into a phyloseq object
ulcers <- phyloseq(OTU, SAMP, TAX, phylotree)


#### Filter data ####
# Remove samples where Spaceflight is "Not Applicable"
ulcers_filt <- subset_samples(ulcers, Spaceflight !="No Applicable")
# Remove non-bacterial sequences, if any
ulcers_final <- subset_taxa(ulcers_filt,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

##### Saving #####
save(ulcers_final, file="ulcers_final.RData")
```

# Core microbiome analysis
```
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("ulcers_final.RData")

#### "core" microbiome ####

# Convert to relative abundance
ulcers_RA <- transform_sample_counts(ulcers_final, fun=function(x) x/sum(x))

# Filter dataset by Spaceflight category
ulcers_space <- subset_samples(ulcers_RA, `Spaceflight`=="Space Flight")
ulcers_ground <- subset_samples(ulcers_RA, `Spaceflight`=="Ground Control")

# Find ASVs that are found in more than 10% of samples in each Spaceflight category
space_ASVs <- core_members(ulcers_space, detection=0, prevalence = 0.1)
ground_ASVs <- core_members(ulcers_ground, detection=0, prevalence = 0.1)

# Label each ASV
tax_table(prune_taxa(space_ASVs,ulcers_final))
tax_table(prune_taxa(ground_ASVs,ulcers_final))
# Very few core microbiome members

# Plot those ASVs' relative abundance
prune_taxa(space_ASVs,ulcers_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Spaceflight`, scales ="free")

# Display results in a Venn diagram showing all the ASVs that showed up in each Spacelight category
space_list <- core_members(ulcers_space, detection=0.001, prevalence = 0.10)
ground_list <- core_members(ulcers_ground, detection=0.001, prevalence = 0.10)

spaceflight_list_full <- list(Space_Flight = space_list, Ground_Control = ground_list)

# Create a Venn diagram using all the ASVs shared and unique to Spaceflight conditions
ulcers_venn <- ggVennDiagram(x = spaceflight_list_full)

ggsave("venn_spaceflight.png", ulcers_venn)
```

# Indcator species analysis
```
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
load("ulcers_final.RData")

#### Indicator Species/Taxa Analysis ####

# glom to Genus (group by Genus because samples were resolved up to Genus level)
ulcers_genus <- tax_glom(ulcers_final, "Genus", NArm = FALSE)
ulcers_genus_RA <- transform_sample_counts(ulcers_genus, fun=function(x) x/sum(x))

#ISA
isa_ulcers <- multipatt(t(otu_table(ulcers_genus_RA)), cluster = sample_data(ulcers_genus_RA)$`Spaceflight`)
summary(isa_ulcers)
taxtable <- tax_table(ulcers_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of 
# anything beyond the glomed taxa level
ulcers_isatable <- isa_ulcers$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05)

ulcers_isatable %>% View()

#### Saving ####
save(ulcers_isatable, file="ulcers_isatable.RData")
```
