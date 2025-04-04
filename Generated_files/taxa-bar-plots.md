# Creating taxa bar plots
```
library(tidyverse)
library(phyloseq)
library(cowplot)

#### Load data ####
load("ulcers_final.RData")

#### Creating taxonomy bar plots ####

# Convert to relative abundance
ulcers_RA <- transform_sample_counts(ulcers_final, fun=function(x) x/sum(x))

# "glom" at taxonomic levels Genus, Family and Phylum; and remove NAs
ulcers_genus <- tax_glom(ulcers_RA, taxrank = "Genus", NArm=FALSE)
ulcers_family <- tax_glom(ulcers_RA, taxrank = "Family", NArm=FALSE)
ulcers_phylum <- tax_glom(ulcers_RA, taxrank = "Phylum", NArm=FALSE)


# Create taxa bar plots
gg_genus <- plot_bar(ulcers_genus, fill="Genus") + 
  facet_wrap(.~Spaceflight, scales = "free_x") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  ylab("Relative Abundance") +
  ggtitle("Distribution of Genera") +
  guides(fill = guide_legend(ncol = 1)) + 
  theme(
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14),
  strip.text = element_text(size = 14),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14) 
) +
  theme(legend.position = "right")
  

gg_family <- plot_bar(ulcers_family, fill="Family") + 
  facet_wrap(.~Spaceflight, scales = "free_x") +
  ylab("Relative Abundance") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  ggtitle("Distribution of Family") +
  guides(fill = guide_legend(ncol = 1)) + 
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14) 
  ) +
  theme(legend.position = "right")

gg_phylum <- plot_bar(ulcers_phylum, fill="Phylum") + 
  facet_wrap(.~Spaceflight, scales = "free_x") +
  ylab("Relative Abundance") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  ggtitle("Distribution of Phyla") +
  guides(fill = guide_legend(ncol = 1)) + 
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14) 
  ) +
  theme(legend.position = "right")


# Combine plots
gg_all <- plot_grid(gg_genus, gg_family, gg_phylum, labels = c("A", "B", "C"), ncol = 3)

ggsave("plot_genus.png"
       , gg_genus
       , height=10, width =10)

ggsave("plot_family.png"
       , gg_family
       , height=10, width =10)

ggsave("plot_phylum.png"
       , gg_phylum
       , height=10, width =10)

ggsave("plot_genus-family-phylum.png"
       , gg_all
       , height=20, width =60)
```
