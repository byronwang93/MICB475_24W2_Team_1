# Agenda
- Discuss concerns w/ dataset
- Discuss specific research Q w this dataset
- Discuss what we'd like to explore (3-4 aims)

## Research paper:
https://pmc.ncbi.nlm.nih.gov/articles/PMC9761711/

## Dataset:
diabetic ulcers dataset

## Concerns
- small sample size, how can we cluster?
- some of the metadata have a lot of NA fields

## Research Question Ideas
Does sex influence the microbial composition of foot ulcers? (exclude the ones in space)
  - Only do the skin swabs
  - health vs ulcer (added element of sex difference)


## Aim Ideas
1. Does sex influence the microbial composition of foot ulcers? (skin swabs only)
2. Does space conditions affect microbial composition of foot ulcers? (tissue samples only, space vs ground control)
3. Does sex influence the microbial composition of foot ulcers in space conditions? (tissue samples only, space vs ground control)
4. Functional analysis (healthy vs DFUs samples)
   - disease progression? or potential?
   - biofilm formation
   - compare with a reference genome and identify differences -> if possible describe general functions/lack of function
     - eg. we noticed that these species were missing from these samples. could affect disease potential/progression


### Meeting notes
- received data from Dr. Sun
- concerned about the small sample size
- how did they cluster -> with such a small amount
- 6 samples that remained on Earth and 6 samples sent to space
- All together -> would have few samples
- Struggling to connect all the different samples
- Separating the different types into one product or question (subset of patients)
- Can't cluster by sex -> too small of a sample size
- Even separating from sample sizes -> categories are the same (no difference)

- Not much to work with
- find something to make it work -> pivot to another dataset
- AIM 1 -> replicating the data that was already shown by the authors (similar) -> replicate their diversity metrics (COULD ADD MORE -> FOCUS ON SPACE VS. EARTH)
- AIM 2 -> (get rid of sex) replicate the abundance analysis but adjust the methods, use a different packaged to look analysis. Overall want 6 samples is enough considering the tough conditions (2 vs 4 -> problem, so want to have at least 3) -> NEW ONES SENT BY HANS
- Aldex2 and MICOM
- They only looked at relative abundance
- Two packages -> not taking relative abundance of a number of reads -> use ratios between the same microbes across samples will remove the bias of the 16s sequences. 16s relative abundance is not a good metric, but biased to abundance. Log of ratios. Take the ratio of 2 phyla and the log of that ratio and calculate it across samples will give better information about the differential abundance. Slighlty less biased information of microbes that are truly differentially abundant between the two samples. Replicating their analysis, but pratically it is a different method. Can directly compare to their findings. Use both tools. Spin it into a methodology problem. (R package)
- AIM 3 - indicator biome and core microbiome
- AIM 4 - functional analysis and aspect of the same 2 groups (space vs nonspace) -> extrapolate who is there and base the analysis over that -> see the next closest genome (picrust2)
- At least 4 figures -> core microbiome and indicator species
- Have four aims at the 6 samples that stayed on Earth and those that sent to Space
- if we wanted to directly replicate results -> would have multiple pairwise comparisons because they are distinct
- A lot of analysis of matched samples
- Good methodology -> expanding that same code to different samples
- Focus on space vs. earth (confirm the diversity and expand it by choosing different diversity metrics)
- Inspiration: https://huttenhower.sph.harvard.edu/tools/ -> all tools on 16s RNA analysis
- Introduction: Brief description of what is known and broad -> planning what we are doing. Brief overview of knowledge and literature
- Research Objective -> Aims 

### Todos

- [ ] RUN CODE THROUGH QIIME (filtering step, include only what we are working with) 
- [ ] HAVE STEPS ON HOW TO ACHIEVE ALL OF THE AIMS (STEP BY STEP WORKFLOW) -> do what we need to do and if there are any changes
