# Agenda
## Interested datasets
- alcohol dependency
- multiple sclerosis
- vaping dataset
- ISS

## Potential research questions
- Is it possible to accurately predict an individual’s smoking / vaping status based on their oral microbiome profile? -> Vaping
- Can the gut microbiome be used to predict symptoms of depression?
- In patients with MS, does smoking or former smoking have an impact on microbiota composition and diversity, could this be driving disease progression in MS?
- How does microbial diversity + composition on the ISS compare to other confined environments such as simulated space habitats (HI-SEAS)? + college dormitories? -> ISS

Prevalence:
- WHY? -> smoking is very prevalent in our society
- Could be a new avenue in preventative medicine 

### Meeting notes

#### Proposed ideas
- MS -> it's being done by 2 other groups?
- smoking -> predictive model (must be able to find signficantly differences in the models for it to work)
  - didn't find a diff in vaping vs non-vaping but did find something in smoking
- depression -> another group is doing predictive model
- microbial diversity -> already done

#### Some notes
- predictive models aren't too hard?
- NASA datasets
  - they publish new datasets every day?
  - space radiation dataset not touched
- can only look at amplicon sequencing
- (take 16S reads -> map them to see who you have -> what do
  - what is function of this
  - do not use `expression` terms, use `representation`
  - finding more genes that could do it

### Potential new question
- NASA potential question 1 - OSD #487 - we can do a functional analysis (they didn't compare space vs earth - just looked at tissue types collectively)
  - they care about them getting diabetic foot ulcers
  - commonly get diabetes if you go to space
  - if diabetic foot ulcers that occur in space differ from those that naturally occuring? are theirs different? unique profile? would need a dataset that's amplicon sequencing for foot ulcers in control
  - NASA dataset has controls?
  - 1 sample sent to earth, 1 sample sent to space -> want to analyze how it progresses
  - Challenges
      - understanding very specific methodology
      - tough to compare to others
  - we'll still do:
      - taxonomic analysis
      - subset by things such as sex, age, etc
      - if have time, add functional component -> figure out how to clean data + what metadata is signficant or interested to us (ex. is space vs on earth actually significant)
      - might want to change continuously variables into categorical (ex. age 20-40 considered a `young` range) -> and compare diff groups
- 2
  - a longitudinal study ? looking at how microbiome progressed overtime - overtime in space + overtime when they returned
      - can predict it decreased overtime in space and then increased again when come back
  - looking at just skin microbiome analysis in space - OSD #572
  - can only do amplicon sequencing (search for `16S` or `amplicon` in search bar from NASA website)
  - there are plant datasets too -> soil
- 3
  - 

### TODOS - 1 more meeting (proposal due on 23rd)
[] find specific research Q w/ this dataset - make sure we really understand the experimental procedures + thinking about what we'd like to explore (3-4 aims) (ex. diversity analysis)
[] read through proposal assignment on canvas
[] be on lookout for team server
[] online meeting next week
[] archive all files on github (server isn't the most stable) - qza files
