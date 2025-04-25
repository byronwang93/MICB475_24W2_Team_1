## Important notes
- a few files such as `ulcers-demux.qza` is too large to be stored on github (2.07GB), will look into potentially storing zipped files onto this github

## Steps

1. created `ulcers` directory in `/data`
2. generating the `ulcers-demux.qza` file
   - a few problems w/ our data
      - manifest file had columns for **paired-end** data but QIIME 2 initially expected a **single-end** manifest format
         - :. we switched to using V2 manifest format (`PairedEndFastqManifestPhred33V2`)
      - our manifest file was in Windows format (CRLF) so QIIME2 couldn't interpret header line
         - we converted the file to Unix line endings (LF) using `sed` after first copying it to a writable location
```bash
# copy to writable location
cp /mnt/datasets/project_2/nasa/ulcers_manifest.tsv /data/ulcers/ulcers_manifest_fixed.tsv

# change headings to be LF instead of CRLF (for QIIME compatability)
sed -i 's/\r$//' /data/ulcers/ulcers_manifest_fixed.tsv

# took raw paired-end FASTQ files -> packaged into a QIIME2 artifact
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /mnt/datasets/project_2/nasa/ulcers_manifest.tsv \
  --output-path ulcers-demux.qza \
  --input-format PairedEndFastqManifestPhred33
```
4. creating visualization of demultiplexed samples - `ulcers-demux.qzv`
```bash
qiime demux summarize \
--i-data ulcers-demux.qza \
--o-visualization ulcers-demux.qzv
```
6. denoising and clustering step - `table.qza` `rep-seqs.qza` `denoising-stats.qza`
   - taking over an hour to run (decided to run in separate screen)
```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ulcers-demux.qza \
  --p-trim-left-f 0 \
  --p-trunc-len-f 205 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 205 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
```
7. visualizing ASV stats - `stats.qzv` `table.qzv` `rep-seqs.qzv`
```bash
# Visualize DADA2 stats
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

# Visualize ASVs stats
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/nasa/ulcers_metadata.tsv
  
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
```
8. taxonomic analysis - `ref-seqs-trimmed-qza` `classifier.qza` `taxonomy.qza` `taxa-bar-plots.qzv`
- background info
    - now that we know what our unique ASVs are, want to know which of these correspond to what taxonomic groups, species, etc
    - `rep-seqs` gives you raw files
        - instead of blasting each separately, will use QIIME2
    - will use `SILVA` database
        - containing full 16S seqs of almost all sequences we have up to date - trim it + truncate it based on part of 16S seq that is part of your dataset
```bash
# use silva-138-99-nb-classifier.qza as your database
# -> represents full 16S sequences
  
# Extract your amplicon of interest from the reference database
#replace the ref-otus.qza with the representative sequence file on the server (e.g. /mnt/datasets/silva_ref_files/silva-138-99-seqs.qza)
#replace primer sequences with your correct sequences
#replace trunc-len with the one you defined in your denoising step
qiime feature-classifier extract-reads \
  --i-sequences /mnt/datasets/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 230 \
  --o-reads ref-seqs-trimmed.qza

# Train classifier with your new ref-seq file
# Replace ref-taxonomy.qza with the representative taxonomy file on the server (e.g. /mnt/datasets/silva_ref_files/silva-138-99-tax.qza)
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-trimmed.qza \
  --i-reference-taxonomy /mnt/datasets/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier classifier.qza

# Use the trained classifier to assign taxonomy to your reads (rep-seqs.qza)
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
  
# Taxonomy barplots
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file /mnt/datasets/project_2/nasa/ulcers_metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```
9. filtering based on taxonomy, frequency or metadata - `table-no-mitochondria-no-chloroplast.qza` `feature-frequency-filtered-table.qza` `spaceflight-filtered-table.qza`
- remove mitochondria + chloroplast
```bash
# remove mitochondria + chloroplast seqs
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

qiime feature-table filter-features \
--i-table table-no-mitochondria-no-chloroplast.qza \
--p-min-frequency 62 \
--o-filtered-table feature-frequency-filtered-table.qza

# NEED TO DO META DATA FILTERING? -> filter for Space Flight and Ground Control
qiime feature-table filter-samples \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --m-metadata-file /mnt/datasets/project_2/nasa/ulcers_metadata.tsv \
  --p-where "[Spaceflight] != 'Not Applicable'" \
  --o-filtered-table spaceflight-filtered-table.qza
```
10. MPT alpha rarefaction step - `table-spaceflight-filtered-and-no-mitochondria-no-chloroplast.qzv` `table-no-mitochondria-no-chloroplast.qzv` `aligned-rep-seqs.qza masked-aligned-rep-seqs.qza` `unrooted-tree.qza` `rooted-tree.qza` `alpha-rarefaction.qzv`
```bash
#1 summarizes feature table w/ filtered data + visualize it alongside sample metadata
qiime feature-table summarize \
  --i-table spaceflight-filtered-table.qza \
  --o-visualization table-spaceflight-filtered-and-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/nasa/ulcers_metadata.tsv

#2 summarizes feature table + visualize it alongside sample metadata
qiime feature-table summarize \
  --i-table feature-frequency-filtered-table.qza \
  --o-visualization table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/nasa/ulcers_metadata.tsv
  
# generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
# alpha rarefaction -> set depth to be 125,000 
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 120000 \
  --m-metadata-file /mnt/datasets/project_2/nasa/ulcers_metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
```