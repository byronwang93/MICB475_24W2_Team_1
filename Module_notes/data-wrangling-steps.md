# Byron's scratch notes

## Important notes
- a few files such as `ulcers-demux.qza` is too large to be stored on github (2.07GB), will look into potentially storing zipped files onto this github

## Steps

1. created `ulcers` directory in `/data`
2. generating the `ulcers-demux.qza` file
```
   qiime tools import \
  --type EMPSingleEndSequences \
  --input-path /mnt/datasets/project_1/moving_pictures/emp-single-end-sequences \
  --output-path emp-single-end-sequences.qza
```
   - a few problems w/ our data
      - manifest file had columns for **paired-end** data but QIIME 2 initially expected a **single-end** manifest format
         - :. we switched to using V2 manifest format (`PairedEndFastqManifestPhred33V2`)
      - our manifest file was in Windows format (CRLF) so QIIME2 couldn't interpret header line
         - we converted the file to Unix line endings (LF) using `sed` after first copying it to a writable location
```
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
```
qiime demux summarize \
--i-data ulcers-demux.qza \
--o-visualization ulcers-demux.qzv
```
6. denoising and clustering step - `table.qza` `rep-seqs.qza` `denoising-stats.qza`
   - taking over an hour to run (decided to run in separate screen)
```
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
```
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
8. 
