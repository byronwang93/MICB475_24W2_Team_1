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
