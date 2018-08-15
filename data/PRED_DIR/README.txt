In order to execute bedtools to measure the overlap either among the CNVs 
predicted by multiple callers or between the CNV predictions and their
corresponding gold-standard CNVs for each sample and CNV type, they
must be separated into smaller files (one for each CNV type per sample).

Such smaller files generated from CNVs predictions are written to this 
directory. This subdirectory will keep the main 'data' directory
clutter-free.
