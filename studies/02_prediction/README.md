# Prediction pilot skeleton

Planned pilot: three traits; mostly-shared and mostly-trait-specific architectures; 10 replicates; a common 70/30 train/test split; ST-BED BayesR, MT-BED BayesR, ST-CSR SBayesR, and MT-CSR SBayesR. Initial metrics are correlation and MSE against true genetic value. Future metrics are phenotype correlation, calibration slope, and MT-minus-ST prediction accuracy.

Every method in a replicate must use the same split, and test phenotypes must not be used for fitting. Oracle genotype-scale validation is mandatory before scoring. Marker, sample, and trait alignment are strict. This study may not make implementation changes in `sblr`.
