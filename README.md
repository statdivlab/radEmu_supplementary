## Code to reproduce the results of 

This repository contains code to reproduce the results and figures of our manuscript "Estimating Fold Changes from Partially Observed Outcomes with Applications in Microbial Metagenomics" by David Clausen and Amy Willis. 

Code to create each figure is contained within the relevant labeled subfolder. 

Is something missing? Let us know by emailing Amy or opening an issue. 

We are grateful to Wirbel et al (the authors of the colorectal cancer metaanalysis) for making their data publicly available and easy to download. Code to download their data is available in `fig456/0-process.R` and reproduced below for your convenience:

```
### mOTU table
download.file(url = "https://zenodo.org/records/3517209/files/species_profiles_g2_l75_motus2.0.0.tsv?download=1", destfile = "species_profiles_g2_l75_motus2.0.0.tsv", method = "libcurl")

### covariate information
download.file(url = "https://zenodo.org/records/3517209/files/meta_all.tsv?download=1", destfile = "meta_all.tsv")


### Wirbel et al's p-values from a blocked Wilcoxon test on proportion-scale mOTU data
download.file(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-019-0406-6/MediaObjects/41591_2019_406_MOESM3_ESM.xlsx", destfile = "orig_wirbel_results.xlsx")
```

