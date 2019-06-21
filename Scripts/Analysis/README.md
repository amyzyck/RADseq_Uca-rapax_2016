# Population Connectivity Analyses

The R Markdown files in this directory contain the code for making files for other programs and then running various analyses. I separated the R Markdown files into All SNP Data (before filtering based on outliers), Neutral SNP Data, and Outlier SNP Data. All VCF files loaded into R stem from `SNP.TRSdp5g5mafMIap9g9dMMsnpDHWEmaf0252Amaf05.recode.vcf`.

*Forewarning: the code is very repetitive across all R Markdown files.

### Making files for other programs

For the ALL SNP dataset (before filtering based on outliers), VCF files were converted the into the following formats:

- genind object for PCA, DAPC, and popgen summary statistics
- hierfstat for calculating F_{ST}
- Inputs for EEMS

For the Neutral SNP dataset, VCF files were converted into the following formats:

- genind object for PCA, DAPC, and popgen summary statistics
- hierfstat for calculating F_{ST}
- Inputs for EEMS

For the Outlier SNP dataset, VCF files were converted into the following formats:

- genind object for PCA, DAPC, and popgen summary statistics
- hierfstat for calculating F_{ST}

### Analyses performed

- Pairwise Fst using [hierfstat](https://github.com/jgx65/hierfstat)
- Popgen summary statistics (Ho, He, overall Fst, Fis) using [hierfstat](https://github.com/jgx65/hierfstat)
- PCA using [adegenet](https://github.com/thibautjombart/adegenet) (plotting using [PCAviz](https://github.com/NovembreLab/PCAviz))
- DAPC using [adegenet](https://github.com/thibautjombart/adegenet)
- [EEMS](https://github.com/dipetkov/eems)
