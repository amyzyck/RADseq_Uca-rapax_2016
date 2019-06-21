## Outlier Detection

[PCAdapt](https://github.com/bcm-uga/pcadapt), [OutFLANK](https://github.com/whitlock/OutFLANK), [BayeScan](https://github.com/mfoll/BayeScan), and [Bayenv2](https://gcbias.org/bayenv/) were used to identify any outliers. Outlier detection programs were run on `SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252A.recode.vcf` with haplotypes. Those results are presented below.

Steps were repeated on `SNP.TRSdp5g5mafMIap9g9dMMsnpDHWEmaf0252Amaf05.recode.vcf` with SNPs. Results are reported on at the end of this document.

### PCAdapt
In RStudio
```javascript
#Load pcadapt library
library(pcadapt)

vcf2pcadapt("SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252A.recode.vcf", output = "tmp.pcadapt", allele.sep = c("/", "|"))
```
```javascript
No variant got discarded.
Summary:

	- input file:				SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252A.recode.vcf
	- output file:				tmp.pcadapt

	- number of individuals detected:	320
	- number of loci detected:		10597
```
```javscript
#Make the pcadapt file "readable"
filename <- read.pcadapt("tmp.pcadapt", type = "pcadapt")
```
```javascript
10597 lines detected.
320 columns detected.
```

```javascript
#Create first PCA
x <- pcadapt(input = filename, K = 20)
```
```javascript
#Plot the likelihoods
plot(x, option = "screeplot")
```
![K20](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/Output/Outlier_Detection/ScreePlotK20.png)

2 or 3 seem to be a good cut off. Now zoomed in:
```javascript
plot(x, option = "screeplot", K = 10)
```
![K10](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Output/Outlier_Detection/ScreePlotK10.png)

I'm going to pick 3.
```javascript
x <- pcadapt(input = filename, K = 3)
```

#### Calculate population designations
In terminal
```javascript
vcftools --vcf SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252A.recode.vcf --missing-indv
```
```javascript
After filtering, kept 320 out of 320 Individuals
Outputting Individual Missingness
After filtering, kept 10597 out of a possible 10597 Sites
Run Time = 1.00 seconds
```
```javascript
cut -f1 out.imiss |grep -v INDV| cut -f1 -d "_" | sort | uniq -c
```
```javascript
     30 FBN
     25 FBS
     19 OBN
     28 OBS
     27 PCN
     25 PCS
     28 SPN
     29 SPS
     25 WC1
     27 WC2
     28 WC3
     29 WC4
```

Back in RStudio
```javascript
# create population designations
poplist.names <- c(rep("FBN", 30),rep("FBS", 25),rep("OBN", 19), rep("OBS",28), rep("PCN",27), rep("PCS",25), rep("SPN",28), rep("SPS",29), rep("WC1",25), rep("WC2",27), rep("WC3",28), rep("WC4",29))
```

```javascript
#Plot the actual PCA (first two PCAs)
plot(x, option = "scores", pop = poplist.names)
```

![PCA](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/Output/Outlier_Detection/PCA1%262.png)

Very little spatial structure :(

### Looking for Outliers
#### Manhattan Plot
```javascript
#Start looking for outliers
#Make Manhattan Plot
plot(x , option = "manhattan")
```
![Manhattan](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Output/Outlier_Detection/ManhattanPlot.png)

#### Q-Q Plot
```javascript
#Make qqplot
plot(x, option = "qqplot", threshold = 0.05)
```

![QQPlot](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Output/Outlier_Detection/QQplot.png)

#### Look at P-value distribution
```javascript
# Look at P-value distribution
plot(x, option = "stat.distribution")
```

![pvalue](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Output/Outlier_Detection/pvalue.png)

```javascript
library(qvalue)
qval1 <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
```

#### Save outliers
```javascript
outliers1 <- which(qval1 < alpha)
outliers1
```
```javascript
  [1]   51   78  113  164  246  247  553  667  701  933 1481 1850 2194 2256 2540
 [16] 2683 2737 2858 2859 2860 3665 3783 3785 3786 3854 3855 3856 3857 3909 3910
 [31] 3912 3913 3914 3916 3917 3918 3919 3920 3921 3923 3924 3925 3926 4274 4275
 [46] 4277 4699 4700 5348 5401 5406 5407 5409 5411 5413 5414 5415 5416 5417 5418
 [61] 5774 5861 5862 5863 5865 5866 5867 5868 5869 5870 5871 5872 5873 5874 5875
 [76] 5876 5877 5878 5880 5881 5882 5884 6009 6010 6076 6077 6078 6079 6080 6081
 [91] 6082 6083 6084 6085 6086 6087 6088 6089 6092 6093 6873 6875 6878 6991 7744
[106] 7997 8606 8616 8619 8620 8621 8622 8623 8624 8634 9200 9201 9203 9204 9208
[121] 9210 9214
```

```javascript
system("rm outliers1.txt", wait=FALSE)
invisible(lapply(outliers1, write, "outliers1.txt", append=TRUE))
```

In terminal
```javascript
head outliers1.txt
```
```javascript
51
78
113
164
246
247
553
667
701
933
```

```javascript
mawk '!/#/' SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252A.recode.vcf | cut -f1,2 > totalloci
NUM=(`cat totalloci | wc -l`)
paste <(seq 1 $NUM) totalloci > loci.plus.index
cat outliers1.txt | parallel "grep -w ^{} loci.plus.index" | cut -f2,3> outlier.loci.txt

head outlier.loci.txt
```
```javascript
dDocent_Contig_735      23
dDocent_Contig_820      43
dDocent_Contig_1033     282
dDocent_Contig_1242     300
dDocent_Contig_1597     285
dDocent_Contig_1597     286
dDocent_Contig_2344     177
dDocent_Contig_2522     164
dDocent_Contig_2565     10
dDocent_Contig_2813     77
```

### Outflank
Work with SNPs with MAF > 0.05

In terminal
```javascript
vcftools --vcf SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252A.recode.vcf --maf 0.05 --recode --recode-INFO-all --out SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05
```
```javascript
After filtering, kept 320 out of 320 Individuals
Outputting VCF file...
After filtering, kept 6786 out of a possible 10597 Sites
Run Time = 6.00 seconds
```

In RStudio
```javascript
library(OutFLANK)  # outflank package
library(vcfR)
library(bigsnpr)   # package for LD pruning
```

```javascript
my_vcf <- read.vcfR("SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05.recode.vcf")
```
```javascript
Scanning file to determine attributes.
File attributes:
  meta lines: 64
  header_line: 65
  variant count: 6786
  column count: 329
Meta line 64 read in.
All meta lines processed.
gt matrix initialized.
Character matrix gt created.
  Character matrix gt rows: 6786
  Character matrix gt cols: 329
  skip: 0
  nrows: 6786
  row_num: 0
Processed variant: 6786
All variants processed
```

```javascript
geno <- extract.gt(my_vcf) # Character matrix containing the genotypes
position <- getPOS(my_vcf) # Positions in bp
chromosome <- getCHROM(my_vcf) # Chromosome information

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

G[is.na(G)] <- 9

head(G[,1:10])
```
```javascript
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
[1,]    0    0    1    0    0    0    0    0    0     0
[2,]    1    0    0    0    0    0    0    0    0     0
[3,]    0    0    1    0    0    0    0    0    0     0
[4,]    1    1    1    0    0    1    0    0    0     0
[5,]    1    0    1    0    0    1    0    0    0     0
[6,]    1    0    1    0    0    0    0    0    0     0
```

```javascript
pop <- as.vector(poplist.names)
```

```javascript
my_fst <- MakeDiploidFSTMat(t(G), locusNames = paste0(chromosome,"_", position), popNames = pop)
```

```javascript
plot(my_fst$He, my_fst$FST)
```

![FST](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Output/Outlier_Detection/FST1.png)

```javascript
plot(my_fst$FST, my_fst$FSTNoCorr)
abline(0,1)
```

![NoCorr](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Output/Outlier_Detection/NoCorr.png)

We need to give OUTFlank a set of quasi-independent SNPs to estimate the neutral FST distribution. To approximate this, we will prune our SNPs to one per RAD contig

In terminal
```javascript
curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/untested/Filter_one_random_snp_per_contig.sh
chmod +x Filter_one_random_snp_per_contig.sh
./Filter_one_random_snp_per_contig.sh SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05.recode.vcf
```
```javascript
Filtered VCF file is saved under name SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05.filtered1SNPper.vcf
```

In RStudio
```javascript
my_vcf_sub <- read.vcfR("SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05.filtered1SNPper.vcf")
```
```javascript
Scanning file to determine attributes.
File attributes:
  meta lines: 64
  header_line: 65
  variant count: 2426
  column count: 329
Meta line 64 read in.
All meta lines processed.
gt matrix initialized.
Character matrix gt created.
  Character matrix gt rows: 2426
  Character matrix gt cols: 329
  skip: 0
  nrows: 2426
  row_num: 0
Processed variant: 2426
All variants processed
```
```javascript
geno_sub <- extract.gt(my_vcf_sub) # Character matrix containing the genotypes
position_sub <- getPOS(my_vcf_sub) # Positions in bp
chromosome_sub <- getCHROM(my_vcf_sub) # Chromosome information

G_sub <- matrix(NA, nrow = nrow(geno_sub), ncol = ncol(geno_sub))

G_sub[geno_sub %in% c("0/0", "0|0")] <- 0
G_sub[geno_sub  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_sub[geno_sub %in% c("1/1", "1|1")] <- 2

G_sub[is.na(G_sub)] <- 9

head(G_sub[,1:10])
```
```javascript
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
[1,]    0    0    1    0    0    0    0    0    0     0
[2,]    0    0    1    0    0    0    0    0    0     0
[3,]    1    0    1    0    0    0    0    0    0     0
[4,]    1    0    0    0    1    0    0    0    0     1
[5,]    1    2    2    2    1    1    2    2    2     0
[6,]    0    0    1    0    0    0    0    0    0     0
```

```javascript
pop <- as.vector(poplist.names)
```

```javascript
my_fst_sub <- MakeDiploidFSTMat(t(G_sub), locusNames = paste0(chromosome_sub,"_", position_sub), popNames = pop)
```

```javascript
plot(my_fst_sub$He, my_fst_sub$FST)
```

![FST2](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Output/Outlier_Detection/FST2.png)

```javascript
out_trim <- OutFLANK(my_fst_sub, NumberOfSamples = 12, qthreshold=0.05, RightTrimFraction=0.05, LeftTrimFraction=0.05,Hmin =0.05)
```
```javascript
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin =0.05, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
```

![NoCorr](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Output/Outlier_Detection/Nosamplesizecorr.png)

```javascript
hist(out_trim$results$pvaluesRightTail)
```

![Histo](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Output/Outlier_Detection/Histo.png)

```javascript
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar,
                                   dfInferred = out_trim$dfInferred, qthreshold = 0.1, Hmin =0.05)
```

```javascript
my_out <- P1$OutlierFlag==TRUE
plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1))
points(P1$He[my_out], P1$FST[my_out], col="blue")
```

![Outlier](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Output/Outlier_Detection/Outlierpoint.png)

```javascript
P1[which(P1$OutlierFlag == TRUE),]
```

```javascript
        LocusName                         He           FST             T1
3727	dDocent_Contig_9324_207	0.3041831	0.05602335	0.008578829
```

### BayeScan

#### Convert VCF file to other outputs

In terminal
```javascript
cp /home/BIO594/DATA/Week7/example/BSsnp.spid

java -jar /usr/local/bin/PGDSpider2-cli.jar -spid BSsnp.spid -inputfile SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05.recode.vcf -outputfile SNP1.BS
```

```javascript
BayeScan2.1_linux64bits SNP1.BS -thin 50 -nbp 30
```
Copy R source file
```javascript
cp /home/azyck/Week7and8/simulated/plot_R.r
```

In RStudio
```javascript
source("plot_R.r")
plot_bayescan("SNP_fst.txt")
```

![Baye1](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Output/Outlier_Detection/Baye1.png)

```javascript
$outliers
[1] 3727

$nb_outliers
[1] 1
```

```javascript
bs <- plot_bayescan("SNP_fst.txt", FDR = 0.1)
```

![Baye2](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Output/Outlier_Detection/Baye2.png)

```javascript
bs$outliers
```
```javascript
[1] 3727
```

```javascript
mawk '!/#/' SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05.recode.vcf | cut -f1,2 > totalloci
NUM=(`cat totalloci | wc -l`)
paste <(seq 1 $NUM) totalloci > loci.plus.index
echo -e "3727" | parallel "grep -w ^{} loci.plus.index" | cut -f2,3> outlier2.loci.txt

head outlier2.loci.txt
```
```javascript
dDocent_Contig_9324     207
```

### Bayenv2
First, convert vcf to BayEnv input
```javascript
java -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05.recode.vcf -outputfile SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05BayEnv.txt -spid SNPBayEnv.spid
```

```javascript
bayenv2 -i SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05BayEnv.txt -p 12 -k 100000 -r 63479 > matrix.out
```

This code generates 100,000 iterations. We only need the last one.
```javascript
tail -13 matrix.out | head -12 > matrix
```

With the matrix we will use our environmental factor file:
```javascript
cat environ
```
```javascript
0.851194237     0.739018209     1.379424387     1.19540284      0.185880782     -0.086000325    -0.862152712    -2.284603823    -0.476198768    -0.331092896    -0.219243369    -0.091628561
0.006091651     -0.090159786    -1.693127051    -1.588203181    1.387471969     1.387471969     0.275339273     -1.069038664    0.189042281     0.294217526     0.396376273     0.504517741
```
The environmental file are standardized environmental data with each line representing an environemtal factor with the value for each population tab delimited.

Next, we calculate the Bayes Factor for each SNP for each environmental variable:
```javascript
calc_bf.sh SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05BayEnv.txt environ matrix 12 10000 2
```

Next, we convert the output into something suitable to input into R
```javascript
paste <(seq 1 981) <(cut -f2,3 bf_environ.environ ) > bayenv.out
cat <(echo -e "Locus\tBF1\tBF2") bayenv.out > bayenv.final
```

In RStudio
```javascript
table_bay <- read.table("bayenv.final",header=TRUE)
plot(table_bay$BF1)

table_bay[which(table_bay$BF1 > 100),]
```

***No Outliers detected***

#### Combine all outlier loci into one file
```javascript
cat outlier*.loci.txt > all.outliers
cut -f1 all.outliers | sort | uniq | wc -l
```
```javascript
38
```

### Create a VCF file with just neutral loci
```javascript
vcftools --vcf SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05.recode.vcf --exclude-positions all.outliers --recode-INFO-all --out neutralloci --recode
```
```javascript
After filtering, kept 320 out of 320 Individuals
Outputting VCF file...
After filtering, kept 6663 out of a possible 6786 Sites
Run Time = 7.00 seconds
```

#### 6757 loci in `neutralloci.recode.vcf`

Create outlier only vcf file
```javascript
vcftools --vcf SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05.recode.vcf --recode --recode-INFO-all --positions all.outliers --out outlierloci
```

#### 123 loci in `outlierloci.recode.vcf`

****

Results for the SNP VCF file `SNP.TRSdp5g5mafMIap9g9dMMsnpDHWEmaf0252Amaf05.recode.vcf`:

*118 outliers were detected from PCAdapt and Outflank (117 and 1, respectively).
*8,117 Neutral loci located in `neutralsnploci.recode.vcf`.
*118 outlier loci are located in `outliersnps.recode.vcf`.

****
#### Filtering one random SNP per contig

Principal Component Analysis needs 1 SNP per contig.

```javascript
curl -L -O https://github.com/jpuritz/dDocent/blob/master/scripts/untested/Filter_one_random_snp_per_contig.sh
chmod +x Filter_one_random_snp_per_contig.sh
./Filter_one_random_snp_per_contig.sh SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05.recode.vcf
```
```javasript
Filtered VCF file is saved under name SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252Amaf05.filtered1SNPper.vcf
```
Repeating on SNP VCF file
```javascript
./Filter_one_random_snp_per_contig.sh SNP.TRSdp5g5mafMIap9g9dMMsnpDHWEmaf0252Amaf05.recode.vcf
```
```javasript
Filtered VCF file is saved under name SNP.TRSdp5g5mafMIap9g9dMMsnpDHWEmaf0252Amaf05.filtered1SNPper.vcf
```
