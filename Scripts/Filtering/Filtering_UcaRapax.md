# Filtering
Author: Amy Zyck

Date: June 21, 2019

Steps for filtering followed the SNP Filtering Tutorial: http://www.ddocent.com/filtering/

We will take the `TotalRawSNPs.vcf` file that dDocent created and conduct additional filtering using [VCFtools](http://vcftools.sourceforge.net/).

In `ddocent_env` directory within `Fiddler_Crab` make a new filtering directory.
```javascript
mkdir filtering
cd filtering
ln -s ../TotalRawSNPs.vcf .
```

### Scripts used for filtering steps are located in the [dDocent repository](https://github.com/jpuritz/dDocent/tree/master/scripts) on GitHub.

Change all genotypes with less than 5 reads to missing data
```javascript
vcftools --vcf TotalRawSNPs.vcf --recode-INFO-all --minDP 5 --out TRSdp5 --recode
```
```javascript
After filtering, kept 376 out of 376 Individuals
Outputting VCF file...
After filtering, kept 711393 out of a possible 711393 Sites
Run Time = 768.00 seconds
```
Filter out all variants that are not successfully genotyped in at least 50% of samples and do not have a minimum quality score of 20.
```javascript
vcftools --vcf TRSdp5.recode.vcf --recode-INFO-all --max-missing 0.5 --minQ 20 --out TRSdp5g5 --recode
```
```javascript
After filtering, kept 376 out of 376 Individuals
Outputting VCF file...
After filtering, kept 386185 out of a possible 711393 Sites
Run Time = 451.00 seconds
```

MAF filtering for FreeBayes Output
```javascript
curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/untested/multi.maf.sh
chmod +x multi.maf.sh
multi.maf.sh TRSdp5g5.recode.vcf 0.001 TRSdp5g5maf
```
```javascript
After filtering, kept 376 out of 376 Individuals
Outputting VCF file...
After filtering, kept 385420 out of a possible 386185 Sites
Run Time = 528.00 seconds
```

Use a custom script called `filter_missing_ind.sh` to filter out bad individuals
```javascript
curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/filter_missing_ind.sh
chmod +x filter_missing_ind.sh
./filter_missing_ind.sh TRSdp5g5maf.recode.vcf TRSdp5g5maf5MIa
```
```javascript
After filtering, kept 376 out of 376 Individuals
Outputting Individual Missingness
After filtering, kept 385420 out of a possible 385420 Sites
Run Time = 39.00 seconds

                                          Histogram of % missing data per individual
     Number of Occurrences
       40 ++---------+--**-----+----------+---------+----------+----------+---------+----------+---------+---------++
          +          +  **     +          +         +          'totalmissing' using (bin($1,binwidth)):(1.0) ****** +
          |             **                                                                                          |
       35 ++            **                                                                                         ++
          |             ***                                                                                         |
          |             ***                                                                                         |
       30 ++           ****                                                                                        ++
          |            ****                                                                                         |
          |            ****                                                                                         |
       25 ++           ****                                                                                        ++
          |            *****                                                                                        |
       20 ++           *****                                                                                       ++
          |         ********                                                                                        |
          |         **********                                                                                      |
       15 ++        **********                                                                                     ++
          |         **********                                                                                      |
          |         ****************                                                                                |
       10 ++        ****************                                                          **                   ++
          |       ******************                                                          **                    |
          |       * ******************                                                        **                    |
        5 ++ **** * **************** ***                                                      ***                  ++
          | ******* **************** ***********                                             ****          *******  |
          + ******* **************** *******  **************************************************************     ***+
        0 ++********************************************************************************************************+
          0         0.1       0.2        0.3       0.4        0.5        0.6       0.7        0.8       0.9         1
                                                       % of missing data

The 85% cutoff would be 0.285997
Would you like to set a different cutoff, yes or no
yes
Please eneter new cutoff.
0.5

After filtering, kept 340 out of 376 Individuals
Outputting VCF file...
After filtering, kept 385420 out of a possible 385420 Sites
Run Time = 497.00 seconds
```

Use a second custom script `pop_missing_filter.sh` to filter loci that have high missing data values in a single population. This step needs a file that maps individuals to populations `popmap`.
```javascript
ln -s ../popmap .
```
```javascript
curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/pop_missing_filter.sh
chmod +x pop_missing_filter.sh
./pop_missing_filter.sh TRSdp5g5maf5MIa.recode.vcf popmap 0.1 1 TRSdp5g5mafMIap9
```
```javascript
After filtering, kept 340 out of 340 Individuals
Outputting VCF file...
After filtering, kept 149880 out of a possible 385420 Sites
Run Time = 274.00 seconds
```
Filter out any sites with less than 95% overall call rate and MAF of 0.001.
```javascript
vcftools --vcf TRSdp5g5mafMIap9.recode.vcf --recode-INFO-all --max-missing 0.95 --maf 0.001 --out TRSdp5g5mafMIap9g9 --recode
```
```javascript
After filtering, kept 340 out of 340 Individuals
Outputting VCF file...
After filtering, kept 120155 out of a possible 149880 Sites
Run Time = 259.00 seconds
```
Use another custom filter script `dDocent_filters`
```javascript
curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/dDocent_filters
chmod +x dDocent_filters
./dDocent_filters TRSdp5g5mafMIap9g9.recode.vcf TRSdp5g5mafMIap9g9d
```
```javascript
This script will automatically filter a FreeBayes generated VCF file using criteria related to site depth,
quality versus depth, strand representation, allelic balance at heterzygous individuals, and paired read representation.
The script assumes that loci and individuals with low call rates (or depth) have already been removed.

Contact Jon Puritz (jpuritz@gmail.com) for questions and see script comments for more details on particular filters

Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth
 76959 of 120155

Are reads expected to overlap?  In other words, is fragment size less than 2X the read length?  Enter yes or no.
yes
Is this from a mixture of SE and PE libraries? Enter yes or no.
no
Number of additional sites filtered based on properly paired status
 1227 of 43196

Number of sites filtered based on high depth and lower than 2*DEPTH quality score
 4189 of 41969

                                               Histogram of mean depth per site

      450 +---------------------------------------------------------------------------------------------------------+
          | +    +     +    +    +     +    +     +    +    +     +   *+ *   +    +     +    +    +     +    +     +|
          |                                               **eandepthpe*si*e' using (bin($1,binwidth)):(1.0) ******* |
      400 |-+                                             ***        *** *                                        +-|
          |                                            ** *** *     **** ** *                                       |
      350 |-+                                   **    *** *** *   ********* *   *                                 +-|
          |                                     **    ******* * *************   *                                   |
          |                                     **    ******* ***************   *                                   |
      300 |-+                                   **   ******** ****************  *                                 +-|
          |                                     **   *************************  *                                   |
      250 |-+                                   **   *************************  ***                               +-|
          |                                     ***********************************                                 |
          |                        *  ** **  * ************************************                                 |
      200 |-+                      * *** **  ****************************************                             +-|
          |                        * ****** ******************************************                              |
      150 |-+                      * **************************************************  ***                      +-|
          |                      * ****************************************************  ***                        |
          |                  *  ************************************************************                        |
      100 |-+                ***************************************************************                      +-|
          |              ** *******************************************************************                     |
       50 |-+            **********************************************************************  **   **          +-|
          |           ***************************************************************************** ****       **   |
          | +  *****************************************************************************************************|
        0 +---------------------------------------------------------------------------------------------------------+
            15   30    45   60   75    90  105   120  135  150   165  180   195  210   225  240  255   270  285   300
                                                          Mean Depth

If distrubtion looks normal, a 1.645 sigma cutoff (~90% of the data) would be 176038.7185
The 95% cutoff would be 277
Would you like to use a different maximum mean depth cutoff than 277, yes or no
no
Number of sites filtered based on maximum mean depth
 2176 of 41969

Number of sites filtered based on within locus depth mismatch
 35 of 39793

Total number of sites filtered
 80397 of 120155

Remaining sites
 39758

Filtered VCF file is called Output_prefix.FIL.recode.vcf

Filter stats stored in TRSdp5g5mafMIap9g9d.filterstats
```
Break complex mutational events (combinations of SNPs and INDELs) into sepearte SNP and INDEL calls, and then remove INDELs.
```javascript
vcfallelicprimitives TRSdp5g5mafMIap9g9d.FIL.recode.vcf --keep-info --keep-geno > TRSdp5g5mafMIap9g9d.prim.vcf
vcftools --vcf TRSdp5g5mafMIap9g9d.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.TRSdp5g5mafMIap9g9d
```
```javascript
After filtering, kept 340 out of 340 Individuals
Outputting VCF file...
After filtering, kept 43983 out of a possible 49258 Sites
Run Time = 93.00 seconds
```

This data set contains technical replicates. I will use a custom script `dup_sample_filter.sh` to automatically remove sites in VCF files that do not have congruent genotypes across duplicate individuals. It will only consider genotypes that have at least 5 reads.

`dup_sample_filter.sh` is located on the dDocent GitHub and also here: https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/blob/master/Scripts/Filtering/dup_sample_filter.sh

The technical replicates are listed in `duplicates.samples.1` with the following format (each name should be separated by a tab):
```javascript
FBN_327a        FBN_327b
FBN_327a        FBN_327c
FBN_327a        FBN_327d
FBN_327b        FBN_327c
FBN_327b        FBN_327d
FBN_327c        FBN_327d
OBN_9b         OBN_9c
OBN_9b         OBN_9d
OBN_9c         OBN_9d
OBS_245a        OBS_245b
OBS_245a        OBS_245c
OBS_245a        OBS_245d
OBS_245b        OBS_245c
OBS_245b        OBS_245d
OBS_245c        OBS_245d
PCN_210a        PCN_210b
PCN_223a        PCN_223b
PCS_361a        PCS_361b
PCS_365a        PCS_365b
SPS_74a        SPS_74b
SPS_92a        SPS_92b
SPS_92a        SPS_92c
SPS_92a        SPS_92d
SPS_92b        SPS_92c
SPS_92b        SPS_92d
SPS_92c        SPS_92d
WC2_301a        WC2_301b
WC2_305a        WC2_305b
WC2_305a        WC2_305c
WC2_305b        WC2_305c
```  
Copy to `filtering`.
```javascript
ln -s ../dup_sample_filter.sh .
ln -s ../duplicates.samples.1 .
```
Run script.
```javascript
./dup_sample_filter.sh SNP.TRSdp5g5mafMIap9g9d.recode.vcf duplicates.samples.1
```
This produces a `mismatched.loci` file.
```javascript
head mismatched.loci
```
```javascript
2       dDocent_Contig_6167     14
2       dDocent_Contig_12111    122
2       dDocent_Contig_5485     184
2       dDocent_Contig_11981    240
14      dDocent_Contig_5512     27
2       dDocent_Contig_11857    114
2       dDocent_Contig_11642    216
2       dDocent_Contig_3841     25
2       dDocent_Contig_11322    100
2       dDocent_Contig_3530     210
```

```javascript
echo -e "Mismatches\tNumber_of_Loci" > mismatch.txt
for i in {2..20}
do
paste <(echo $i) <(mawk -v x=$i '$1 > x' mismatched.loci | wc -l) >> mismatch.txt
done
```
In RStudio
```javascript
library(ggplot2)
mismatch <- read.table("mismatch.txt", header = TRUE)
df=data.frame(mismatch)

p <- ggplot(df, aes(x=Mismatches, y=Number_of_Loci)) + geom_point() +theme_bw() + scale_x_continuous(minor_breaks = seq(1,20,by=1))
p
```
![mismatched](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/blob/master/Output/Filtering/mismatched.png)

I decided to filter out loci with mismatched values > 6.
```javascript
mawk '$1 > 6' mismatched.loci | cut -f2,3 > mismatchedloci

vcftools --vcf SNP.TRSdp5g5mafMIap9g9d.recode.vcf --exclude-positions mismatchedloci --recode --recode-INFO-all --out SNP.TRSdp5g5mafMIap9g9dMM
```
```javascript
After filtering, kept 340 out of 340 Individuals
Outputting VCF file...
After filtering, kept 43532 out of a possible 43983 Sites
Run Time = 97.00 seconds
```


Use the script `rad_haplotyper.pl` written by [Chris Hollenbeck](https://github.com/chollenbeck). This tool takes a VCF file of SNPs and will parse through BAM files looking to link SNPs into haplotypes along paired reads.

First, copy the most recent VCF file to the directory with the BAM files and `reference.fasta`. In my case, it is `ddocent_env`.
```javascript
cp SNP.TRSdp5g5mafMIap9g9dMM.recode.vcf ../
cd ../
```

```javascript
curl -L -O https://raw.githubusercontent.com/chollenbeck/rad_haplotyper/e8bdc79f1d1d9ce3d769996315fc1ffd3a7d0e4e/rad_haplotyper.pl
chmod +x rad_haplotyper.pl
rad_haplotyper.pl -p popmap -v SNP.TRSdp5g5mafMIap9g9dMM.recode.vcf -r reference.fasta -g SNPTRSdp5g5mafMIap9g9dMM -mp 5 -x 40 -z 0.1 -e
```
Output will resemble
```javascript
Building haplotypes for FBN_306
Building haplotypes for FBN_307
Building haplotypes for FBN_308
Building haplotypes for FBN_310
Building haplotypes for FBN_311
...
Filtered 182 loci below missing data cutoff
Filtered 187 possible paralogs
Filtered 0 loci with low coverage or genotyping errors
Filtered 2 loci with an excess of haplotypes
Writing Genepop file: SNPTRSdp5g5mafMIap9g9dMM
```
All additional loci to be removed are stored in `stats.out`
```javascript
head stats.out
```
```javascript
Locus                   Sites   Haplotypes      Inds_Haplotyped Total_Inds      Prop_HaplotypedStatus   Poss_Paralog    Low_Cov/Geno_Err        Miss_Geno       Comment
dDocent_Contig_10004    19      23              339             340             0.997                   PASSED          0                        01
dDocent_Contig_10012    1       2               323             340             0.950                   PASSED          0                        017
dDocent_Contig_10017    18      21              340             340             1.000                   PASSED          0                        00
dDocent_Contig_10019    1       2               340             340             1.000                   PASSED          0                        00
dDocent_Contig_10030    32      28              334             340             0.982                   PASSED          1                        32
dDocent_Contig_10039    13      15              338             340             0.994                   PASSED          0                        02
dDocent_Contig_10054    32      64              338             340             0.994                   PASSED          0                        20
dDocent_Contig_10057    10      10              332             340             0.976                   PASSED          0                        17                       
```

Use this file to create a list of loci to filter.
```javascript
mawk '/FILT/' stats.out | cut -f1 > bad.hap.dp3.loci
```

Remove the bad RAD loci using the script `remove.bad.hap.loci.sh`
```javascript
curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/remove.bad.hap.loci.sh
chmod +x remove.bad.hap.loci.sh
./remove.bad.hap.loci.sh bad.hap.dp3.loci SNP.TRSdp5g5mafMIap9g9dMM.recode.vcf
```
This generates the filtered VCF file `SNP.TRSdp5g5mafMIap9g9dMM.filtered.vcf`

```javascript
vcftools --vcf SNP.TRSdp5g5mafMIap9g9dMM.filtered.vcf --missing-indv
```
```javascript
After filtering, kept 340 out of 340 Individuals
Outputting Individual Missingness
After filtering, kept 38121 out of a possible 38121 Sites
Run Time = 2.00 seconds
```

```javascript
head out.imiss
```
```javascript
INDV    N_DATA  N_GENOTYPES_FILTERED    N_MISS  F_MISS
FBN_306 38121   0       405     0.0106241
FBN_307 38121   0       158     0.0041447
FBN_308 38121   0       164     0.00430209
FBN_309 38121   0       683     0.0179166
FBN_310 38121   0       117     0.00306917
FBN_311 38121   0       73      0.00191496
FBN_312 38121   0       91      0.00238714
FBN_313 38121   0       117     0.00306917
FBN_314 38121   0       246     0.00645314
```

```javascript
mawk '/FBN_327/ || /OBN_9/ || /OBN_009/ || /OBS_245/ || /PCN_210/ || /PCN_223/ || /PCS_361/ || /PCS_365/ || /SPC_92/ || /WC2_301/ || /WC2_305/ || /SPS_74/ || /SPS_92/' out.imiss > dup.imiss
cat dup.imiss
```
```javascript
FBN_327a        38121   0       17      0.000445948
FBN_327b        38121   0       7       0.000183626
FBN_327c        38121   0       10      0.000262323
FBN_327d        38121   0       43      0.00112799
OBN_009a        38121   0       94      0.00246583
OBN_9b          38121   0       3       7.86968e-05
OBN_9c          38121   0       10      0.000262323
OBN_9d          38121   0       51      0.00133785
OBS_245a        38121   0       33      0.000865665
OBS_245b        38121   0       6       0.000157394
OBS_245c        38121   0       212     0.00556124
OBS_245d        38121   0       11      0.000288555
PCN_210a        38121   0       92      0.00241337
PCN_210b        38121   0       77      0.00201988
PCN_223a        38121   0       146     0.00382991
PCN_223b        38121   0       112     0.00293801
PCS_361a        38121   0       233     0.00611212
PCS_361b        38121   0       110     0.00288555
PCS_365a        38121   0       133     0.00348889
PCS_365b        38121   0       183     0.0048005
SPS_74a         38121   0       37      0.000970594
SPS_74b         38121   0       69      0.00181003
SPS_92a         38121   0       6       0.000157394
SPS_92b         38121   0       34      0.000891897
SPS_92c         38121   0       1       2.62323e-05
SPS_92d         38121   0       84      0.00220351
WC2_301a        38121   0       96      0.0025183
WC2_301b        38121   0       486     0.0127489
WC2_305a        38121   0       350     0.00918129
WC2_305b        38121   0       88      0.00230844
WC2_305c        38121   0       269     0.00705648
```

```javascript
mawk '!/FBN_327b/ && !/OBN_9b/ && !/OBS_245b/ && !/PCN_210b/ && !/PCN_223b/ && !/PCS_361b/ && !/PCS_365a/ && !/SPS_74a/ && !/SPS_92c/ && !/WC2_301a/ && !/WC2_305b/' dup.imiss > duplicate.samples.to.remove
vcftools --vcf SNP.TRSdp5g5mafMIap9g9dMM.filtered.vcf --recode --recode-INFO-all --out SNP.TRSdp5g5mafMIap9g9dMM --remove duplicate.samples.to.remove
```
```javascript
After filtering, kept 320 out of 340 Individuals
Outputting VCF file...
After filtering, kept 38121 out of a possible 38121 Sites
Run Time = 36.00 seconds
```

Filter out loci that are out of HWE in more than half the populations, using `filter_hwe_by_pop.pl` written by [Chris Hollenbeck](https://github.com/chollenbeck)
```javascript
curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/filter_hwe_by_pop.pl
chmod +x filter_hwe_by_pop.pl
./filter_hwe_by_pop.pl -v SNP.TRSdp5g5mafMIap9g9dMM.recode.vcf -c 0.5 -p popmap -o SNP.TRSdp5g5mafMIap9g9dMMHWE
```
```javascript
Processing population: FBN (33 inds)
Processing population: FBS (30 inds)
Processing population: OBN (31 inds)
Processing population: OBS (32 inds)
Processing population: PCN (33 inds)
Processing population: PCS (32 inds)
Processing population: SPN (32 inds)
Processing population: SPS (33 inds)
Processing population: WC1 (28 inds)
Processing population: WC2 (32 inds)
Processing population: WC3 (29 inds)
Processing population: WC4 (31 inds)
Outputting results of HWE test for filtered loci to 'filtered.hwe'
Kept 38003 of a possible 38121 loci (filtered 118 loci)
```
```javascript
vcftools --vcf SNP.TRSdp5g5mafMIap9g9dMMHWE.recode.vcf --recode --recode-INFO-all --out SNP.TRSdp5g5mafMIap9g9dMMHWEmaf0252A --maf 0.025 --max-alleles 2
```
```javascript
After filtering, kept 320 out of 320 Individuals
Outputting VCF file...
After filtering, kept 10597 out of a possible 38003 Sites
Run Time = 11.00 seconds
```

****
**I repeated the filtering, but omitted the `rad_haplotyper` steps, to generate a VCF file of SNPs: `SNP.TRSdp5g5mafMIap9g9dMMsnpDHWEmaf0252Amaf05.recode.vcf`.**
