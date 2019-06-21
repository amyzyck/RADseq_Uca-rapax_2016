## Read Trimming, Assembly, Mapping, and SNP Calling
----

Last Updated: June 21, 2019

Data uploaded and analyzed on KITT (made by J. Puritz)

Data is located: ```PATH /home/azyck/Fiddler_Crab```


## Counting raw reads
```javascript
for fq in *.fq.gz
> do
> echo $fq
> zcat $fq | echo $((`wc -l`/4))
> done
```

#### [Output](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/blob/master/Scripts/Assembly/RawReads_Counts)
****
## dDocent will be used for Quality Filtering, _De Novo_ Assembly, Read Mapping, and SNP Calling

dDocent: http://www.ddocent.com/

## **Quality control of raw sequencing reads (FASTQC)**

FASTQC: https://github.com/s-andrews/FastQC

Before running dDocent, I would like to check the quality of the raw reads

```javascript
mkdir fastqc_results
cd fastqc_results
conda install -c bioconda fastqc

ln -s ../*.fq.gz .

#fastqc takes awhile to run on larger data sets, so I recommend running it overnight
fastqc *

#fastqc generates a html and zip file for all files
#a full report is created using multiqc

multiqc .
```
This generates an .html file of the MultiQC report. To view, I copied the file from KITT to my computer using [FileZilla](https://filezilla-project.org/).

### Multiqc report summary of raw reads

#### Sequence Quality
![SeqQual](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/blob/master/Output/raw_reads/raw_MeanQualityScores.png)

This graph shows the mean quality at each position across the read. The quality looks pretty good. Quality is lower at the front and end of the sequence, but still within the green.

#### Per Sequence Quality Scores

![PerSeq](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/blob/master/Output/raw_reads/raw_PerSequenceQualityScores.png)

This graph looks at the average quality scores per sequence. Again looks pretty good. Any low quality scores will be trimmed.

#### Per Sequence GC Content

![GC](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/blob/master/Output/raw_reads/raw_GCcontent.png)

If all individuals are from the same species, then this graph will show that the sequences are roughly normally distributed. It is possible that there are individuals from another species mixed in. These individuals (if any) will be identified using Principal Components Analysis (PCA).

## Now for dDocent!

As this is my first time running dDocent on my own data, I went through the [Quick Start Guide](http://www.ddocent.com/quick/). I recommend:

1. Reading through the [User Guide](http://www.ddocent.com/UserGuide/).
2. Completing the [Assembly Tutorial](http://www.ddocent.com/assembly/), using the simulated dataset.


#### First, create a dDocent conda environment

*Back in ```Fiddler_Crab``` directory

```javascript
mkdir ddocent_env
cd ddocent_env

ln -s ../*.fq.gz .

conda create -n ddocent_env ddocent
```

#### Activate the dDocent environment
```javascript
source activate ddocent_env
```

#### Run dDocent once with read trimming to test install
```javascript
dDocent
```

```javascript
dDocent 2.7.7

Contact jpuritz@gmail.com with any problems


Checking for required software

All required software is installed!

dDocent run started Wed Apr 24 19:47:38 EDT 2019

376 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes

dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
20

dDocent detects 503 gigabytes maximum memory available on this system.
Please enter the maximum memory to use for this analysis in gigabytes.
For example, to limit dDocent to ten gigabytes, enter 10.
This option does not work with all distributions of Linux. If runs are hanging at variant calling, enter 0
Then press [ENTER]
0

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
yes

Do you want to perform an assembly?
Type yes or no and press [ENTER].
no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
no

Mapping will not be performed

Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.
```

#### Before continuing with dDocent, I am checking the quality of the trimmed reads with FastQC

*Back in ```Fiddler_Crab``` directory
```javascript
mkdir cleaned_reads
cd cleaned_reads

ln -s ../ddocent_env/*.fq.gz .

conda install -c bioconda fastqc

#fastqc takes awhile to run on larger data sets, so I recommend running it overnight
fastqc *.R1.fq.gz
fastqc *.R2.fq.gz

#a full report is created using multiqc
multiqc .
```

This generates an .html file of the MultiQC report. To view, I copied the file from KITT to my computer using [FileZilla](https://filezilla-project.org/).

### Multiqc report summary of cleaned reads

#### Sequence Quality

![cleanSeqQual](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/blob/master/Output/cleaned_reads/cleaned_MeanQualityScores.png)

#### Per Sequence Quality scores

![cleanPerSeq](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/blob/master/Output/cleaned_reads/cleaned_PerSequenceQualityScores.png)

First two are slightly better, but not much improvement was needed to begin with.

#### Per Sequence GC content

![cleanGC](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/blob/master/Output/cleaned_reads/cleaned_GCcontent.png)

GC content is still bad, but hopefully I can identify any individuals from other species with PCA.
*****
## **Assembly**

Back in ```ddocent_env``` directory
#### Created a new folder called ```RefOpt```
```javascript
mkdir RefOpt
cd RefOpt
```

#### Place a subset of individuals of the total data set in ```RefOpt```
```javascript
ln -s ../FBN_318.F.fq.gz
ln -s ../FBN_318.R.fq.gz
ln -s ../FBN_318.R1.fq.gz
ln -s ../FBN_318.R2.fq.gz

# Repeated this for FBN_329; FBS_40; FBS_62; OBN_023; OBN_2; OBS_241; OBS_258; PCN_209; PCN_220; PCS_347; PCS_364; SPN_378; SPN_399; SPS_86; SPS_93; WC1_408; WC1_419; WC2_282; WC2_302; WC3_469; WC3_485; WC4_437; WC4_460
# 2 individuals from each locality chosen randomly
```

### Scripts used for assembly steps are located in the [dDocent repository](https://github.com/jpuritz/dDocent/tree/master/scripts) on GitHub.

#### Run `ReferenceOpt.sh`
```javascript
curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/ReferenceOpt.sh
chmod +x ReferenceOpt.sh
```
```javascript
bash ./ReferenceOpt.sh 2 4 2 4 PE 16
```
```javascript

                                                   Histogram of number of reference contigs

                         3 +----------------------------------------------------------------------------------------+
                           |          +       *  + ****     +           +          +          +          +          |
                           |                  *    ****     'plot.kopt.data' using (bin($1,binwidth)):(1.0) ******* |
                           |                  *    ****                                                             |
                       2.5 |-+                *    ****                                                           +-|
                           |                  *    ****                                                             |
                           |                  *    ****                                                             |
                           |                  *    ****                                                             |
                         2 |-+        **   *****  ********           **       ***                                 +-|
                           |          **   *****  ****** *           **       * *                                   |
                           |          **   *****  ****** *           **       * *                                   |
                       1.5 |-+        **   *****  ****** *           **       * *                                 +-|
                           |          **   *****  ****** *           **       * *                                   |
                           |          **   *****  ****** *           **       * *                                   |
                           |          **   *****  ****** *           **       * *                                   |
                         1 |-+       ******************* ********************** ************************************|
                           |         ******************* ******** * ****** *  * *****  *  *  * **** * *  *   *   *  |
                           |         ******************* ******** * ****** *  * *****  *  *  * **** * *  *   *   *  |
                           |         ******************* ******** * ****** *  * *****  *  *  * **** * *  *   *   *  |
                       0.5 |-+       ******************* ******** * ****** *  * *****  *  *  * **** * *  *   *   *+-|
                           |         ******************* ******** * ****** *  * *****  *  *  * **** * *  *   *   *  |
                           |         ******************* ******** * ****** *  * *****  *  *  * **** * *  *   *   *  |
                           |         ******************* ******** * ****** *  * *****  *  *  *+**** * *  *   *   *  |
                         0 +----------------------------------------------------------------------------------------+
                         10000      15000      20000      25000       30000      35000      40000      45000      50000
                                                          Number of reference contigs

Average contig number = 25964.8
The top three most common number of contigs
X       Contig number
1       49899
1       47803
1       45897
The top three most common number of contigs (with values rounded)
X       Contig number
2       33400
2       22300
2       21700
```

#### Visualize data in `kopt.data` in R-Studio
```javasript
library(ggplot2)

data.table <- read.table("kopt.data", header = FALSE, col.names= c("k1","k2","Similarity", "Contigs"))

data.table$K1K2 <- paste(data.table$k1, data.table$k2, sep=",")

df=data.frame(data.table)
df$K1K2 <- as.factor(df$K1K2)

p <- ggplot(df, aes(x=Similarity, y=Contigs, group=K1K2)) + scale_x_continuous(breaks=seq(0.8,0.98,0.01)) + geom_line(aes(colour = K1K2))
p
```

![kopt.data](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/blob/master/Output/Assembly/koptdata.png)

Picked 0.9 as the similarity threshold as this is the inflection point in the curve.

#### Run `RefMapOpt.sh` using the similarity threshold of 0.9
Note- You will need to have the trimmed reads files *.R1.fq.gz and *.R2.fq.gz included to run this script
```javascript
curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/RefMapOpt.sh
chmod +x RefMapOpt.sh
```
```javascript
./RefMapOpt.sh 2 4 2 4 0.9 PE 20
```

This loops across cutoffs of 2-4 using a similarity of 90% for clustering, parellized across 20 processors, using PE assembly technique.

The output is stored in a file called mapping.results
```javascript
cat mapping.results
```
```javascript
Cov     Non0Cov Contigs MeanContigsMapped       K1      K2      SUM Mapped     SUM Properly     Mean Mapped     Mean Properly   MisMatched
84.101  136.566 43287   25752.5                 2       2       65530184        51280936       3.64057e+06      2.84894e+06     593921
133.911 166.359 25782   20031.3                 2       3       62147416        52853141       3.45263e+06      2.93629e+06     315434
153.002 177.66  21173   17459.8                 2       4       58314090        48027482       3.23967e+06      2.66819e+06     418717
99.8526 150.196 35357   22767.8                 3       2       63550548        49332981       3530586          2.74072e+06     556878
151.695 180.5   21808   17738                   3       3       59549846        46078804       3.30832e+06      2.55993e+06     558437
180.241 203.798 17950   15247.6                 3       4       58239203        44554680       3.23551e+06      2475260         584175
113.48  163.508 30308   20447                   4       2       61910409        50773700       3.43947e+06      2.82076e+06     434541
165.088 190.109 18865   15754.7                 4       3       56061748        46065600       3.11454e+06      2559200         387594
217.761 241.803 15353   13308.1                 4       4       60182907        53964523       3.34349e+06      2.99803e+06     207248
```
I chose 3 and 3. I ran dDocent for assembly on the subset.
```javascript
curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/dDocent
chmod +x dDocent
./dDocent
```
```javascript
dDocent 2.8.7

Contact jpuritz@uri.edu with any problems


Checking for required software

All required software is installed!

dDocent version 2.8.7 started Wed May 22 22:18:41 EDT 2019

24 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 24 individuals
dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
20

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
yes
What type of assembly would you like to perform?  Enter SE for single end, PE fo                                                                                        r paired-end, RPE for paired-end sequencing for RAD protocols with random sheari                                                                                        ng, or OL for paired-end sequencing that has substantial overlap.
Then press [ENTER]
PE
Reads will be assembled with Rainbow
CD-HIT will cluster reference sequences by similarity. The -c parameter (% similarity to cluster) may need to be changed for your taxa.
Would you like to enter a new c parameter now? Type yes or no and press [ENTER]
yes
Please enter new value for c. Enter in decimal form (For 90%, enter 0.9)
0.9
Do you want to map reads?  Type yes or no and press [ENTER]
no
Mapping will not be performed
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.


dDocent will require input during the assembly stage.  Please wait until prompt says it is safe to move program to the background.

                       Number of Unique Sequences with More than X Coverage (Counted within individuals)

  2.2e+06 +---------------------------------------------------------------------------------------------------------+
          |           +           +          +           +           +           +          +           +           |
    2e+06 |*+                                                                                                     +-|
          |*                                                                                                        |
          | *                                                                                                       |
  1.8e+06 |-*                                                                                                     +-|
          |  *                                                                                                      |
  1.6e+06 |-+ *                                                                                                   +-|
          |   *                                                                                                     |
  1.4e+06 |-+  *                                                                                                  +-|
          |    *                                                                                                    |
  1.2e+06 |-+   ***                                                                                               +-|
          |        **                                                                                               |
          |          *                                                                                              |
    1e+06 |-+         *****                                                                                       +-|
          |                *                                                                                        |
   800000 |-+               ******                                                                                +-|
          |                       *****                                                                             |
   600000 |-+                          ******                                                                     +-|
          |                                  ******************                                                     |
          |                                                    *****************************                        |
   400000 |-+                                                                               ************************|
          |           +           +          +           +           +           +          +           +           |
   200000 +---------------------------------------------------------------------------------------------------------+
          2           4           6          8           10          12          14         16          18          20
                                                           Coverage

Please choose data cutoff.  In essence, you are picking a minimum (within individual) coverage level for a read (allele) to be used in the reference assembly
3

                                 Number of Unique Sequences present in more than X Individuals

   160000 +---------------------------------------------------------------------------------------------------------+
          |*                   +                    +                     +                    +                    |
          | *                                                                                                       |
   140000 |-+*                                                                                                    +-|
          |   *                                                                                                     |
          |    *                                                                                                    |
   120000 |-+   *                                                                                                 +-|
          |      *                                                                                                  |
          |       *                                                                                                 |
   100000 |-+      *                                                                                              +-|
          |         *                                                                                               |
    80000 |-+        *****                                                                                        +-|
          |               ****                                                                                      |
          |                   *                                                                                     |
    60000 |-+                  *********                                                                          +-|
          |                             **                                                                          |
          |                               **********                                                                |
    40000 |-+                                       ***********                                                   +-|
          |                                                    ***********                                          |
          |                                                               *********************                     |
    20000 |-+                                                                                  *********************|
          |                                                                                                         |
          |                    +                    +                     +                    +                    |
        0 +---------------------------------------------------------------------------------------------------------+
          2                    4                    6                     8                    10                   12
                                                     Number of Individuals

Please choose data cutoff.  Pick point right before the assymptote. A good starting cutoff might be 10% of the total number of individuals
3
At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type 'bg' without the quotes and press enter
Type 'disown -h' again without the quotes and press enter

Now sit back, relax, and wait for your analysis to finish

dDocent assembled 83663 sequences (after cutoffs) into 21731 contigs

dDocent has finished with an analysis in /home/azyck/Fiddler_Crab/ddocent_env/RefOpt

dDocent started Wed May 22 22:18:41 EDT 2019

dDocent finished Wed May 22 22:24:47 EDT 2019
```

#### Run dDocent on the subset of data for Read Mapping
```javascript
dDocent 2.7.8

Contact jpuritz@gmail.com with any problems


Checking for required software

All required software is installed!

dDocent version 2.7.8 started Fri May 3 11:44:38 EDT 2019

24 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes

dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
20

dDocent detects 503 gigabytes maximum memory available on this system.
Please enter the maximum memory to use for this analysis in gigabytes.
For example, to limit dDocent to ten gigabytes, enter 10.
This option does not work with all distributions of Linux. If runs are hanging at variant calling, enter 0
Then press [ENTER]
0

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
yes

BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa.
Would you like to enter a new parameters now? Type yes or no and press [ENTER]
yes
Please enter new value for A (match score).  It should be an integer.  Default is 1.
1
Please enter new value for B (mismatch score).  It should be an integer.  Default is 4.
4
Please enter new value for O (gap penalty).  It should be an integer.  Default is 6.
6

Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.
```

#### Use ```samtools flagstat``` to investigate ```bam``` files.
https://github.com/bahlolab/bioinfotools/blob/master/SAMtools/flagstat.md
```javascript
samtools flagstat FBN_318-RG.bam
```

#### Output
```javascript
2213833 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
2213833 + 0 mapped (100.00% : N/A)
2213833 + 0 paired in sequencing
1134875 + 0 read1
1078958 + 0 read2
1992898 + 0 properly paired (90.02% : N/A)
2205126 + 0 with itself and mate mapped
8707 + 0 singletons (0.39% : N/A)
176622 + 0 with mate mapped to a different chr
167350 + 0 with mate mapped to a different chr (mapQ>=5)
```

I looked at a few other ```bam``` files. All had 100% mapped. A low percentage of the mappings were singletons (only one read from a pair), however quite a few reads had mates mapped to different reference contigs. Percentage of proper pairings was around 90%.

#### Re-run dDocent for Read Mapping, with a match score of ```1```, mismatch score of ```3```, and a gap penalty of ```5```.

#### Then check the stats again with ```samtools flagstat```.
```javascript
samtools flagstat FBN_318-RG.bam
```

#### Output
```javascript
2264622 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
2264622 + 0 mapped (100.00% : N/A)
2264622 + 0 paired in sequencing
1159182 + 0 read1
1105440 + 0 read2
2065116 + 0 properly paired (91.19% : N/A)
2255762 + 0 with itself and mate mapped
8860 + 0 singletons (0.39% : N/A)
154393 + 0 with mate mapped to a different chr
140298 + 0 with mate mapped to a different chr (mapQ>=5)
```

1, 3, and 5 produced slightly better stats than 1, 4, and 6, so I will use these parameters with the full dataset.

Copy the `reference.fasta` file from this `RefOpt` directory to the main working directory.
```javascript
cp reference.fasta ../
cd ../
```

Run dDocent on the full data set, skipping trimming and assembly.
```javascript
./dDocent
```
```javascript
dDocent 2.8.7

Contact jpuritz@uri.edu with any problems


Checking for required software

All required software is installed!

dDocent version 2.8.7 started Wed May 22 22:31:24 EDT 2019

376 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 376 individuals
dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
20

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
yes
BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa.
Would you like to enter a new parameters now? Type yes or no and press [ENTER]
yes
Please enter new value for A (match score).  It should be an integer.  Default is 1.
1
Please enter new value for B (mismatch score).  It should be an integer.  Default is 4.
3
Please enter new value for O (gap penalty).  It should be an integer.  Default is 6.
5
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
yes

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.
amaeliazyck@gmail.com


At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type 'bg' without the quotes and press enter
Type 'disown -h' again without the quotes and press enter

Now sit back, relax, and wait for your analysis to finish
```
```javascript
After filtering, kept 228044 out of a possible 711393 Sites
```

`dDocent` does minimal filtering. Using VCFtools, SNPs are filtered to only those that are called in 90% of all individuals. Further filtering must be completed using VCFtools.

dDocent produces a VCF file `TotalRawSNPs.vcf`, which will be used for additional [filtering](https://github.com/amaeliazyck/RADseq_Uca-rapax_2016/tree/master/Scripts/Filtering).
