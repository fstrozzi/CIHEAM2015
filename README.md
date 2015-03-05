# CIHEAM2015
Code and tutorials for CIHEAM course. QC, Variation calling and filtering sessions

# Index

[QC of NGS reads](https://github.com/fstrozzi/CIHEAM2015#qc-of-ngs-reads)

[Variation Calling](https://github.com/fstrozzi/CIHEAM2015#variation-calling)

[VCF Filtering](https://github.com/fstrozzi/CIHEAM2015#vcf-filtering)


QC of NGS reads
===============

### Exercise: Run FastQC to check your data

First of all prepare a directory for the QC practical:

```shell
mkdir QC
```

and create one sub-directory for FastQC:

```shell
mkdir -p QC/FastQC
```

Now enter the QC/FastQC directory

```shell
cd QC/FastQC
```

And run the FastQC analysis: 

```shell
perl /home/formacion/COMUNES/IAMZ/soft/FastQC-0.11.2/fastqc \
/home/formacion/COMUNES/IAMZ/data/CIHEAM/reads_FORTRIMMING/*.gz \
-o ./ --noextract -t 8
```

Command line explanation:

* the ```-o``` option specifiy the output folder where FastQC will write the results. In this case it is the current directory, specified with ```./```
* the ```--noextract``` tells FastQC to not automatically extract the zip files it creates
* the ```-t 8``` istruct FastQC to use 8 CPUs for this analysis

Now a look at the results:

[R1 reads](http://htmlpreview.github.io/?https://raw.githubusercontent.com/fstrozzi/CIHEAM2015/master/FastQC/lane1_NoIndex_L001_R1_001_fastqc.html)

[R2 reads](http://htmlpreview.github.io/?https://raw.githubusercontent.com/fstrozzi/CIHEAM2015/master/FastQC/lane1_NoIndex_L001_R2_001_fastqc.html)

If you want, you can also download all the FastQC reports just generated and view them locally:

```shell
scp -P 2222 "formacion09@calendula.fcsc.es:/home/formacion/formacionXX/QC/FastQC/lane*_NoIndex_L001_R*_001_fastqc.html" .
```

Now a look at some bad examples:

[Low quality reads](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)

[Adapter contamination](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/RNA-Seq_fastqc.html#M9)

### Exercise: Run Trimmomatic to filter your data

Trimmomatic is a powerful tool, which allows you process raw reads and remove low quality bases and adapter sequences.

Before running Trimmomatic, create a directory where it will create the output files and enter it:

```shell
mkdir -p QC/Trimmomatic
cd QC/Trimmomatic
```

To run Trimmomatic on one Sample, use this command line:

```shell
java -jar /home/formacion/COMUNES/IAMZ/soft/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 8 -phred33 \
/home/formacion/COMUNES/IAMZ/data/CIHEAM/reads_FORTRIMMING/lane1_NoIndex_L001_R1_001.fastq.gz \ 
/home/formacion/COMUNES/IAMZ/data/CIHEAM/reads_FORTRIMMING/lane1_NoIndex_L001_R2_001.fastq.gz Sample_1_R1_paired.fastq.gz \
Sample_1_R1_unpaired.fastq.gz Sample_1_R2_paired.fastq.gz Sample_1_R2_unpaired.fastq.gz LEADING:5 TRAILING:5 \ 
SLIDINGWINDOW:4:20 MINLEN:36 \ 
ILLUMINACLIP:/home/formacion/COMUNES/IAMZ/soft/Trimmomatic-0.33/adapters/TruSeq2-PE.fa:2:30:10
```

Command line explanations:

* ```-threads 8``` tells Trimmomatic to use 8 CPUs to speed up the analysis
* ```LEADING:5``` and ```TRAILING:5``` specify to remove bases having a quality score below 5 at the very beginning or end of the read
* ```SLIDINGWINDOW:4:20``` means Trimmomatic will scan the reads using a window of 4 bases and will trim the window if the average quality is below 20
* ```MINLEN:36``` discard reads shorter than 36nt after the trimming process
* ```ILLUMINACLIP``` this line specify which adapter file will be used to search and remove adaperts sequences within the reads

When it completes, look at the files generated:

```shell
$ ls
Sample_1_R1_paired.fastq.gz  Sample_1_R1_unpaired.fastq.gz  Sample_1_R2_paired.fastq.gz  Sample_1_R2_unpaired.fastq.gz
```

There are 4 files, two for paired-end reads and two for unpaired reads, i.e. reads that lost their mate due to trimming and filtering.

Variation Calling
=================================

### Exercise: Prepare a BAM file only with the reads mapped on chromosome 18

Before running the variation callers we need to subset to prepare a working folder:

```shell
mkdir variation_calling
cd  variation_calling
```

Now we need to prepare a BAM file with only the reads mapped on chromosome 18, to make the variation calling process quicker. This is quite simple and we are using Samtools for that

```shell
samtools view -b Sample_1.bam 18 > Sample_1.chr18.bam
```

Command line explanation (quite simple):

* ```-b``` specifies that the input file is a BAM file
* ```18``` just tell Samtools to print out only the reads on the chromosome 18 region

And of course we need to index this BAM file, for the callers to be able to extract the reads.

```shell
samtools index Sample_1.chr18.bam
```

### Exercise: Run FreeBayes

To run FreeBayes, we create a folder and we move there:

```shell
mkdir freebayes
cd freebayes
```
And now we can run the command

```shell
freebayes -b Sample_1.chr18.bam --min-mapping-quality 30 --min-alternate-count 5 --min-coverage 5 --min-base-quality 30 -f /home/formacion/COMUNES/IAMZ/data/CIHEAM/ReferenceGenome/bt_umd31/Bos_taurus.UMD3.1.fa -v Sample_1.fb.vcf -r 18
```

Command line explanation:

* ```-b``` specifies the input BAM file (in case of multiple samples, just repeat ```-b``` for every BAM file)
* ```--min-mapping-quality 30``` the minimum mapping quality in Phred scale to evaluate a position
* ```--min-alternate-count 5``` the minimum number of observation supporting an alternative allele to evaluate a position
* ```--min-coverage 5``` we want at least 5 reads mapped to evaluate a position
* ```--min-base-quality 30``` the minimum base quality to evaluate a position
* ```-f``` this is the Fasta file with the reference genome
* ```-v``` the name of the output VCF file
* ```-r 18``` just restrict the analysis to chromsome 18 


### Exercise: Run GATK UnifiedGenotyper

We create a folder for GATK analysis and we move there:

```shell
cd $HOME/variant_calling 
mkdir gatk
cd gatk
```

And we run the command line for GATK UnifiedGenotyper:

```shell
java -Xmx8G -jar /home/formacion/COMUNES/IAMZ/soft/GATK-3.3.0/GenomeAnalysisTK.jar -T UnifiedGenotyper \
-R /home/formacion/COMUNES/IAMZ/data/CIHEAM/ReferenceGenome/bt_umd31/Bos_taurus.UMD3.1.fa \
-I Sample_1.chr18.bam \
-o Sample_1.ug.vcf \
--min_base_quality_score 30 \
-stand_call_conf 30 \
-stand_emit_conf 30 \
-L 18 \
-nt 8 \
-glm BOTH
```

Command line explanation:

* ```-T``` this specify the tool to use from GATK, in this case UnifiedGenotyper
* ```-R``` the Fasta file with the reference genome
* ```-I``` the input BAM file
* ```-o``` the output VCF file
* ```--min_base_quality_score 30``` minimum base quality to evaluate a position
* ```-stand_call_conf 30``` minimum quality to call a variant
* ```-stand_emit_conf 30``` minimum quality to emit (i.e. write in the VCF file) a variant
* ```-L 18``` just process reads on chromosome 18
* ```-nt 8``` uses 8 CPUs to speed up the calculations
* ```-glm BOTH``` Genotype likelihoods model to use. BOTH means call SNPs and InDels

### Exercise: Run Samtools Variation Calling

We create a folder for Samtools and we move there:

```shell
cd $HOME/variant_calling
mkdir samtools
cd samtools
```

To run Samtools variation calling we use this command line:

```shell
samtools mpileup -t DV -q30 -Q30 \
-u -f /home/formacion/COMUNES/IAMZ/data/CIHEAM/ReferenceGenome/bt_umd31/Bos_taurus.UMD3.1.fa Sample_1.chr18.bam\
| bcftools call -m -v > Sample_1.st.vcf
```

Command line explanation:
* ```-t``` tells Samtools to include also the DV field into the output VCF file
* ```-u``` means Samtools will emit output as uncompressed. Since we are piping this output to Bcftools we do not waste time in compressing and uncompressing again the stream of data
* ```-f``` this is the Fasta file with the reference genome
* ```bcftools call -m -v``` Bcftools takes the output of Samtools and call the variants, with options: 
    * ```-m``` multiallelic caller (new algorithm)
    * ```-v``` emits only variants sites


VCF Filtering
=============

### Exercise: Select only SNPs from a VCF file

We first create a working directory for the filtering activities:

```shell
cd $HOME/variant_calling
mkdir filtering
```

As a first exercise, we filter the raw VCF files to retain only the SNPs:

```shell
java -Xmx4G -jar /home/formacion/COMUNES/IAMZ/soft/GATK-3.3.0/GenomeAnalysisTK.jar -T SelectVariants \
-R /home/formacion/COMUNES/IAMZ/data/CIHEAM/ReferenceGenome/bt_umd31/Bos_taurus.UMD3.1.fa \
--variant Sample_1.ug.vcf \
-o Sample_1.ug.SNP.vcf \
-selectType SNP
```

Command line explanation:

* ```-T``` specifies the GATK tool to use, in this case SelectVariants
* ```-R``` the Fasta file with the reference genome
* ```--variant``` is the input VCF file
* ```-o``` the output VCF file with filtered variants
* ```-selectType``` specifies the type of variant to be retained, SNP in this case


After this exercise is completed, do the same thing for the VCF files generated with FreeBayes and Samtools. Remeber to keep the same convention for file names.

### Exercise: Compress and index a VCF file

After having filtered the raw VCF files to retain only SNPs, we will proceed with the VCF compression and indexing, which are necessary for the variations filtering steps:


```shell
bgzip Sample_1.ug.SNP.vcf
```

```shell
tabix -p vcf Sample_1.ug.SNP.vcf.gz
```

These command lines are quite simple:

* ```bgzip``` is a compression tool designed specifically for VCF files and other common file formats (e.g. GFF)
* ```tabix``` is an indexing tool that process a bgzip compressed file. The only option needed here is the ```-p``` which specifies the file format, VCF in this case

After this exercise is completed, do the same thing for the VCF files generated with FreeBayes and Samtools. Remeber to keep the same convention for file names.

### Exercise: Filter the VCF file according to different parameters

Now the real filtering. To do this we are going to use the Bcftools software and we experiment with different filters combinations:

**Simple filters using BCFtools**

```shell
bcftools filter -i'DP>100' ../freebayes/Sample_1.fb.SNP.vcf.gz > Sample_1.fb.SNP.filtered.vcf 
```

```shell
bcftools filter -i'%QUAL>100' ../freebayes/Sample_1.fb.SNP.vcf.gz > Sample_1.fb.SNP.filtered.vcf
```

```shell
bcftools filter -i'%QUAL>100 && DP>10' ../freebayes/Sample_1.SNP.vcf.gz > Sample_1.fb.SNP.filtered.vcf
```

In this case we are filtering using the quality and depth of coverage parameters:

After each fileter, use this command line to assess the number of remaing SNPs:

```shell
grep -c -v "#" Sample_1.fb.SNP.filtered.vcf
```

A more detailed information can also be generated using another utility from Bcftools

```shell
bcftools stats Sample_1.fb.SNP.filtered.vcf
```

This is a long report on various metrics, including the number of nucleotide transitions and transversions, the quality and coverage distributions etc.

Try to experiment for a while using ```bcftools filter``` also with the other VCF files from GATK and Samtools. Perform some filtering using also specifics parameters for each caller, e.g the DV field for Samtools or the QD field for GATK.

At the end of this part, you need to have one filtered VCF file for each caller.

### Exercise: Run vcf-compare to get statistics on different VCF files

Now we perform a comparison of the SNPs that are in common between the 3 filtered datasets generated with FreeBayes, GATK and Samtools:

```shell
vcf-compare Sample_1.st.SNP.filtered.vcf Sample_1.ug.SNP.filtered.vcf Sample_1.fb.SNP.filtered.vcf | grep "^VN" | cut -f 2-
```

Command line explanation:

* ```vcf-compare``` just takes one or more VCF file and output information on the comparison
* ```| grep "^VN" | cut -f 2-``` these commands are just used to extract the information suitable to generate a Venn diagram

### Exercise: Plot the results of vcf-compare using R

**This part is done locally on your PC and NOT on the cluster**

Open a terminal and run R:

```shell
$ R

R version 3.1.1 (2014-07-10) -- "Sock it to Me"
Copyright (C) 2014 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin13.1.0 (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

>
```

From now on all the commmands will be issued on the R shell

First, install the ```venneuler``` package:

```R
install.packages("venneuler")
```
if you are promped with the selection of a location where to download the package, just select "Spain" and confirm.

Then we create and plot a Venn Diagram:

```R
library(venneuler)
v = venneuler(c(FB=66,"FB&ST"=181,"FB&UG"=1799,ST=2727,UG=8330,"ST&UG"=47455,"FB&ST&UG"=125625))
plot(v)
```

This is just an exmaple to show how to pass data to this package. For a good practice, you need to take the numbers generated from your ```vcf-compare``` analysis. For simplicity just use short codes for each caller:

* FB: FreeBayes
* UG: GATK UnifiedGenotyper
* ST: Samtools

As you can see from the example, the union of different callers results is defined using an ```&```.

At the end of this exercise you should get a Venn diagram similar to this one:

![alt text](http://i.imgur.com/gQznCm1.jpg "Venn Diagram of variation calls")







