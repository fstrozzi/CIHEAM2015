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
perl /home/formacion/COMUNES/IAMZ/soft/FastQC-0.11.2/fastqc /home/formacion/COMUNES/IAMZ/data/CIHEAM/reads_FORTRIMMING/*.gz -o ./ --noextract -t 8
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
SLIDINGWINDOW:4:20 MINLEN:36 ILLUMINACLIP:/home/formacion/COMUNES/IAMZ/soft/Trimmomatic-0.33/adapters/TruSeq2-PE.fa:2:30:10
```

When if finish, look at the files generated:

```shell
$ ls
Sample_1_R1_paired.fastq.gz  Sample_1_R1_unpaired.fastq.gz  Sample_1_R2_paired.fastq.gz  Sample_1_R2_unpaired.fastq.gz
```

There are 4 files, two for paired-end reads and two for unpaired reads, i.e. reads that lost their mate due to trimming and filtering.

Variation Calling
=================================

### Exercise: Prepare a BAM file only with the reads mapped on chromosome 18

```shell
samtools view -b Sample_1.md.sort.bam 18 > Sample_1.chr18.bam
```

```shell
samtools index Sample_1.chr18.bam
```

### Exercise: Run FreeBayes

```shell
freebayes -b Sample_1.chr18.bam --no-mnps --min-mapping-quality 30 --min-alternate-count 5 --min-coverage 5 --min-base-quality 30 -f /home/formacion/COMUNES/IAMZ/data/CIHEAM/ReferenceGenome/bt_umd31/Bos_taurus.UMD3.1.fa -v Sample_1.fb.vcf -r 18
```

### Exercise: Run GATK UnifiedGenotyper

```shell
java -Xmx8G -jar /home/formacion/COMUNES/IAMZ/soft/GATK-3.3.0/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/formacion/COMUNES/IAMZ/data/CIHEAM/ReferenceGenome/bt_umd31/Bos_taurus.UMD3.1.fa -I Sample_1.chr18.bam -o Sample_1.ug.vcf --min_base_quality_score 30 -stand_call_conf 30 -stand_emit_conf 30 -L 18 -nt 8 -glm BOTH
```

### Exercise: Run Samtools Variation Calling

```shell
samtools mpileup -q30 -Q30 -uf /home/formacion/COMUNES/IAMZ/data/CIHEAM/ReferenceGenome/bt_umd31/Bos_taurus.UMD3.1.fa Sample_1.chr18.bam | bcftools call -mv > Sample_1.st.vcf
```

VCF Filtering
=============

### Exercise: Select only SNPs from a VCF file

```shell
java -Xmx4G -jar /home/formacion/COMUNES/IAMZ/soft/GATK-3.3.0/GenomeAnalysisTK.jar -T SelectVariants -R /home/formacion/COMUNES/IAMZ/data/CIHEAM/ReferenceGenome/bt_umd31/Bos_taurus.UMD3.1.fa --variant Sample_1.ug.vcf -o Sample_1.ug.SNP.vcf -selectType SNP

java -Xmx4G -jar /home/formacion/COMUNES/IAMZ/soft/GATK-3.3.0/GenomeAnalysisTK.jar -T SelectVariants -R /home/formacion/COMUNES/IAMZ/data/CIHEAM/ReferenceGenome/bt_umd31/Bos_taurus.UMD3.1.fa --variant Sample_1.fb.vcf -o Sample_1.fb.SNP.vcf -selectType SNP

java -Xmx4G -jar /home/formacion/COMUNES/IAMZ/soft/GATK-3.3.0/GenomeAnalysisTK.jar -T SelectVariants -R /home/formacion/COMUNES/IAMZ/data/CIHEAM/ReferenceGenome/bt_umd31/Bos_taurus.UMD3.1.fa --variant Sample_1.st.vcf -o Sample_1.st.SNP.vcf -selectType SNP
```

### Exercise: Compress and index a VCF file

```shell
bgzip Sample_1.fb.SNP.vcf
bgzip Sample_1.ug.SNP.vcf
bgzip Sample_1.st.SNP.vcf
```

```shell
tabix -p vcf Sample_1.fb.SNP.vcf.gz
tabix -p vcf Sample_1.ug.SNP.vcf.gz
tabix -p vcf Sample_1.st.SNP.vcf.gz
```

### Exercise: Filter the VCF file according to different parameters


### Exercise: Run vcf-compare to get statistics on different VCF files

```shell
vcf-compare Sample_1.st.SNP.vcf.gz Sample_1.ug.SNP.vcf.gz Sample_1.fb.SNP.vcf.gz | grep "^VN" | cut -f 2-
```

### Exercise: Plot the results of vcf-compare using R

```R
library(venneuler)
v = venneuler(c(FB=66,"FB&HC"=181,"FB&UG"=1799,HC=2727,UG=8330,"HC&UG"=47455,"FB&HC&UG"=125625)
plot(v)
```

