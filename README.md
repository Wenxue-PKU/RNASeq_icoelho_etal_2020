# Repository for the RNASeq Analysis of the Paper "Trem-2 promotes emergence of restorative macrophages and endothelial cells during recovery from hepatic tissue damage" (Coelho et al., 2020)

This repository contains the raw data (counts files) and scripts for the analysis of the RNASeq for the referenced paper.

The first thing performed after getting the fastq files from the genomic facility, we performed Quality Control.

##### Quality Control

This quality control is done using the program **Fastqc** (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/); considering this program produces a report per file, we can also use a program to merge all the reports into one - **MultiQC** (https://multiqc.info/).

The reports for *the files corresponding to the 4 lanes separated* are saved in *fastqc_reports_perlane*. A summary of all reports is under the name *mQC_raw_perlane*.

The reports for *the files corresponding to the 4 lanes merged* are saved in *fastqc_reports_merged*. A summary of all reports is under the name *mQC_raw_merged*.

Looking at the overall reports, the samples look like to have a good overall quality, sequencing depth and profile (some of the plots will give out warnings, such as in %GC and Overrepresented sequences: these warnings are not uncommon in RNA-Seq datasets but, are not very worrying).

##### Alignment

For the alignment steps, we're going to use **STAR Alignment** (https://github.com/alexdobin/STAR)

###### 1) Create the Reference Genome Indexes

First of all, create a new folder to store the genome indexes and also create a shortcut for it:

```
mkdir gen_index
export gen_index=./gen_index/
```

Then, you can assess it by just using:

```
cd $gen_index
```

To align our samples, we need the reference genome indexes, that are just the corresponding reference genome (.fasta) and its corresponding annotation (.gtf). To get both of these files, we go to the Ensembl website (https://www.ensembl.org/info/data/ftp/index.html) and get the corresponding files:

```
wget ftp://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

wget ftp://ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.gtf.gz
gunzip Mus_musculus.GRCm38.95.gtf.gz
```

We are going to use the most recent version of the annotation (95 release). Now that we have the files, we proceed to use STAR with the option of  `genomeGenerate`

```
STAR --runThreadN 10 --runMode genomeGenerate --genomeDir $gen_index/ --genomeFastaFiles $gen_index/Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile $gen_index/Mus_musculus.GRCm38.95.gtf 

```

Options explained:

```
--runThreadN 10 \ #Nr. of Cores used
--runMode genomeGenerate \ #argument for the program to know what is going to run
--genomeDir $gen_index/ \ #path to the genome indexes
--genomeFastaFiles $gen_index/Mus_musculus.GRCm38.dna.primary_assembly.fa \ #path to the fasta file of the reference genome
--sjdbGTFfile $gen_index/Mus_musculus.GRCm38.95.gtf \ #path to the gtf file of the reference genome
```

Afterwards, you can proceed to the alignment step. For this, we will perform a cycle but, before, we need to create the several required folders to store the necessary files:

```
mkdir aligned_merged
mkdir qualimap_aligned_merged

```

Now, we go to the directory where the merged files are and proceed with the cycle:

```
cd raw_merged/
export gen_index=../../gen_index #this step helps in case of using a different window
```

```
for f in *.fastq.gz;
do STAR --genomeDir $gen_index --readFilesIn $f --readFilesCommand zcat --sjdbGTFfile $gen_index/Mus_musculus.GRCm38.95.gtf --quantMode GeneCounts --runThreadN 10 --outFileNamePrefix ../aligned_merged/"$f"_ --outSAMtype BAM SortedByCoordinate;
qualimap rnaseq -bam ../aligned_merged/"$f"_*.bam -gtf $gen_index/Mus_musculus.GRCm38.95.gtf -outdir ../qualimap_aligned_merged/$f --java-mem-size=20G;
rm -rf ../aligned/$f*.bam; done
```

This cycle will perform the alignment and perform the corresponding report for the file of the alignment.

After assessing the overall quality of the alignment using Qualimap & log files produced by STAR, we concluded that we could proceed with the DESeq2 Analysis.
