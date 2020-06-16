<!-- dx-header -->
# RHtyper 

Predict *RHD*/*RHCE* allele using Whole-Genome Sequencing data

<!-- /dx-header -->

<!-- Insert a description of your app here -->
> **RHtyper - prediction of RH alleles**

**Please read this important information before running the app**

## **Introduction**

The *RHD* and *RHCE* genes encode Rh blood group antigens and exhibit extensive nucleotide polymorphisms and chromosome structural changes, often resulting in expression of Rh variant antigens. RH variation is common in patients with sickle cell disease (SCD) and can drive loss of antigen epitopes or expression of new epitopes, predisposing them to Rh alloimmunization. Serologic antigen typing is limited to common Rh antigens, necessitating a genetic approach to detect variant antigen expression. RHtyper is a novel algorithm developed for comprehensive RH genotyping from whole-genome sequencing (WGS) data.

## **What does RHtyper do?**

This app implements RHtyper, a novel algorithm to predict RH allele pairs using WGS short reads data. It takes as input WGS alignment as BAM format. The BAM alignemnt can be generated using BWA aligner.


## **What are typical use cases for RHtyper?**

Use this app to perform RH allele typing of samples sequenced using Illumina short read sequencing.

## **How does RHtyper work?**

RHtyper performs the following main steps:

Variant calling in the *RHD* and *RHCE* gene regions.
RHD and allele zygosity determination using coverage profiles and variant information.
Likelihood scoring and ranking to predict the RH allele pairs


## **What data are required for RHtyper to run?**

1. RHtyper requires:
   * WGS short reads alignment in BAM format (bam)
   * Index of the BAM alignment file (bai)
   * Version of the reference genome
   * Estimated coverage depth of the WGS data
   * Output prefix. The prefix added to the output files
2. Other adjustable parameters/options:
   * Gene for typing. User can select to type only *RHD* or *RHCE* gene
   * Reference sequence. User can provide customized genome reference sequence to be consistent with the one used for alignment
   * Cutoff of alternative read number. The lower bound of read number with alternative nucleotide to call variant
   * Allele linking database. User can provide customized *RHD*-*RHCE* allele linking information to improve the accuracy of the typing. If not specified, RHtyper will use linking information from ISBT database
   * Allele population frequency database. User can provide allele population frequency information to improve the accuracy of the typing. If not specified, RHtyper will use internal database
   * Verbose level. Level of log details to report.

## **What data are the format for Allele linking database and Allele population frequency database?**

Both databases are in tab-delimited text format and can be created using Excel (with headers). Examples are shown below.

1. Allele linking database

| Gene | Allele name | Allele detail | Alias | Linked | comment |
| :---: | :---: | :---: | :---: | :---: | :---: |
| RHCE | RHCE\*01.02.01 RHCE\*ce.02.01 | RHCE\*ce48C,1025T | RHCE\*ceTI | RHD\*04.01_RHD\*DIVa ||
| RHCE | RHCE\*01.05.01 RHCE\*ce.05.01 | RHCE\*ce48C,712G,787G,800A | RHCE\*ceEK | RHD\*DAR ||
| RHD | RHD\*04.01 RHD\*DIVa | RHD\*186T,410T,455C,1048C | | RHCE\*01.02.01 RHCE\*ce.02.01 ||
| ... | ... | ... | ... | ... | ... |

2. Allele population frequency database

| Allele name	| Nucleotide | Allele detail | Alias | PopFreq |
| :---: | :---: | :---: | :---: | :---: |
| RHD\*10.00 RHD\*DAU0 | c.1136C>T | RHD\*1136T | RHD\*1136T | 0.1651 |
| RHD\*09.03.01 RHD\*DAR3.01 | c.602C>G; c.667T>G; c.819G>A |	RHD\*602G,667G | RHD\*602G,667G | 0.0298 |
| RHCE\*02 or RHCE\*Ce RHCE\*C RHCE\*e | c.48G>C | | RHCE\*Ce | 0.119 |
| ... | ... | ... | ... | ... |

## **What does RHtyper output?**

The outputs incldue the following files:

1. Main output:
   * *bloodtyping.pdf* - the PDFs contained information, including typed alleles, called variants and coverage profiles for the sample
   * *bloodtyping.xlsx* - Information of typed allele pairs and allele pair ranking in Excel format
2. Other outputs can be used for cohort level analyses:
   * *bloodtyping.txt* - Information of typed allele pairs in tab-delimited text format
   * *exonCNV.txt* - Copy number variation status per exon in tab-delimited text format
   * *final.variant.txt* - variant identified in tab-delimited text format

## **How to run RHtyper?**
**NOTE**: the latest DNAnexus web interface is used for the following workflow

1. Create a new project. If you have an existing project, please skip to step 5.

![Step1](tutorial/S1.png)

2. Name the New Project and click Create Project

![Step2](tutorial/S2.png)

3. Upload BAM and BAI files: Under MANAGE, click Add, and select Upload Data

![Step3](tutorial/S3.png)

4. Upload BAM and BAI files following the instructions on the screen 

![Step4](tutorial/S4.png)

5. Find RHtyper: Click the search icon and type “RHtyper”

![Step5](tutorial/S5.png)

6. Select RHtyper

![Step6](tutorial/S6.png)

7. Install RHtyper for the first-time user. **Note**: introduction and parameter settings of RHtyepr can be found on this page

![Step7](tutorial/S7.png)

8. Click "Run" button to set up a new analysis

![Step8](tutorial/S8.png)

9. Select the project created at step1 or other projects containing the alignment BAM files for analyses.

![Step9](tutorial/S9.png)

10. Select a folder for output results: Under Analysis Settings, click Execution Output Folder

![Step10](tutorial/S10.png)

11. Create an output folder: Create a New Folder or select an exiting folder, then click Select Folder

![Step11](tutorial/S11.png)









10. Output is saved in the folder specified at step8

![Step10](tutorial/output.png)
