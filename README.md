<!-- dx-header -->
# RHtyper (DNAnexus Platform App)

Predict RHD/RHCE allele using Whole-Genome Sequencing data

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://wiki.dnanexus.com/.
<!-- /dx-header -->

<!-- Insert a description of your app here -->
> **RHtyper - prediction of RH alleles**

**Please read this important information before running the app**

**Introduction**

The *RHD* and *RHCE* genes encode Rh blood group antigens and exhibit extensive nucleotide polymorphisms and chromosome structural changes, often resulting in expression of Rh variant antigens. RH variation is common in patients with sickle cell disease (SCD) and can drive loss of antigen epitopes or expression of new epitopes, predisposing them to Rh alloimmunization. Serologic antigen typing is limited to common Rh antigens, necessitating a genetic approach to detect variant antigen expression. RHtyper is a novel algorithm developed for comprehensive RH genotyping from whole-genome sequencing (WGS) data.

**What does RHtyper do?**

This app implements RHtyper, a novel algorithm to predict RH allele pairs using WGS short reads data. It takes as input WGS alignment as BAM format. The BAM alignemnt can be generated using BWA aligner.


**What are typical use cases for RHtyper?**

Use this app to perform RH allele typing of samples sequenced using Illumina short read sequencing.


**How does RHtyper work?**
