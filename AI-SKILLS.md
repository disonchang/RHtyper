🧬 Tool Overview: RHtyper
RHtyper is designed to automate the typing of RH alleles (RHD and RHCE) using Whole Genome Sequencing (WGS) or Whole Exome Sequencing (WES) data. This document outlines how AI agents and LLMs should interact with, interpret, and troubleshoot this tool within bioinformatics workflows.

🤖 LLM Compatibility & Integration
This tool is designed to be LLM-Agnostic. Whether using GPT-4, Claude 3.5, or specialized models like BioGPT, the following integration standards apply:
1. Prompt Engineering: LLMs can generate CLI commands for RHtyper based on user-provided file paths (.bam, .cram, .vcf).
2. Workflow Automation: Compatible with agentic frameworks (LangChain, Nextflow) to chain RHtyper output into downstream blood group phenotyping reports.
3. Standardized I/O: Supports JSON and TSV outputs for easy parsing by LLM "tools" or "functions."

✨ Core AI Capabilities
1. Read Depth Analysis: The AI can interpret RHtyper’s copy number variation (CNV) calls to distinguish between RHD deletion and RHD pseudogenes.
2. Nomenclature Translation: Translates RHtyper raw allele calls into standard ISBT (International Society of Blood Transfusion) terminology.
3. RH Allele Prediction: Assists users in determining the most likely allele pairs based on allele likelihood provided by the tool. 
4. Zygosity Prediction: Assists users in determining the likelihood of homozygous vs. heterozygous states based on allele frequency data provided by the tool.
5. Quality Control (QC): Identifies low-coverage regions in WES data that may lead to "No Calls" in specific RH exons.

🚧 Limitations & Safety
1. Clinical Validation: AI-generated interpretations of RHtyper results must be validated by a molecular biologist before clinical use.
2. WES Blind Spots: AI should warn users that WES data may miss deep intronic variants that RHtyper might otherwise detect in WGS data.
3. Hallucination Check: Always verify allele names against the Rhesus Base or similar databases; LLMs may occasionally hallucinate non-existent allele strings.

🛠️ AI Troubleshooting Guide
| Symptom | AI Diagnostic Logic | Recommended Fix |
| :--- | :--- | :--- |
| **Low Exonic Coverage** | Detect if input is WES vs WGS. | WGS has lower coverage; decrease `--coverage`, or use WGS mode |
| **Ref Mismatch** | Verify if BAM matches GRCh37 or GRCh38. | Match the tool reference using the `--genome_build` flag. |
| **Zygosity Errors** | Analyze B-Allele Frequency (BAF) ratios. | Adjust the `--altN` parameter (e.g., 3 or 5), or the --heterozyoteCov (e.g. 0.1, 0.2). |
