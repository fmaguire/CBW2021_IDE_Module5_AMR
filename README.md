### TO BE EDITED

## Table of contents
1. [Introduction](#intro)
2. [RGI for Genome Analysis](#rgi)
3. [RGI at the Command Line](#rgicommand)
4. [RGI for Merged Metagenomics Reads](#rgimerged)
5. [Metagenomics Data and the Burrows-Wheeler Transform](#bwt)
6. [Kmers and Pathogen-of-Origin for Metagenomics](#pathID)

<a name="intro"></a>
## Introduction

This module gives an introduction to prediction of antimicrobial resistome and phenotype based on comparison of genomic or metagenomic DNA sequencing data to reference sequence information. While there is a large diversity of reference databases and software, this tutorial is focused on the Comprehensive Antibiotic Resistance Database ([CARD](http://card.mcmaster.ca)) for genomic AMR prediction.

The relationship between AMR genotype and AMR phenotype is complicated and no tools for complete prediction of phenotype from genotype exist. Instead, analyses focus on prediction or catalog of the AMR resistome – the collection of AMR genes and mutants in the sequenced sample. While BLAST and other sequence similarity tools can be used to catalog the resistance determinants in a sample via comparison to a reference sequence database, interpretation and phenotypic prediction are often the largest challenge. To start the tutorial, we will use the Comprehensive Antibiotic Resistance Database ([CARD](http://card.mcmaster.ca)) to examine the diversity of resistance mechanisms, how they influence bioinformatics analysis approaches, and how CARD’s Antibiotic Resistance Ontology (ARO) can provide an organizing principle for interpretation of bioinformatics results.

CARD’s website provides the ability to: 

* Browse the Antibiotic Resistance Ontology (ARO) and associated knowledgebase.
* Browse the underlying AMR detection models, reference sequences, and SNP matrices.
* Download the ARO, reference sequence data, and indices in a number of formats for custom analyses.
* Perform integrated genome analysis using the Resistance Gene Identifier (RGI).

In this part of the tutorial, your instructor will walk you through the following use of the CARD website to familiarize yourself with its resources:

* re-evaluate

<a name="#rgi"></a>
## RGI for Genome Analysis

The diversity of antimicrobial resistance mechanisms requires a diversity of detection algorithms and a diversity of detection limits. CARD’s Resistance Gene Identifier (RGI) currently integrates four CARD detection models: **Protein Homolog Model**, **Protein Variant Model**, **rRNA Variant Model**, and **Protein Overexpression Model**. Unlike naïve analyses, CARD detection models use curated cut-offs, currently based on BLAST/DIAMOND bitscore cut-offs. Many other available tools are based on BLASTN or BLASTP without defined cut-offs and avoid resistance by mutation entirely. 

In this part of the tutorial, we will walk you through the following use of CARD’s Resistome Gene Identifier with “Perfect and Strict hits only”:

* Resistome prediction for the multidrug resistant *Acinetobacter baumannii* MDR-TJ, complete genome (NC_017847)
* Resistome prediction for the plasmid isolated from *Escherichia coli* strain MRSN388634 plasmid (KX276657)
* Explain the difference in triclosan resistance between two clinical strains of *Pseudomonas aeruginosa* that appear clonal based on identical MLST (Pseudomonas1.fasta, Pseudomonas2.fasta; these are SPAdes assemblies)
 
<a name="rgicommand"></a>
## RGI at the Command Line

RGI is a command line tool as well, so we’ll do a demo analysis of 33 . We’ll additionally try RGI’s heatmap tool for visualization of the results.

Login into your course account’s working directory, make a module 5 directory, and start the conda environment with RGI installed (`amr`):

```bash
mkdir module5
cd module5
conda activate amr
```

Take a peak at the list of three species samples and the options for the RGI software:

```bash
ls /home/ubuntu/CourseData/IDGE_data/module5/neisseria
ls /home/ubuntu/CourseData/IDGE_data/module5/psuedomonas
ls /home/ubuntu/CourseData/IDGE_data/module5/salmonella
rgi -h
```

First we need to acquire the latest AMR reference data from CARD:
Note: `rgi auto_load` to replace (`--local`)

```bash
rgi load -h
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
less card.json
rgi load -i card.json --local
ls
```

We don’t have time to analyze all 33 samples, so let’s analyze 1 as an example:

Pick one of the contigs. Or select one set of reads from either Neisseria, Pseudomonas, or Salmonella

```bash

# shovill --> if running from reads otherwise,

# using contig
rgi main –h
rgi main -i /home/ubuntu/CourseData/IDGE_data/module4/ecoli/SAM*.fasta -o single_sample -t contig -a BLAST -n 4 --local --clean
ls
less single_sample.json
less single_sample.txt
```

Using pre-compiled results for all 33 samples, so let’s try RGI’s beta heat map tool (pre-compiled images can be downloaded from the course GitHub repo):

```bash
ls /home/ubuntu/CourseData/IDGE_data/module4/ecoli_json
rgi heatmap –h
rgi heatmap -i precompiled_resources/ -cat gene_family -o genefamily_samples -clus samples
rgi heatmap -i  precompiled_resources -cat drug_class -o drugclass_samples -clus samples
rgi heatmap -i precompiled_resources  -o cluster_both -clus both
rgi heatmap -i precompiled_resources  -o cluster_both_frequency -f -clus bothls
```
