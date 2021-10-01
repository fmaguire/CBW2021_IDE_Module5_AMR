### TO BE EDITED

## Table of contents
0. [Learning Objectives](#obj)
1. [Introduction](#intro)
2. [Using CARD](#card)
3. [AMR Detection from Unassembled Reads](#rgibwt)
4. [AMR Detection from Assembled Contigs](#rgi)
5. [hAMRonizing AMR Detection Results](#hamronization)
6. [Predicting Resistance Phenotype](#heatmap)

<a name="obj"></a>
## Learning Objectives

By the end of this practical you will be able to:

* Use the Comprehensive Antibiotic Resistance Database ([CARD](http://card.mcmaster.ca)] to find information about AMR genes using the Antibiotic Resistance Ontology (ARO).

* Identify AMR genes in unassembled bacterial reads using the Resistance Gene Identifier's (RGI) read-mapping mode

* Assemble bacterial genome reads into contigs and use RGI to identify AMR genes

* Use [hAMRonization](https://github.com/pha4ge/hAMRonization) to compare and report the results of AMR gene detection in a standardised format

* Understand the difficulties of predicting phenotype from AMR genomics.

<a name="intro"></a>
## Introduction

The relationship between AMR genotype and AMR phenotype is complicated and no tools for complete prediction of phenotype from genotype exist. 
Instead, analyses focus on prediction or catalog of the AMR resistome – the collection of AMR genes and mutants in the sequenced sample. 
While BLAST and other sequence similarity tools can be used to catalogue the resistance determinants in a sample via comparison to a reference sequence database, interpretation and phenotypic prediction are often the largest challenge. 

There are several databases which try and organise and collect information about AMR to try and help with these challenges.
Many of these are either specialised on a specific type of resistance gene (e.g., [beta-lactamases](http://bldb.eu/)), organism (e.g., [_Mycobacterium tuberculosis_](https://github.com/jodyphelan/tbdb)), or are an automated amalgamation of other databases (e.g., [MegaRes](https://megares.meglab.org/)). 
There are also many tools for detecting AMR genes each with their own strengths and weaknesses (see [this paper] for a non-comprehensive list!).

The "big 3" databases that are both comprehensive (involving many organisms, genes, and types of resistance), regularly updated, have their own gene identification tool(s), and are carefully maintained and curated are: 

1. Comprehensive Antibiotic Resistance Database ([CARD](https://card.mcmaster.ca)) with the Resistance Gene Identifier ([RGI](https://github.com/arpcard/rgi))
2. National Center for Biotechnology Information's National Database of Antibiotic Resistant Organisms ([NDARO](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/)) with [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/)
3. [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/) database with its associated [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) tool

There are efforts to try to ensure [co-ordination of curation efforts](https://github.com/arpcard/amr_curation/) between these databases.

In this practical we are going to focus on CARD and the associated RGI tool because: the Antibiotic Resistance Ontology is an effective way to organise information about AMR, we are biased as CARD is Canadian and pretty much all the CBW faculty collaborate (or are part of) the group that develops CARD!

<a name="card"></a>
## Using CARD

To start this tutorial, we will use [CARD](http://card.mcmaster.ca) to examine the diversity of resistance mechanisms, how they influence bioinformatics analysis approaches, and how CARD's Antibiotic Resistance Ontology (ARO) can provide an organizing principle for interpretation of bioinformatics results.

CARD's website provides the ability to: 

* Browse the Antibiotic Resistance Ontology (ARO) and associated knowledge base.
* Browse the underlying AMR gene detection models, reference sequences, and SNP matrices.
* Download the ARO, reference sequence data, and indices in a number of formats for custom analyses.
* Detect AMR genes in a genome using a web-based RGI.

In this part of the tutorial, your instructor will walk you through the following use of the CARD website to familiarize yourself with its resources.

1. What are the mechanism of resistance described in the Antibiotic Resistance Ontology?
2. Examine the NDM-1 beta-lactamase protein, it’s mechanism of action, conferred antibiotic resistance, it’s prevalence, and it’s detection model. 
3. Examine the aac(6')-Iaa aminoglycoside acetyltransferase, it’s mechanism of action, conferred antibiotic resistance, it’s prevalence, and it’s detection model. 
4. Examine the fluoroquinolone resistant gyrB for M. tuberculosis, it’s mechanism of action, conferred antibiotic resistance, and it’s detection model. 
5. Examine the glycopeptide resistance gene cluster VanA, it’s mechanism of action, conferred antibiotic resistance, and it’s detection model(s). 
6. Examine the MexAB-OprM efflux complex, it’s mechanism of action, conferred antibiotic resistance, it’s prevalence, and it’s detection model(s). 

<a name="#rgibwt"></a>
## AMR Detection from Unassembled Reads

The diversity of antimicrobial resistance mechanisms requires a diversity of detection algorithms and a diversity of detection limits. CARD’s Resistance Gene Identifier (RGI) currently integrates four CARD detection models: **Protein Homolog Model**, **Protein Variant Model**, **rRNA Variant Model**, and **Protein Overexpression Model**. Unlike naïve analyses, CARD detection models use curated cut-offs, currently based on BLAST/DIAMOND bitscore cut-offs. Many other available tools are based on BLASTN or BLASTP without defined cut-offs and avoid resistance by mutation entirely. 

Using these cut-offs RGI assigns a criteria to the AMR genes its detects: 

* Perfect: detected gene identically matches the reference gene in CARD
* Strict: detected gene BLAST/DIAMOND bitscore is more than or equal to the curated cut-off for that detection model
* Loose (if run with `--include_loose`): potentially questionable hits below curated bitscore cut-offs

RGI was originally designed to only work on assembled genomic data (using BLAST or DIAMOND) but now incudes a feature to run on unassembled reads (using read-mapping approaches).

There are many features and options in the RGI tool. You can explore these by looking at `rgi --help` in the terminal or reading the [documentation](https://github.com/arpcard/rgi).

First thing we are going to do is set up RGI and then use it to identify AMR genes in unassembled bacterial reads.

### Setting-up

Most of the tool's we will use in this session have already been installed in a conda environment in your machine images.

So first let's create a folder to contain analyses from today and activate that conda environment:

```bash
mkdir module5
cd module5
conda activate amr
```

AMR is a fast moving field, and CARD uses a combination of auomated literature mining and manual expert curation to keep up.
If possible you always want to be using the most recent database versions for analyses.
The exception to this rule is when you are trying to reproduce an earlier analysis or ensure consistency in predictions.

Either way, it is always a good idea to check what database has been installed with RGI (or any other tool) before you run it.

Running this command will tell you the database versions currently being used by RGI:

```bash
rgi database --all --version
``` 
You can then go to the [CARD download page](https://card.mcmaster.ca/download) to double-check the latest `card_canonical` database (called "CARD Data" on the download page) and `card_variants` (called "CARD Prevalence, Resistomes, & Variants data" on the download page).

You shouldn't need to do this here but if you ever want to load a new database you can either do it manually with `rgi load` (instructions in the RGI [documentation](https://github.com/arpcard/rgi#rgi-databases)) or via the `rgi autoload` feature.

So we've got our tool, we've got our reference database, now we just need input data! 

We've pre-loaded your instance with data from 3 different bacterial species: _Neisseria gonorrhoeae_, _Pseudomonas aeruginosa_, and _Salmonella enterica_.

You can look at this data by typing:

```bash
ls /home/ubuntu/CourseData/IDGE_data/module5
ls /home/ubuntu/CourseData/IDGE_data/module5/neisseria
ls /home/ubuntu/CourseData/IDGE_data/module5/pseudomonas
ls /home/ubuntu/CourseData/IDGE_data/module5/salmonella
```

### Using RGI-BWT

We are going to start by picking just a small set of reads derived from an unknown plasmid to start. 
While RGI-BWT is pretty fast it would still take a while to run a full set of genomic (or metagenomic reads)!

```bash
mkdir rgi_bwt_results
cd rgi_bwt_results
```

Then we can run rgi bwt using the `bowtie2` short read mapper

```bash
rgi bwt --read_one ~/CourseData/IDE_data/module5/reads/unknown_plasmid/unknown_plasmid_1.fq.gz --read_two ~/CourseData/IDE_data/module5/reads/unknown_plasmid/unknown_plasmid_2.fq.gz --output_file unknown_plasmid_rgi_bwt --threads 2 --aligner bowtie2
```
RGI bwt is still in beta and produces a LOT of output files.

However, the file we are most interested in for now is `unknown_plasmid_rgi_bwt.gene_mapping_data.txt`

If you open this file (either by downloading it and opening it in a spreadsheet program) you can see AMR genes that have been identified in these reads.

* Which AMR genes were found?

Read-based analyses has advantages and disadvantages: 
* Higher sensitivity (we find as many AMR genes as possible) 
* Lower specificity (we are more likely to make mistakes when identifying AMR genes)
* Incomplete data (we are likely to find fragments of genes instead of whole genes, this can lead to confusion between similar genes)
* No genomic context (we don't know where a gene we detect comes from in the genome, is it associated with a plasmid?)

* Based on the results from this result which of these may be happening?


<a name="rgi"></a>
## AMR Detection from Assembled Contigs

That's the theory but how do read-based analyses really compare to assembly based analyses?

First, let's assemble the same set of reads we just analysed into contigs.
For this we are going to use `shovill` as it is a convenient way to get quick and reasonable good assemblies in one go.
Again, this is a relatively quick process in the grand scheme for a real set of bacterial genomic reads (5-10 million short-reads typically) but ~10-15 minutes is a lot of time for a tutorial.
Therefore, we are still going to use this small test set.

```bash
shovill --R1 ~/CourseData/IDE_data/module5/reads/unknown_plasmid/unknown_plasmid_1.fq.gz --R2 ~/CourseData/IDE_data/module5/reads/unknown_plasmid/unknown_plasmid_2.fq.gz --outdir unknown_plasmid --cpus 2
```

This should result in an assembled contig file to be created in `unknown_plasmid/contig.fa`.

RGI can then be run using these contigs as input.  This will involve RGI automatically detecting genes in the contigs (using `prodigal`) and comparing those genes against the CARD database using either BLAST or DIAMOND.
We'll use DIAMOND here as it faster than BLAST (at the cost of a very slightly decreased accuracy).

```bash
rgi main --input_sequence unknown_plasmid/contigs.fa --alignment_tool diamond --num_threads 2 --output_file unknown_plasmid_contig_rgi --clean
```

The output file we are interested in will be called `unknown_plasmid_contig_rgi.txt`



_INSERT MATERIAL ABOUT RUNNING ON ANOTHER SET OF CONTIGS FOR REAL_


<a name="hamronization"></a>
## hAMRonizing AMR Detection Results

Comparing the results between the two runs.


```bash
conda activate hamronization
hamronize rgi unknown_plasmid_contig_rgi.txt --analysis_software_version 5.2.0 --reference_database_version 3.1.3  --input_file_name unknown_plasmid --output unknown_plasmid_contig_hamronize.tsv
hamronize rgi unknown_plasmid_rgi_bwt.gene_mapping_data.txt --analysis_software_version 5.2.0 --reference_database_version 3.1.3  --input_file_name unknown_plasmid --output unknown_plasmid_reads_hamronize.tsv

hamronize summarize --summary_type interactive --output unknown_plasmid_comparison.html unknown_plasmid_contig_hamronize.tsv unknown_plasmid_reads_hamronize.tsv
```
Download `unknown_plasmid_comparison.html` and look at it



_INSERT MATERIAL ABOUT REPEATING THIS FOR A BUNCH OF THE CONTIG RESULTS_ 


<a name="heatmap"></a>
## Predicting Resistance Phenotype


_INSERT MATERIAL ABOUT HEATMAPS_

_LOOKING AT DRUG CLASS COMPARISON ESPECIALLY_

_WHY IS PSEUDOMONAS RESISTANT TO EVERYTHING_





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
