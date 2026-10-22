PMDD Gut Microbiome Analysis

Overview

This repository contains the analysis scripts used to identify associations between premenstrual dysphoric disorder (PMDD), menstrual phase, gut microbiome composition, predicted microbial functions, cytokines, and psychometric surveys.

The study includes longitudinal gut microbiome samples collected during the follicular and luteal menstrual phases.

Software and Dependencies

Analyses were performed using:

QIIME2: Version 2024.2

PICRUSt2: Version 2.6.2

R: Version 4.5.2

Major R packages include:

phyloseq

vegan

lme4

lmerTest

emmeans

ANCOMBC

ggplot2

dplyr

ComplexHeatmap

ape

Additional dependencies are specified in individual scripts.



PMDD-Microbiome/
	README.md
	
	QIIME2/
 		qiime2_pipeline.sh

	 PICRUSt2/
		picrust2_pipeline.sh

	 R_Analyses/
		alphadiversity.R
   		betadiversity.R
   		Taxa.R
   		Cytokine.R
   		Psychometric.R
		Steroid.R
   		PicrustCytokine.R
		KEGG.R
		
		
	 Figures/
	 	SuppFig1.R
		Figure1.R
		Figure2.R
		Figure3.R
		Figure4.R
		Figure5.R
		Figure6.R
		Figure7.R
		Figure8.R
    		

Analysis Workflow

1. Sequence Processing

16S rRNA gene sequencing data were processed using QIIME2. The pipeline includes quality control, feature-table generation, taxonomic classification, and alpha- and beta-diversity calculations. 

Commands are provided in QIIME2/qiime2_pipeline.sh.

2. Alpha and Beta Diversity

Alpha diversity was assessed using Shannon diversity, observed features, Faith's phylogenetic diversity, and Pielou's evenness.

Linear mixed-effects models were used to examine PMDD status, menstrual phase, and their interaction, accounting for repeated observations within participants, with the alphadiversity.R script.

Beta diversity was assessed using Bray–Curtis and unweighted UniFrac distances.

Associations were identified between beta diversity distances and PMDD, menstrual phase, and their interaction with the betadiversity.R script. 

3. Taxonomic Analysis

Genus-level differential abundance analyses were performed using ANCOM-BC2.

Associations with PMDD and menstrual phase were evaluated using the Taxa.R script.

4. Cytokine Analysis

Cytokine associations with alpha, beta, and taxonomic abundances were made with and without PMDD interaction in the Cytokine.R script.

5. Psychometric Analysis

Psychometric measures were evaluated using the Psychometric.R script.

6. Fecal Neuroactive Steroid Analysis

Fecal steroids were associated with alpha and beta diversity as well as taxonomic abundances in the Steroid.R script.

7. Functional Prediction

Predicted microbial functional profiles were generated using PICRUSt2 in the picrust2_pipeline.sh.

MetaCyc functional profiles were analyzed using ANCOM-BC2 against PMDD, phase, and cytokines in the PicrustCytokine.R script.

8. KEGG Pathways

KEGG pathways were identified and associated with the participant metadata in the KEGG.R script.

9. Figure Generation

Figure scripts generate the visualizations used in the manuscript, including diversity plots, principal coordinates analyses, taxonomic associations, and predicted functional pathway heatmaps.


Data Availability

Sequencing data: [repository accession or access restrictions]

Participant-level metadata: Available upon request to authors


Reproducibility- Input file paths may need to be adjusted according to the user's computing environment.


Citation: Not yet published


Contact- Liisa Hantsoo, lhantso1@jhmi.edu