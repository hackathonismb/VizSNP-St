## Visualizing NGS-generated SNPs in iCn3D

## Table of Contents
- [Introduction](#Introduction)
- [Methodology](#Methodology)
- [Results](#Results)
- [Dependencies](#Dependencies)
- [Team](#Team)
- [Acknowledgment](#Acknowledgment)
- [References](#References)
- [License](#License)

## Introduction

**The Problem**
Recent advancement in NGS technologies have led to the generation of a high volume of variant data.However,establishing the structural consequences of SNPs from the increasing volume of sequencing data remains a major challenge [1,2]. Most of the existing methods require highly skilled personnel to deduce meaningful structural and functional effects of the variants. This has led to development of tools to predict deleterious effects of SNPS. The SNP-analysis tools facilitate comprehensive evaluation of SNPs to non-structural biologists. For instance, VEP and OpenCRAVAT) are among sequence-based methods that are used to predict deleterious or missense SNPs [3,4]. Structural analysis of variants have been improved and annotation tools. Among these tools is the iCn3D structural viewer.However there is still need to improve not only visualization techniques but also intergration of SNPS and their effect on protein structure.

**What is  iCn3D?**
iCn3D is a web-based tool for structural analysis and annotation studies. It does not only analyze structures interactively but also in a batch mode. iCn3D structural viewer synchronizes  1D,2D and 3D structural displays. iCn3D can not only visualize SNPS from the existing databases,clinvar and dbsnp but also import custom stacks. For SNP analysis in iCn3D, a VCF can be initially read by a python script, then submitted to the VEP. An alpha folder is open in iCn3D or a list of PDB’s from Uniprot are used [2].

**Project goals**
In this project, we would like to have SNPs from VCF files to be opened directly in iCn3D without manual intervention. Additionally, we would like to improve on visualization of SNP effects on the structure by highlighting the SNPs. This will allow for quick click-through a provided list of SNPs, evaluate and make predictions on deleterious effects of SNPS in the structural context.



## Methodology


## Dependencies


## Team 
- Bonface Onyango
- Christian Cruz
- Manoj M Wagle
- Michael Sierk
- Pranavathiyani G

## Acknowledgment


## References
[1]. Gonzaga-Jauregui, Claudia, James R. Lupski, and Richard A. Gibbs. "Human genome sequencing in health and disease." Annual review of medicine 63 (2012): 35.

[2]. Wang, J., Youkharibache, P., Marchler-Bauer, A., Lanczycki, C., Zhang, D., Lu, S., ... & Ge, Y. (2022). iCn3D: From Web-Based 3D Viewer to Structural Analysis Tool in Batch Mode. Frontiers in Molecular Biosciences, 102.

[3]. McLaren, William, Laurent Gil, Sarah E. Hunt, Harpreet Singh Riat, Graham RS Ritchie, Anja Thormann, Paul Flicek, and Fiona Cunningham. "The ensembl variant effect predictor." Genome biology 17, no. 1 (2016): 1-14.

[4]. Pagel, Kymberleigh A., Rick Kim, Kyle Moad, Ben Busby, Lily Zheng, Collin Tokheim, Michael Ryan, and Rachel Karchin. "Integrated informatics analysis of cancer-related variants." JCO clinical cancer informatics 4 (2020): 310-317.


## License
Licensed under MIT License - Copyright (c) 2022 hackathonismb (Refer LICENSE file for more details)
