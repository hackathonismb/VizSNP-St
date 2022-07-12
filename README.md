## Visualizing NGS-generated SNPs in iCn3D


![Proteins](https://user-images.githubusercontent.com/92297911/178541105-d86b8198-7c54-445c-a6cd-6183dede708c.png)



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

Single nucleotide polymorphisms (SNPs) are characterized by a change in the single-nucleotide of the DNA sequence occurring at a specific position. Rapid advancements in next-generation sequencing (NGS) technologies have eased the identification of these point mutations. They account for approximately 90% of genetic variations in humans and are found to be associated with several diseases. Tools such as VEP and OpenCRAVAT have been developed to predict the deleterious effects of these variants. However, the consequences of SNPs on protein structure and function are poorly understood. Hence, there is a great need to integrate this data with the available structural data. Mapping of SNPs to the protein structure and its visualization can significantly help us in studying various disease-causing mutations.

- **iCn3D** is a web-based tool for structural analysis and annotation studies. It does not only analyze structures interactively but also in a batch mode. iCn3D structural viewer synchronizes  1D,2D and 3D structural displays. iCn3D can not only visualize SNPS from the existing databases,clinvar and dbsnp but also import custom stacks. For SNP analysis in iCn3D, a VCF can be initially read by a python script, then submitted to the VEP. An alpha folder is open in iCn3D or a list of PDB’s from Uniprot are used [2].

- **Project goals**: In this project, we would like to have SNPs from VCF files to be opened directly in iCn3D without manual intervention. Additionally, we would like to improve on visualization of SNP effects on the structure by highlighting the SNPs. This will allow for quick click-through a provided list of SNPs, evaluate and make predictions on deleterious effects of SNPS in the structural context.



## Methodology


## Dependencies


## Team 
- Bonface Onyango
- Christian Cruz
- Manoj M Wagle
- Michael Sierk
- Pranavathiyani G
-Brenda Kiage
## Acknowledgment


## References
[1] Gonzaga-Jauregui, Claudia, James R. Lupski, and Richard A. Gibbs. "Human genome sequencing in health and disease." Annual review of medicine 63 (2012): 35.

[2] Wang, J., Youkharibache, P., Marchler-Bauer, A., Lanczycki, C., Zhang, D., Lu, S., ... & Ge, Y. (2022). iCn3D: From Web-Based 3D Viewer to Structural Analysis Tool in Batch Mode. Frontiers in Molecular Biosciences, 102.

[3] McLaren, William, Laurent Gil, Sarah E. Hunt, Harpreet Singh Riat, Graham RS Ritchie, Anja Thormann, Paul Flicek, and Fiona Cunningham. "The ensembl variant effect predictor." Genome biology 17, no. 1 (2016): 1-14.

[4] Pagel, Kymberleigh A., Rick Kim, Kyle Moad, Ben Busby, Lily Zheng, Collin Tokheim, Michael Ryan, and Rachel Karchin. "Integrated informatics analysis of cancer-related variants." JCO clinical cancer informatics 4 (2020): 310-317.


## License
Licensed under MIT License - Copyright (c) 2022 hackathonismb (Refer LICENSE file for more details)
