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

- **Single nucleotide polymorphisms (SNPs)** are characterized by a change in the single-nucleotide of the DNA sequence occurring at a specific position. Rapid advancements in next-generation sequencing (NGS) technologies have eased the identification of these point mutations. They account for approximately 90% of genetic variations in humans and are found to be associated with several diseases. Tools such as VEP and OpenCRAVAT have been developed to predict the deleterious effects of these variants. However, the consequences of SNPs on protein structure and function are poorly understood. Hence, there is a great need to integrate this data with the available structural data. Mapping of SNPs to the protein structure and its visualization can significantly help us in studying various disease-causing mutations.

- **iCn3D** is a powerful, web-based 3D viewer used for representing biomolecular structures. It facilitates visualization in 1D, 2D, and complex 3D while allowing synchronization of the selection across various structural displays. One of the important features of iCn3D is its ability to generate a shareable link that includes all the custom filters and user-provided labels. Anyone with this link can reproduce the same display of the structure. Annotations can either be directly extracted and mapped to the protein structure from the various NCBI databases (such as dbSNP, ClinVar, Conserved Domain Database, and others), or users can submit their own annotations via custom tracks. Finally, it facilitates the display of both sequence-structure and structure-structure alignments with corresponding superposition.

- **Project goals**: In this project, we would like to have SNPs from VCF files to be opened directly in iCn3D without manual intervention. Additionally, we would like to improve on visualization of SNP effects on the structure by highlighting the SNPs. This will allow for quick click-through a provided list of SNPs, evaluate and make predictions on deleterious effects of SNPS in the structural context.



## Methodology


## Dependencies


## Team 
- Bonface Onyango
- Brenda Kiage
- Christian Cruz
- Manoj M Wagle
- Michael Sierk
- Pranavathiyani G

## Acknowledgment
- We would like to thank the **International Society for Computational Biology/Intelligent Systems for Molecular Biology (ISCB/ISMB)** and the **National Center for Biotechnology Information (NCBI)** for their support and for providing all the required computational resources during the codeathon.


## References
- Gonzaga-Jauregui, Claudia, James R. Lupski, and Richard A. Gibbs. "Human genome sequencing in health and disease." Annual review of medicine 63 (2012): 35.

- McLaren, William, Laurent Gil, Sarah E. Hunt, Harpreet Singh Riat, Graham RS Ritchie, Anja Thormann, Paul Flicek, and Fiona Cunningham. "The ensembl variant effect predictor." Genome biology 17, no. 1 (2016): 1-14.

- Pagel, Kymberleigh A., Rick Kim, Kyle Moad, Ben Busby, Lily Zheng, Collin Tokheim, Michael Ryan, and Rachel Karchin. "Integrated informatics analysis of cancer-related variants." JCO clinical cancer informatics 4 (2020): 310-317.

- Wang, J., Youkharibache, P., Marchler-Bauer, A., Lanczycki, C., Zhang, D., Lu, S., ... & Ge, Y. (2022). iCn3D: From Web-Based 3D Viewer to Structural Analysis Tool in Batch Mode. Frontiers in Molecular Biosciences, 102.


## License
Licensed under MIT License - Copyright (c) 2022 hackathonismb (Refer LICENSE file for more details)
