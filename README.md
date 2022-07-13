***

<p align = "center">
  <img src = "https://user-images.githubusercontent.com/74168582/178574993-f3225e4e-63c5-4f67-8804-2aa4480ff24d.png" width="700" height="450">
</p>

***

</br>

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

- **iCn3D** is a powerful, web-based 3D viewer used for representing biomolecular structures. It facilitates visualization in 1D, 2D, and complex 3D while allowing synchronization of the selection across various structural displays. One of the important features of iCn3D is its ability to generate a shareable link that includes all the custom filters and user-provided labels. Anyone with this link can reproduce the same display of the structure. Annotations can either be directly extracted and mapped to the protein structure from various NCBI databases (such as dbSNP, ClinVar, Conserved Domain Database, and others), or users can submit their own annotations via custom tracks. Finally, it facilitates the display of both sequence-structure and structure-structure alignments with corresponding superposition.

- **Project goals**: In this project, we would like to retrieve SNPs from genomics data and show them in iCn3D with 1D/2D/3D representations on interaction networks. We will design a pipeline that automatically extracts SNPs from the variant call format (VCF) file and generates an iCn3D link annotated with variant effect predictions (SIFT, PolyPhen). Additionally, we would like to improve the visualization of SNPs mapped to the structure. This would allow us to study the deleterious effects of SNPs in the structural context.



## Materials and methods
# Materials:
- vcf file
- Compressed vcf file (vcf.gz)
- Indexed vcf file (vcf.gz.tbi)
- Gene ID
- 
**NOTE**: All files indicated above should be in the same directory

## Methodology
1. Input vcf.gz file and gene ID 
2. Extract Swissprot ID of the gene ID
- Connect to Uniprot database https://rest.uniprot.org/idmapping
- Convert from Ensembl ID to UniProtKB-Swiss-Prot
3. Extract gene information from Ensebl.org 
- Gene ID,start position,stop positions and chromosome number
4. Extract vcf file variants that match the given gene from the vcf.gz file 
- Start,start position, stop position  and chromosome number
5. Run VEP with the identified variants and capture sift and polyphen scores 
- use rest.ensembl.org
6. Generates the iCn3D link  based on the variants
- Extract the variants that are deleterious from sift and polyphen ducts and returns a combined string per prediction

# Test gene  ID 
- ENSG00000093072

# Test Files
- sample.vcf 
- t.vcf.gz  
- t.vcf.gz.tbi  

## Command Line Use:
```
python draft.py -v t.vcf.gz -g ENSG00000093072
```

## Results
Link: https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?afid=Q9NZK5&date=20220713&v=3.11.5&command=view

 Screenshort: 
![Screenshot 2022-07-13 at 18-00-57 Q9NZK5(AlphaFold) in iCn3D](https://user-images.githubusercontent.com/92297911/178769428-da5c1d34-fff0-4bc2-9d1b-e6937b75c5e0.png)

## Forthcoming features
- Filter SNPS from coding regions
- Highlighting the effect of each variant on the structure of the protein
- 
## Dependencies
- Python version >= 3
- Required modules:
  * pysam
- You can easily install this with conda (Ex: conda install -c bioconda pysam)
- Most other standard core modules should already be available on your system


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
- Collins, Francis S., Lisa D. Brooks, and Aravinda Chakravarti. "A DNA polymorphism discovery resource for research on human genetic variation." Genome research 8.12 (1998): 1229-1231.

- Gonzaga-Jauregui, Claudia, James R. Lupski, and Richard A. Gibbs. "Human genome sequencing in health and disease." Annual review of medicine 63 (2012): 35.

- McLaren, William, Laurent Gil, Sarah E. Hunt, Harpreet Singh Riat, Graham RS Ritchie, Anja Thormann, Paul Flicek, and Fiona Cunningham. "The ensembl variant effect predictor." Genome biology 17, no. 1 (2016): 1-14.

- Pagel, Kymberleigh A., Rick Kim, Kyle Moad, Ben Busby, Lily Zheng, Collin Tokheim, Michael Ryan, and Rachel Karchin. "Integrated informatics analysis of cancer-related variants." JCO clinical cancer informatics 4 (2020): 310-317.

- Wang, J., Youkharibache, P., Marchler-Bauer, A., Lanczycki, C., Zhang, D., Lu, S., ... & Ge, Y. (2022). iCn3D: From Web-Based 3D Viewer to Structural Analysis Tool in Batch Mode. Frontiers in Molecular Biosciences, 102.


## License
Licensed under MIT License - Copyright (c) 2022 hackathonismb (Refer LICENSE file for more details)
