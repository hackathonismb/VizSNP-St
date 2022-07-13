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

Single nucleotide polymorphisms (SNPs) are characterized by a change in the single-nucleotide of the DNA sequence occurring at a specific position. Rapid advancements in next-generation sequencing (NGS) technologies have eased the identification of these point mutations. They account for approximately 90% of genetic variations in humans and are found to be associated with several diseases. Tools such as VEP and OpenCRAVAT have been developed to predict the deleterious effects of these variants. However, the consequences of SNPs on protein structure and function are poorly understood. Hence, there is a great need to integrate this data with the available structural data. Mapping and visualization of SNPs to the protein structure can significantly help us in studying various disease-causing mutations.

- **iCn3D** is a web-based tool for structural analysis and annotation studies. It does not only analyze structures interactively but also in a batch mode. iCn3D structural viewer synchronizes  1D,2D and 3D structural displays. iCn3D can not only visualize SNPS from the existing databases,clinvar and dbsnp but also import custom stacks. For SNP analysis in iCn3D, a VCF can be initially read by a python script, then submitted to the VEP. An alpha folder is open in iCn3D or a list of PDB’s from Uniprot are used [2].

- **Project goals**: In this project, we would like to have SNPs from VCF files to be opened directly in iCn3D without manual intervention. Additionally, we would like to improve on visualization of SNP effects on the structure by highlighting the SNPs. This will allow for quick click-through a provided list of SNPs, evaluate and make predictions on deleterious effects of SNPS in the structural context.



## Materials and methods
# Materials:
- vcf file
- Compressed vcf file (vcf.gz)
- Indexed vcf file (vcf.gz.tbi)
- Gene ID
- 
**NOTE**: All files indicated above should be in the same directory

## Methodology
1. Input vcf.gz file and gene ID  Compressed and indexed VCF files
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





## Dependencies


## Team 
- Bonface Onyango
- Christian Cruz
- Manoj M Wagle
- Michael Sierk
- Pranavathiyani G

## Acknowledgment


## References
[1] Gonzaga-Jauregui, Claudia, James R. Lupski, and Richard A. Gibbs. "Human genome sequencing in health and disease." Annual review of medicine 63 (2012): 35.

[2] Wang, J., Youkharibache, P., Marchler-Bauer, A., Lanczycki, C., Zhang, D., Lu, S., ... & Ge, Y. (2022). iCn3D: From Web-Based 3D Viewer to Structural Analysis Tool in Batch Mode. Frontiers in Molecular Biosciences, 102.

[3] McLaren, William, Laurent Gil, Sarah E. Hunt, Harpreet Singh Riat, Graham RS Ritchie, Anja Thormann, Paul Flicek, and Fiona Cunningham. "The ensembl variant effect predictor." Genome biology 17, no. 1 (2016): 1-14.

[4] Pagel, Kymberleigh A., Rick Kim, Kyle Moad, Ben Busby, Lily Zheng, Collin Tokheim, Michael Ryan, and Rachel Karchin. "Integrated informatics analysis of cancer-related variants." JCO clinical cancer informatics 4 (2020): 310-317.


## License
Licensed under MIT License - Copyright (c) 2022 hackathonismb (Refer LICENSE file for more details)
