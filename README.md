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
- [Future prospects](#Future-prospects)
- [Dependencies](#Dependencies)
- [Team](#Team)
- [Acknowledgment](#Acknowledgment)
- [References](#References)
- [License](#License)

## Introduction

- **Single nucleotide polymorphisms (SNPs)** are characterized by a change in the single-nucleotide of the DNA sequence occurring at a specific position. Rapid advancements in next-generation sequencing (NGS) technologies have eased the identification of these point mutations. They account for approximately 90% of genetic variations in humans and are found to be associated with several diseases. Tools such as VEP and OpenCRAVAT have been developed to predict the deleterious effects of these variants. However, the consequences of SNPs on protein structure and function are poorly understood. Hence, there is a great need to integrate this data with the available structural data. Mapping of SNPs to the protein structure and its visualization can significantly help us in studying various disease-causing mutations.

- **iCn3D** is a powerful, web-based 3D viewer used for representing biomolecular structures. It facilitates visualization in 1D, 2D, and complex 3D while allowing synchronization of the selection across various structural displays. One of the important features of iCn3D is its ability to generate a shareable link that includes all the custom filters and user-provided labels. Anyone with this link can reproduce the same display of the structure. Annotations can either be directly extracted and mapped to the protein structure from various NCBI databases (such as dbSNP, ClinVar, Conserved Domain Database, and others), or users can submit their own annotations via custom tracks. Finally, it facilitates the display of both sequence-structure and structure-structure alignments with corresponding superposition.

- **Project goals**: In this project, we would like to retrieve SNPs from genomics data and show them in iCn3D with 1D/2D/3D representations on interaction networks. We will design a pipeline that automatically extracts SNPs from the variant call format (VCF) file and generates an iCn3D link annotated with variant effect predictions (SIFT, PolyPhen). Additionally, we would like to improve the visualization of SNPs mapped to the structure. This would allow us to study the deleterious effects of SNPs in the structural context.


## Methodology 
**1] Input Gene ID, annotation, and mapping**
- The script takes as input an Ensembl gene ID and a VCF file.
- Annotation: The following information is retrieved for the gene ID using the [Ensembl Rest API] (https://rest.ensembl.org):
    * Corresponding Ensembl protein ID
    * Chromosomal coordinates
- Mapping: 
    * The input gene ID is also mapped to a UniProt primary accession (UPrimAC) using the UniProtKB/Swiss-Prot database.
    * The UPrimAC is used to retrieve the protein structure.
    * If a PDB structure is available, then structure with the best resolution is chosen. When a PDB structure is unavailable, an Alphafold structure is retrieved instead.

**2] Extract VCF file variants**
- Variants that match the input gene are extracted from the VCF file

**3] Predict the functional effects of variants**
- The deleterious effect of the identified variants is predicted using the Ensembl Variant Effect Predictor (VEP) Rest API.
- Currently, we utilize two pathogenicity prediction scores: SIFT and PolyPhen.

**4] Generate the iCn3D path based on variants of interest**
- Finally, an `iCn3D link` is generated using the protein structure retrieved earlier (see step 1) with the added variant annotation information (those predicted to be `deleterious` by SIFT and PolyPhen).
- Automatically open the link in your favorite web browser.

</br>
<img width="1391" alt="Flowchart_VizSNP-St" src="https://user-images.githubusercontent.com/74168582/178969212-40c68667-9b2c-46cc-bfe1-2c7a6f38de27.png">


## Results

- For running the pipeline, you need a VCF file (gzipped and tabix) from which the variants are to be extracted. 
- Here we tested our pipeline on a sample VCF file (available in the `example` directory)

  ```
  python3 VizSNPSt.py -g ENSG00000141867 -v example.vcf.gz 
  ```

- The output should look something like this:

  ```
  UniProt Primary Accession: O60885

  Here is your iCn3D link:
  https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?afid=O60885&date=20220713&v=3.12.7&command=view annotations; set annotation cdd; set view detailed view;add track | chainid O60885_A | title SIFT_predict | text 517 P;scap interaction O60885_A_517_P
  ```

- The structure with the SNPs will be opened automatically in iCn3D.
- The figures and animation below shows the output in iCn3D

</br>

<p align = "center">
  <img width="1789" alt="Main" src="https://user-images.githubusercontent.com/74168582/180653738-b4bf63d6-a42f-4058-9046-b78e1687919d.png">
 </p>

</br>

<p align="center">
  <img src="https://user-images.githubusercontent.com/74168582/180654007-9da19c54-5793-4ab9-a3b3-1800e16fd942.gif" width="600" height="450">
</p>

## Future prospects:
- Remove requirement for knowing the Ensembl Gene ID (just submit a VCF)
- Perform more sophisticated filtering on the VCF file. (Right now it just selects deleterious mutations.)
- Deal with multiple SNPs, either from a single gene or multiple genes.

## Dependencies
- Python version >= 3
- Required modules:
  * pysam
- You can easily install this with conda (Ex: conda install -c bioconda pysam)
- Most other standard core modules should already be available on your system

## Team 
- [Bonface Onyango](https://github.com/bonfaceonyango)
- [Manoj M Wagle](https://github.com/manojmw)
- [Michael Sierk](https://github.com/msierk)
- [Pranavathiyani G](https://github.com/pranavathiyani)

## Acknowledgment
- We would like to thank the **International Society for Computational Biology/Intelligent Systems for Molecular Biology (ISCB/ISMB)** and the **National Center for Biotechnology Information (NCBI)** for their support and for providing all the required computational resources during the codeathon.
- We would also like to thank Shashi Ranayake from the Center for Biomedical Informatics and Information Technology, Computational Genomics and Bioinformatics Branch at the National Cancer Institute for drafting the original version of the script. Subsequent improvements and new features were added by Manoj M Wagle and Michael Sierk.


## References
- Collins, Francis S., Lisa D. Brooks, and Aravinda Chakravarti. "A DNA polymorphism discovery resource for research on human genetic variation." Genome research 8.12 (1998): 1229-1231.

- Gonzaga-Jauregui, Claudia, James R. Lupski, and Richard A. Gibbs. "Human genome sequencing in health and disease." Annual review of medicine 63 (2012): 35.

- McLaren, William, Laurent Gil, Sarah E. Hunt, Harpreet Singh Riat, Graham RS Ritchie, Anja Thormann, Paul Flicek, and Fiona Cunningham. "The ensembl variant effect predictor." Genome biology 17, no. 1 (2016): 1-14.

- Pagel, Kymberleigh A., Rick Kim, Kyle Moad, Ben Busby, Lily Zheng, Collin Tokheim, Michael Ryan, and Rachel Karchin. "Integrated informatics analysis of cancer-related variants." JCO clinical cancer informatics 4 (2020): 310-317.

- Wang, J., Youkharibache, P., Marchler-Bauer, A., Lanczycki, C., Zhang, D., Lu, S., ... & Ge, Y. (2022). iCn3D: From Web-Based 3D Viewer to Structural Analysis Tool in Batch Mode. Frontiers in Molecular Biosciences, 102.


## License
Licensed under MIT License - Copyright (c) 2022 hackathonismb (Refer LICENSE file for more details)
