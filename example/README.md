## Running the example
You need a VCF file that has been compressed with bgzip, and then run through tabix to create the index file (.tbi)
You also need to know the Ensembl Gene ID that you are looking at (in this case use ENSG00000141867).

Run the example like so:
```
python SNP2iCn3D.py -v t3.vcf.gz -g ENSG00000141867
```

The output should look like this:
```
UniProt Primary Accession: O60885

Here is your iCn3D link:
https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?afid=O60885&date=20220713&v=3.12.7&command=view annotations; set annotation cdd; set view detailed view;add track | chainid O60885_A | title SIFT_predict | text 517 P;scap interaction O60885_A_517_P
```

The structure will be opened automatically in iCn3D.