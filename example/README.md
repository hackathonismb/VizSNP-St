## Running the example
You need a VCF file that has been compressed with bgzip, and then run through tabix to create the index file (.tbi)
You also need to know the Ensembl Gene ID that you are looking at (in this case use ENSG00000141867).

Run the example like so:
```
python SNP2iCn3D.py -v test.vcf.gz 
```

The output will be in iCn3D_links.html, and a .csv file with the same name as the input .vcf.gz file 
