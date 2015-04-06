# Neighboring genes

This program checks for the neighbor genes by syncing and downloading the desired files from ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/ and by performing BLASTp searches.

Python is required.

Launch the script
=================

`python checkNeighbours.py -i <Input file> -b <target bacteria> -r <dir for dwn files> -f <file type> -u <number genes upstream> -d <number genes downstream> -o <results dir>`

Arguments
=========

* -i `Input fasta file -> Single or multiple fasta file to be used as query for the BLASTp searches`
* -b `Target bacteria name -> Specie to be searched on the ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/ and to be used as reference for the BLASTp searches`
* -r `Target directory for downloaded files -> Directory were all files downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/ will be placed`
* -f `File type to be downloaded -> file extension to be downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/. Currently only .faa is supported`
* -u `Number of genes upstream -> Number of result genes upstream of the query sequence`
* -d `Number of genes downstream -> Number of result genes downstream of the query sequence`
* -o `Target directory for the results -> Directory where the files with the genes upstream and downstream of the query genes will be placed`

Current the (Target file type) only supports .faa.