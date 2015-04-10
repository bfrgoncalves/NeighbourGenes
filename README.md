# Neighboring genes

This program checks for the neighbor genes by syncing and downloading the desired files from ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/ and by performing BLASTp searches.

If files to be downloaded are compressed, its components are extracted and then used as reference for the BLASTp searches.

Python, Numpy, Biopython and HTSeq are required.

Launch the script
=================

python checkNeighbors.py -i `<Input file>` -b `<target bacteria>` -r `<dir for dwn files>` -f `<file type>` -u `<number genes upstream>` -d `<number genes downstream>` -o `<results dir>`

Arguments
=========

* `-i` Input fasta file -> Single or multiple fasta file to be used as query for the BLASTp searches
* `-t` Directory to look for -> Directory to look for inside ftp://ftp.ncbi.nlm.nih.gov/genomes/ . Example: `/Bacteria_DRAFT`
* `-b` Target bacteria name -> Specie to be searched on the ftp://ftp.ncbi.nlm.nih.gov/genomes/ and to be used as reference for the BLASTp searches
* `-r` Target directory for downloaded files -> Directory where all downloaded files will be placed
* `-f` File type to be downloaded -> file extension to be downloaded. Currently only .faa is supported
* `-u` Number of genes upstream -> Number of result genes upstream of the query sequence
* `-d` Number of genes downstream -> Number of result genes downstream of the query sequence
* `-o` Target directory for the results -> Directory where the files with the genes upstream and downstream of the query genes will be placed


Test Commands
============

On the cmd line at the application directory type:

`python checkNeighbors.py -i test.fasta -t /Bacteria -b Streptococcus_pneumoniae_670_6B -r downloads -f .faa -u 3 -d 3 -o results`

`python checkNeighbors.py -i test.fasta -t /Bacteria_DRAFT -b Streptococcus_pneumoniae_7286_06 -r downloads -f .faa -u 3 -d 3 -o results`