# Neighboring genes
This program checks for the neighbor genes by syncing and downloading files from ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/ and by peeforming BLASTp searches.

python checkNeighbours.py -i '<Input fasta file with proteic sequences>' -b '<Target bacteria>' -r '<Target directory for downloaded files>' -f '<Target file type>' -u '<Number of genes upstream>' -d '<Number of genes downstream>' -o '<Folder to write the results>'

Current the <Target file type> only supports .faa.
