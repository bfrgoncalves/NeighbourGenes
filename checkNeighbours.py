#!/usr/bin/python
import ftplib
import os
from datetime import datetime
import sys
import argparse
import shutil

from Bio import SeqIO
from BCBio import GFF
import HTSeq
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline

from CommonFastaFunctions import runBlastParser

def main():
	parser = argparse.ArgumentParser(description="This program checks for the neighbour genes of a <multi fasta file> by syncing and downloading files from a <local target dir> of ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/ and by peforming BLASTp searches.")
	parser.add_argument('-i', nargs='?', type=str, help="Input fasta file with proteic sequences", required=True)
	parser.add_argument('-b', nargs='?', type=str, help="Target bacteria", required=True)
	parser.add_argument('-r', nargs='?', type=str, help="Target directory for downloaded files", required=True)
	parser.add_argument('-f', nargs='?', type=str, help='Target file type. Currently only .faa can be used.', required=True)
	parser.add_argument('-u', nargs='?', type=int, help='Number of genes upstream.', required=True)
	parser.add_argument('-d', nargs='?', type=int, help='Number of genes downstream.', required=True)
	parser.add_argument('-o', nargs='?', type=str, help='Folder to write the results.', required=True)

	args = parser.parse_args()

	queryFile = args.i
	target_bug = args.b
	target_dir = args.r
	file_type =args.f
	numberUpstream =args.u
	numberDownstream =args.d
	PathToWrite = args.o

	downloadedFiles = SearchAndDownload(target_bug,target_dir,file_type)

	DataBaseFolder = os.path.join(target_dir,'Databases')
	if not os.path.isdir(DataBaseFolder):
		os.makedirs(DataBaseFolder)
	if not os.path.isdir(PathToWrite):
		os.makedirs(PathToWrite)
	
	DatabaseName = os.path.join(DataBaseFolder,target_bug)

	queryNames, querySequences = ReadFASTAFile(queryFile)

	countFiles = 0
	countNames = 0
	for Names in queryNames:
		countFiles = 0
		queryFileToUse = Create_FASTAquery(queryFile, Names, querySequences[countNames])
		
		PathFolderName = os.path.join(PathToWrite, Names)
		if not os.path.isdir(PathFolderName):
			os.makedirs(PathFolderName)

		for fl in downloadedFiles:
			print '-------------------------------------------'
			print 'Search target gene on ' + fl
			dirName = fl.split('/')[len(fl.split('/'))-2]
			fileName = fl.split('/')[len(fl.split('/'))-1].split('.')[0]
			PathW = os.path.join(PathFolderName,dirName)
			if not os.path.isdir(PathW):
				os.makedirs(PathW)
			PathW = os.path.join(PathW,fileName+'result.faa')
			NameGenes, SequenceGenes = ReadFASTAFile(fl)
			Create_Blastdb(fl, True, DatabaseName)
			print 'Performing BLASTp search'
			matchGene = BLASTp(queryFileToUse, DatabaseName, DataBaseFolder)
			ToWrite = GetNeighbours(NameGenes,SequenceGenes,matchGene, Names, querySequences[countNames], numberUpstream, numberDownstream)
			print 'Writing the results'
			WriteResults(ToWrite, PathW)

			countFiles+=1
		countNames +=1
		os.remove(queryFileToUse)

	shutil.rmtree(DataBaseFolder)
	print 'DONE!'




#######################################################################################################################################################################

def Create_FASTAquery(queryFile, geneName,geneSequence):
	newName = queryFile.replace('.','prov.')
	of_handle=open(newName, 'w')
	of_handle.write(">"  + geneName + "|\n" + geneSequence +"\n")
	of_handle.close()
	return newName

#Create a Blast DB:
def Create_Blastdb( questionDB, dbtypeProt, dbName ):
    isProt=dbtypeProt

    if not os.path.isfile(dbName + ".nin") and not os.path.isfile(dbName + ".nhr") and not os.path.isfile(dbName + ".nsq"):
        
        if not isProt:
            os.system( "makeblastdb -in " + questionDB + " -out " + dbName + " -dbtype nucl -logfile " + dbName + "_blast.log" )
        else:
            os.system( "makeblastdb -in " + questionDB + " -out " + dbName + " -dbtype prot -logfile " + dbName + "_blast.log" )

    elif overwrite:
        if not isProt:
            os.system( "makeblastdb -in " + questionDB + " -out " + dbName + " -dbtype nucl -logfile " + dbName + "_blast.log" )
        else:
            os.system( "makeblastdb -in " + questionDB + " -out " + dbName + " -dbtype prot -logfile " + dbName + "_blast.log" )

    else:
        print "BLAST DB files found. Using existing DBs.."  
    return( dbName )

def BLASTp(queryFile, dbName, blast_out_path):
	blast_out_file = os.path.join(blast_out_path,'blastOut.xml')
	cline = NcbiblastpCommandline(query=queryFile, db=dbName, out=blast_out_file, outfmt=5)
	blast_records = runBlastParser(cline,blast_out_file, "")
	matchGene = ''
	score = -1
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for match in alignment.hsps:
				if score < match.score:
					matchGene = alignment.hit_def
					score = match.score

	return matchGene


def GetNeighbours(arrayOfGenes, arrayOfSequences, targetGene, queryName, querySequence, numberUpstream, numberDownstream):
	targetIndex = arrayOfGenes.index(targetGene)
	ToWrite = []
	for x in range(targetIndex-numberUpstream, targetIndex):
		try:
			ToWrite.append('>'+arrayOfGenes[x]+"\n"+arrayOfSequences[x]+'\n')
		except IndexError:
			print "There are no genes Upstream " + arrayOfGenes[x]

	ToWrite.append('>'+queryName+'\n'+querySequence)
	
	for x in range(targetIndex, targetIndex+numberDownstream+1):
		try:
			ToWrite.append('>'+arrayOfGenes[x]+"\n"+arrayOfSequences[x]+'\n')
		except IndexError:
			print "There are no genes Downstream " + arrayOfGenes[x]

	return ToWrite

def WriteResults(ToWrite, PathToWrite):
	with open(PathToWrite, "w") as r:
		for result in ToWrite:
			r.write(result)

def ReadFASTAFile(FASTAfile):
	NameSeq = []
	Sequence = []
	handle = open(FASTAfile, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	for record in records:
		NameSeq.append(str(record.description))
		Sequence.append(str(record.seq))
	return NameSeq,Sequence

def SearchAndDownload(target_bug,target_dir,file_type):

	print "Connecting to ftp.ncbi.nih.gov..."	
	f=ftplib.FTP('ftp.ncbi.nih.gov')
	f.login()
	f.cwd('/genomes/Bacteria/')
	listing=[]
	dirs=f.nlst();
	print "Connected and Dir list retrieved."

	arrayOfFiles = []
	targetNames = []

	print "Searching for :"+ target_bug
	ct=0;
	for item in dirs:
		if item.find(target_bug)>-1:
			print
			print "----------------------------------------------"
			print "Dir: " + item
			targetNames.append(item)
			#create the dir
			if not os.path.isdir(os.path.join(target_dir,item)):
				print "Dir not found. Creating it..."
				os.makedirs(os.path.join(target_dir,item))
			#1) change the dir
			f.cwd(item)
			#2) get  files from file_type in dir
			try:
				files=f.nlst()
				for fi in files:
					if file_type in fi:
						local_file = os.path.join(target_dir,item,fi)
						arrayOfFiles.append(local_file)
						if os.path.isfile(local_file):
							print "Dir:" + item	
							print "File " + local_file + " already exists."
							#get remote modification time			
							mt = f.sendcmd('MDTM '+ fi)
							#converting to timestamp
							nt = datetime.strptime(mt[4:], "%Y%m%d%H%M%S").strftime("%s")

							if int(nt)==int(os.stat(local_file).st_mtime):
								print fi +" not modified. Download skipped"
							else:
								print "New version of "+fi
								ct+=1
								DownloadAndSetTimestamp(local_file,fi,nt,f)
								print "NV Local M timestamp : " + str(os.stat(local_file).st_mtime)
								print "NV Local A timestamp : " + str(os.stat(local_file).st_atime)

						else:
							print "New file: "+fi
							ct+=1
							mt = f.sendcmd('MDTM '+ fi)
							#converting to timestamp
							nt = datetime.strptime(mt[4:], "%Y%m%d%H%M%S").strftime("%s")
							DownloadAndSetTimestamp(local_file,fi,nt,f)
			except ftplib.error_temp,  resp:
				if str(resp) == "450 No files found":
					print "No "+ file_type +" files in this directory. Skipping"
			f.cwd('..')
	f.quit()
	print "# of "+target_bug+" new files found and downloaded: " + str(ct)
	return arrayOfFiles


def DownloadAndSetTimestamp(local_file,fi,nt,f):
	lf=open(local_file,'wb')
	f.retrbinary("RETR " + fi, lf.write, 8*1024)
	lf.close()
	print fi + " downloaded!"

	#set the modification time the same as server for future comparison
	os.utime(local_file,( int(nt) , int(nt) ))



if __name__ == "__main__":
    main()