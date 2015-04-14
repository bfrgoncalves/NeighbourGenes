#!/usr/bin/python
import ftplib
import os
from datetime import datetime
import sys
import argparse
import shutil
import urllib
import tarfile
from os import listdir
from os.path import isfile, join

from Bio import SeqIO
import HTSeq
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline

from CommonFastaFunctions import runBlastParser

def main():
	parser = argparse.ArgumentParser(description="This program checks for the neighbour genes of a <multi fasta file> by syncing and downloading files from a <local target dir> of ftp://ftp.ncbi.nlm.nih.gov/genomes/ and by peforming BLASTp searches.")
	parser.add_argument('-i', nargs='?', type=str, help="Input fasta file with proteic sequences", required=True)
	parser.add_argument('-t', nargs='?', type=str, help="Target URL inside ftp://ftp.ncbi.nlm.nih.gov/genomes/", required=True)
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
	targetURL = args.t

	downloadedFiles = SearchAndDownload(target_bug,target_dir,file_type, targetURL)

	DataBaseFolder = os.path.join(target_dir,'Databases')
	if not os.path.isdir(DataBaseFolder):
		os.makedirs(DataBaseFolder)
	if not os.path.isdir(PathToWrite):
		os.makedirs(PathToWrite)
	
	DatabaseName = os.path.join(DataBaseFolder,target_bug)

	queryNames, querySequences, sequenceLengths = ReadFASTAFile(queryFile)

	countFiles = 0
	countNames = 0
	for Names in queryNames:
		countFiles = 0
		queryFileToUse = Create_FASTAquery(queryFile, Names, querySequences[countNames])
		
		PathFolderName = os.path.join(PathToWrite, Names)
		if not os.path.isdir(PathFolderName):
			os.makedirs(PathFolderName)

		for fl in downloadedFiles:
			print
			print '-------------------------------------------'
			dirName = fl.split('/')[len(fl.split('/'))-2]
			fileName = fl.split('/')[len(fl.split('/'))-1].split('.')[0]
			PathW = os.path.join(PathFolderName,dirName)
			if not os.path.isdir(PathW):
				os.makedirs(PathW)
			PathW = os.path.join(PathW,fileName+'result.faa')
			NameGenes, SequenceGenes, lenRefSeq = ReadFASTAFile(fl)
			Create_Blastdb(fl, True, DatabaseName)
			print 'Performing BLASTp search'
			matchGene = BLASTp(queryFileToUse, DatabaseName, DataBaseFolder, queryNames, sequenceLengths)
			ToWrite = GetNeighbours(NameGenes,SequenceGenes,matchGene, Names, querySequences[countNames], numberUpstream, numberDownstream)
			WriteResults(ToWrite, PathW, target_dir+'/'+fileName+file_type)

			countFiles+=1
		countNames +=1
		os.remove(queryFileToUse)

	shutil.rmtree(DataBaseFolder)
	print '*************************************'
	print '---------------DONE!-----------------'




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

def BLASTp(queryFile, dbName, blast_out_path, queryNames, sequenceLengths):
	blast_out_file = os.path.join(blast_out_path,'blastOut.xml')
	cline = NcbiblastpCommandline(query=queryFile, db=dbName, out=blast_out_file, outfmt=5)
	blast_records = runBlastParser(cline,blast_out_file, "")
	matchGene = ''
	score = -1
	for blast_record in blast_records:
		queryGeneIndex = queryNames.index(blast_record.query.strip('|'))
		querySequenceLength = sequenceLengths[queryGeneIndex]
		for alignment in blast_record.alignments:
			for match in alignment.hsps:
				identity_length_ratio = float(match.identities)/float(querySequenceLength)
				if identity_length_ratio >= 0.8:
					if score < match.score:
						matchGene = alignment.hit_def
						score = match.score

	return matchGene


def GetNeighbours(arrayOfGenes, arrayOfSequences, targetGene, queryName, querySequence, numberUpstream, numberDownstream):
	try:
		targetIndex = arrayOfGenes.index(targetGene)
		ToWrite = []
		for x in range(targetIndex-numberUpstream, targetIndex):
			try:
				ToWrite.append('>'+arrayOfGenes[x]+"\n"+arrayOfSequences[x]+'\n')
			except IndexError:
				print "There are no results "+ str(x) + "levels upstream " + targetGene

		#ToWrite.append('>'+queryName+'\n'+querySequence+'\n') # INSERT QUERY ON THE RESULTS
		
		for x in range(targetIndex, targetIndex+numberDownstream+1):
			try:
				ToWrite.append('>'+arrayOfGenes[x]+"\n"+arrayOfSequences[x]+'\n')
			except IndexError:
				print "There are no results" + str(x) + "levels downstream " + targetGene
	except ValueError:
		ToWrite = ["There were no matches for the query sequence"]

	return ToWrite

def WriteResults(ToWrite, PathToWrite, referenceFile):
	if 'There were no matches' not in ToWrite[0]:
		print 'Writing the results'
		with open(PathToWrite, "w") as r:
			for result in ToWrite:
				r.write(result)
	else:
		print "There were no matches for the query sequence at " + referenceFile

def ReadFASTAFile(FASTAfile):
	NameSeq = []
	Sequence = []
	sequenceLength = []
	handle = open(FASTAfile, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	for record in records:
		NameSeq.append(str(record.description))
		Sequence.append(str(record.seq))
		sequenceLength.append(len(str(record.seq)))
	return NameSeq,Sequence,sequenceLength

def SearchAndDownload(target_bug,target_dir,file_type,dirToUse):

	ftp = 'ftp.ncbi.nih.gov'
	subURL = '/genomes' + dirToUse
	print "Connecting to ftp.ncbi.nih.gov..."	
	f=ftplib.FTP(ftp)
	f.login()
	f.cwd(subURL)
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
						mainURL = 'ftp://'+ftp
						currentURL = os.path.join(subURL,item,fi)

						if '.tgz' in currentURL:
							ftpstream = urllib.urlopen(mainURL+currentURL)
							print 'Extracting ' + fi
							tar = tarfile.open(fileobj=ftpstream, mode="r|gz")
							tardir = os.path.join(target_dir,item)
							tar.extractall(tardir)
							onlyfiles = [ h for h in listdir(tardir) if isfile(join(tardir,h)) ]
							for files in onlyfiles:
								local_file = os.path.join(target_dir,item,files)
								if '.faa' not in files:
									os.remove(local_file)
								else:
									ct+=1
									arrayOfFiles.append(local_file)

						else:
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