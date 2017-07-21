## This code reads Olivers ORFs and transcript fasta files, and protein and peptide identification files.
## Filters ORFs and transcript fasta files using fastaFileFiltering.py. It also filters the peptide csv.

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import re
import argparse
import pandas as pd

## Read Identified Proteins

def readIdentifiedProteinPeptide(filename,sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj
    
def pepThresholding(prots, pepTh, protRevStr, protContStr):
    print("1. prots dim:"+str(prots.shape))
    prots=prots[~prots['protein accession'].str.contains(protRevStr)]
    print("2. prots dim:"+str(prots.shape))
    prots=prots[~prots['protein accession'].str.contains(protContStr)]
    print("3. prots dim:"+str(prots.shape))
    prots=prots[prots['distinct peptide sequences']>pepTh]
    print("4. prots dim:"+str(prots.shape))
    return prots
    
def identifiedProteinList(prots, columnName):
    identifiedProts=prots[columnName].tolist()
    return identifiedProts
    
def filterFasta(fastaFile, overlapping_ids, outFile):
	fastahandle = open(fastaFile, "rU")
	outputHandle = open(outFile, "w")
	#print(overlapping_ids)
	overlapping_ids[:] = [line.strip() for line in overlapping_ids]
	overlapping_ids[:] = [line.replace('"','') for line in overlapping_ids]
	#print("oids:"+str(len(overlapping_ids)))
	#print(overlapping_ids[0])
	records = list(SeqIO.parse(fastahandle, "fasta"))
	#print("Total=")
	#print(len(records))
	count1=0
	count2=0
	for record in records: # SeqIO.parse(handle, "fasta")
	    #print(record.description)
	    #print("TESTTEST")
	    if record.description in overlapping_ids:
	        overlapping_ids.remove(record.description)
	        count1=count1+1
	        outputHandle.write(">"+record.description+"\n")
	        outputHandle.write(str(record.seq)+"\n")
	        #continue
	        #print("found:"+record.description)
	fastahandle.close() 
	outputHandle.close()

def extractProtIdsFromDescription(prots, sep):
    ##prots is a list of strings, generally the description column of protein identification file
    protIds=[t.partition(sep)[0] for t in prots]
    return protIds

def filterPeptide(peptideObj, protIds, pepColumn):
    peptideObj['Is decoy']=peptideObj['Is decoy'].astype(str)
    #print(peptideObj['Is decoy'])
    protIdsEsc=[re.escape(prot+"_") for prot in protIds]
    protStr="|".join(protIdsEsc)
    ##removing decoy hits
    peptideObj=peptideObj[~peptideObj['Is decoy'].str.contains("TRUE",case=False)]
    
    ##removing PSMs only mapping to Contaminents
    peptideObj=peptideObj[~peptideObj[pepColumn].str.contains("^(CONT.*;?)+$")]
    filteredPeptideObj = peptideObj[peptideObj[pepColumn].str.contains(protStr)]
    return filteredPeptideObj

def writeResult(filteredObj, filename):
    filteredObj.to_csv(filename, index=False)
    
def main(args):
    protColName="description"
    pepColName="proteinacc_start_stop_pre_post_;"
    protObj1=readIdentifiedProteinPeptide(args.proteins,",")
    pepObj=readIdentifiedProteinPeptide(args.peptides,",")
    protObj=pepThresholding(protObj1, pepTh, protRevStr, protContStr)
    writeResult(protObj, args.proteinOutFile)
    ### Filter ORFs fasta file
    identifiedProt=identifiedProteinList(protObj, protColName)
    filterFasta(args.ORFs, identifiedProt, args.ORFsOutFile)
    
    ##Filter transcripts fasta file
    ##we are getting the ids again because the list has changed in the filterFasta function
    identifiedProt=identifiedProteinList(protObj, protColName)
    #print(identifiedProt)
    identifiedTrans=extractProtIdsFromDescription(identifiedProt, "|")
    unqTrans=list(set(identifiedTrans))
    
    filterFasta(args.transcripts, unqTrans, args.transcriptsOutFile)
    
    ###Filter PeptideObj
    ##we are getting the ids again because the list has changed in the filterFasta function
    identifiedProt=identifiedProteinList(protObj, protColName)
    protForPep=extractProtIdsFromDescription(identifiedProt, " ")
    fPeptide=filterPeptide(pepObj, protForPep, pepColName)
    writeResult(fPeptide, args.peptideOutFile)
    

parser = argparse.ArgumentParser(description='Fasta file filtering based on a header list given')
parser.add_argument("--ORFs", required=True, help="ORFs fasta file name")
parser.add_argument("--transcripts", required=True, help="Transcript fasta file name")
parser.add_argument("--proteins", required=True, help="mzIdentML-lib protein export file")
parser.add_argument("--peptides", required=True, help="mzIdentML-lib PSM export file")
parser.add_argument("--ORFsOutFile", required=True, help="ORFs output fasta file name")
parser.add_argument("--transcriptsOutFile", required=True, help="Transcripts output fasta file name")
parser.add_argument("--proteinOutFile", required=True, help="Protein csv out file name")
parser.add_argument("--peptideOutFile", required=True, help="Peptide csv out file name")

args = parser.parse_args()
protRevStr="XXX_"
protContStr="CONT_"
pepTh=1
main(args)

    
## Test Command
## python PIT-DBData_processing.py --ORFs D:\data\Oliver\ORFs\G10.assemblies.fasta.transdecoder.pep --transcripts D:\data\Oliver\transcripts\G10.assemblies.fasta --proteins D:\data\Oliver\PASA\G10\G10.assemblies.fasta.transdecoder.pep+fdr+th+grouping+prt.csvDBIds.csv --peptides D:\data\Oliver\PASA\G10\G10.assemblies.fasta.transdecoder.pep+fdr+th+grouping.csv
##--ORFsOutFile D:\data\Oliver\ORFs\G10.assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsOutFile D:\data\Oliver\transcripts\G10.assemblies.fasta.identified.fasta --peptideOutFile D:\data\Oliver\PASA\G10\G10.assemblies.fasta.transdecoder.pep+fdr+th+grouping_filtered.csv
