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

def writeResult(filteredObj, filename):
    filteredObj.to_csv(filename, index=False)
    
def main(args):
    
    protColName="description"
    pepColName="proteinacc_start_stop_pre_post_;"
    protObj1=readIdentifiedProteinPeptide(args.protein,",")
    
    pepObj=readIdentifiedProteinPeptide(args.peptide,",")
    protObj=pepThresholding(protObj1, pepTh, protRevStr, protContStr)
    #temporarily commented
    writeResult(protObj, args.protOut)
    ### Filter ORFs fasta file
    identifiedProt=identifiedProteinList(protObj, protColName)
    
    ###Filter PeptideObj
    ##we are getting the ids again because the list has changed in the filterFasta function
    identifiedProt=identifiedProteinList(protObj, protColName)
    protForPep=extractProtIdsFromDescription(identifiedProt, " ")
    fPeptide=filterPeptide(pepObj, protForPep, pepColName)
    writeResult(fPeptide, args.pepOut)
    
    filteredPrt=readIdentifiedProteinPeptide(args.protOut,',')
    protIds=filteredPrt["description"].tolist()
    filterFasta(args.fasta, protIds, args.fastaOut)

protRevStr="XXX_"
protContStr="CONT_"
pepTh=1

parser = argparse.ArgumentParser(description='MSGF+ protein peptide identification filtering for standard proteome search/correcting protein peptide filtering from a previous run')

parser.add_argument("--protein", required=True, help="Identified Protein file name")
parser.add_argument("--peptide", required=True, help="Identified Peptide file name")
parser.add_argument("--protOut", required=True, help="Identified Protein file name")
parser.add_argument("--pepOut", required=True, help="Identified Peptide file name")
parser.add_argument("--fasta", required=True, help="Standard fasta file name")
parser.add_argument("--fastaOut", required=True, help="Standard identified sequence output fasta file name")

args = parser.parse_args()
main(args)

##Human adeno-virus
#proteinFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/PASA/human_adeno_mydb_pasa.standard+fdr+th+grouping+prt.csv"
#peptideFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/PASA/human_adeno_mydb_pasa.standard+fdr+th+grouping.csv"