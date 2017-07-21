### This code reads the identified protein/ORF list and filters SAMPLE.transdecoder.genome.gff3 and SAMPLE.transdecoder.gff3. This script
### goes on finding annotation of ORFs based on alt_splice events, sequence similary and location.

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import pandas as pd
import argparse
import re
from AminoAcidVariation import AminoAcidVariation

def readFile(filename, sep, headerFlag):
    if headerFlag==0:
        fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    elif headerFlag==1:
        fileDFObj = pd.read_table(filename, sep=sep, header=None, keep_default_na=False, na_values=[''])
    else:
        print("Unrecognized Header Flag")
    return fileDFObj;

def writeFile(dfObj, filename,sep,headerFlag):
    if headerFlag==0:
        dfObj.to_csv(path_or_buf=filename, sep=sep, na_rep='', float_format=None, header=True, index=False, quoting=None, quotechar='"')
    elif headerFlag==1:
        dfObj.to_csv(path_or_buf=filename, sep=sep, na_rep='', float_format=None, header=False, index=False, quoting=None, quotechar='"')
    else:
        print("Unrecognized Header Flag")

def filterPeptide(peptideObj, protIds, pepColumn):
    peptideObj['Is decoy']=peptideObj['Is decoy'].astype(str)
    #print(peptideObj['Is decoy'])
    protStr="|".join(protIds)
    ##removing decoy hits
    peptideObj=peptideObj[~peptideObj['Is decoy'].str.contains("TRUE",case=False)]
    ##removing PSMs only mapping to Contaminents
    peptideObj=peptideObj[~peptideObj[pepColumn].str.contains("^(CONT.*;?)+$")]
    ##removing PSMs only mapping to Decoy
    peptideObj=peptideObj[~peptideObj[pepColumn].str.contains("^(XXX.*;?)+$")]
    filteredPeptideObj = peptideObj[peptideObj[pepColumn].str.contains(protStr)]
    return filteredPeptideObj

def pepThresholding(prots, pepTh, protRevStr, protContStr):
    print("prots column names:"+str(prots.columns.values))
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

def filterFasta(fastaFile, overlapping_ids):
    fastahandle = open(fastaFile, "rU")
    overlapping_ids[:] = [line.strip() for line in overlapping_ids]
    overlapping_ids[:] = [line.replace('"','') for line in overlapping_ids]
    #print("oids:"+str(len(overlapping_ids)))
    #print(overlapping_ids[0])
    records = list(SeqIO.parse(fastahandle, "fasta"))
    #print("Total=")
    #print(len(records))
    count1=0
    count2=0
    filteredRecords=list()
    for record in records: # SeqIO.parse(handle, "fasta")
        #print(record.description)
        #print("TESTTEST")
        if record.description in overlapping_ids:
            overlapping_ids.remove(record.description)
            count1=count1+1
            filteredRecords.append(record)
    return filteredRecords

def extractProtIdsFromDescription(prots, sep):
    ##prots is a list of strings, generally the description column of protein identification file
    protIds=[t.partition(sep)[0] for t in prots]
    return protIds

def filterAltSpliceLabel(labelFile, transcriptIds):
    labels=readFile(labelFile,'\t',1)
    identifiedCluster=labels[labels[2].isin(set(labels[labels[1].isin(transcriptIds)][2].tolist()))]
    return identifiedCluster

def filterTransdecoderGenome(genomeGFFFile, geneIds):
    annotation=readFile(genomeGFFFile,'\t',1)
    annotation.columns=['seqid','source','type','start','end','score','strand','phase','attributes']
    geneIDs=[re.escape("ID="+g+";") for g in geneIds]
    ids="|".join(geneIDs)
    identifiedGenes=annotation[annotation.attributes.str.contains(ids)].index.tolist()
    identifiedAnnotation=list()
    for i in identifiedGenes:
        print("i="+str(i))
        identifiedAnnotation.append(i)
        if len(annotation)>=i+1:
            for j in range(i+1,len(annotation)):
                #print("J loop")
                if annotation['type'][j]=='gene':
                    break
                else:
                    #print("j="+str(j))
                    identifiedAnnotation.append(j)
    return annotation.loc[identifiedAnnotation]


#prtIds has to be a list
#transcriptIds=extractProtIdsFromDescription(prtIds, sep)



protRevStr="XXX_"
protContStr="CONT_"
pepTh=1
pepColumn='proteinacc_start_stop_pre_post_;'
parser = argparse.ArgumentParser(description='This code reads the identified protein/ORF list and filters SAMPLE.transdecoder.genome.gff3 and SAMPLE.transdecoder.gff3. This script goes on finding annotation of ORFs based on alt_splice events, sequence similary and location.')
parser.add_argument("-p", "--protein", nargs=1, required=True, help="full path of protein csv file", metavar="PATH")
parser.add_argument("-g", "--gff3", nargs=1, required=True, help="full path of the transdecoder genome gff3 file", metavar="PATH")
#parser.add_argument("-f", "--fasta", nargs=1, required=True, help="full path to the prf sequence fasta file", metavar="PATH")
#parser.add_argument("-j", "--isovar", nargs=1, required=True, help="full path to the isoforms with variation file", metavar="PATH")
#parser.add_argument("-p", "--pep", nargs=1, required=True, help="full path to the peptide identification file", metavar="PATH")
parser.add_argument("-o", "--output", nargs=1, required=True, help="full path to the output GFF3 file", metavar="PATH")
args = parser.parse_args()
prots=readFile(args.protein[0], ',', 0)
protsFiltered=pepThresholding(prots, pepTh, protRevStr, protContStr)
protsFiltered.description=protsFiltered.description.str.replace('\s+',' ')
description=pd.DataFrame(protsFiltered.description.str.split(' ').tolist(),columns=['QueryID','GeneID','ORF','GeneID2','QueryID2','Type','Length','Strand','Location'])
filteredGFF=filterTransdecoderGenome(args.gff3[0], description['GeneID'].tolist())

writeFile(filteredGFF, args.output[0],'\t',0)
