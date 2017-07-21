##This script calculates peptide coverage of 'known' proteins.
import argparse
import re
import pandas as pd

def readFile(filename,sep, header):
    #reads delimited file and create a panda object
    if header==1:
    	fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    else:
    	fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''], header=None)
    return fileDFObj
def peptideCoverage(knownPeptide, knownProteins):
	#find peptides that has K or R except at the last position
	miscleavedPeptide=knownPeptide[knownPeptide['Sequence'].str.contains('[K|R]+.+?[K|R]$')]
	splitMiscleaved=miscleavedPeptide['Sequence'].str.extractall("(.+?[K|R|$])")
	splitMiscleavedFiltered=splitMiscleaved[splitMiscleaved[0].str.len()>=8][0].unique()
	##check how many of these miscleaved peptide contains tryptic peptide, miscleavage=0, that have been identified.
	pepForCoverage=knownPeptide[~knownPeptide['Sequence'].isin(splitMiscleavedFiltered)]
	pepForCoverage['PepLength']=pepForCoverage['Sequence'].str.len()
	pepForCoverageUnique=pepForCoverage.drop_duplicates(['Sequence'])
	knownProteinLength=pd.Series(knownProteins['ORF Id'].str.extract("len:(\d+)", expand=False))
	proteinAAsCount=pd.to_numeric(knownProteinLength).sum()
	peptideAAsCount=pepForCoverageUnique['PepLength'].sum()
	return peptideAAsCount/proteinAAsCount
def main(pepFile, annotFile):
	##reads annotation and MSGF+ peptide csv file. This function call peptideCoverage function to calculate 
	##peptide coverage of known proteins.
	peptides=readFile(pepFile,',',1)
	annot=readFile(annotFile,',',1)
	knownProteins=annot[annot['Class']=='known']
	knownProteinsStr="|".join([re.escape(x) for x in knownProteins['ORF Id'].str.extract("([^ ]+)", expand=False).tolist()])
	knownPeptide=peptides[peptides['proteinacc_start_stop_pre_post_;'].str.contains(knownProteinsStr)]

pepFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/PSMs-Peptides-ORFs/human_adeno+fdr+th+grouping_filtered.csv"
annotFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Annotation/human_adeno_mydb_pasa.assemblies.fasta.transdecoder.pep_details_annotation.csv"