## this code merge the known, known var, iso and iso with var csv files produced by IdentifyProteinIsoformSAP.py for identified ORFs


import os
import argparse
import pandas as pd
import glob
import re

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;

def countVariation(vcf, orfID):
    var=vcf.loc[vcf['QueryID'] == orfID]
    return len(var)

parser = argparse.ArgumentParser(description='This python code merge ORF annotation files from IdentifyProteinIsoformSAP.py')
parser.add_argument("-b", "--blast", nargs=1, required=True, help="full path of blast folder", metavar="PATH")
parser.add_argument("-o", "--output", nargs=1, required=True, help="full path of the ouput folder", metavar="PATH")
parser.add_argument("-k", "--known", nargs=1, required=True, help="full path of to the known proteins csv file from the identified ORFs", metavar="PATH")
parser.add_argument("-s", "--knownvariation", nargs=1, required=True, help="full path of to the known with variation  csv file from identified ORFs", metavar="PATH")
parser.add_argument("-i", "--isoform", nargs=1, required=True, help="full path of to the isofrms csv file from the identified ORFs", metavar="PATH")
parser.add_argument("-j", "--isoformvariation", nargs=1, required=True, help="full path of to the isofrms with variations csv file from the identified ORFs", metavar="PATH")
parser.add_argument("-v", "--vcf", nargs=1, required=True, help="full path to the vcf file", metavar="PATH")

args = parser.parse_args()
#print(args)
#onlyfiles = [ f for f in os.listdir(args.blast[0]) if os.path.isfile(os.path.join(args.blast[0],f)) ]
#onlyfiles=glob.glob(args.blast[0]+"/*.assemblies.fasta.transdecoder.pep.identified.loc.csv")
#onlyfiles=glob.glob(args.blast[0]+"/*.assemblies.fasta.transdecoder.pep.csv")
#print(args.blast[0])
f=args.blast[0]
fBase=str(os.path.basename(f)).replace(".identified.loc.csv","")
print(fBase)
with open(args.output[0]+fBase+"_annotation.csv", 'w', newline='') as outfile:
	outfile.write("ORF Id,Protein ID,Class,Variation,Species\n");
	blast=readFile(f,',')
	known=readFile(args.known[0],',')
	kVar=readFile(args.knownvariation[0],',')
	iso=readFile(args.isoform[0],',')
	isoVar=readFile(args.isoformvariation[0],',')
	##identified ORF's vcf file should be read, not the peptide evidence one.
	vcf=readFile(args.vcf[0],'\t')
	info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','QueryLength','QueryStart','QueryEnd','SubjectLength','SubjectStart','SubjectEnd','Type','QPOS'])
	vcf=vcf.drop('INFO',1)
	vcf=vcf.join(info)
	vcf.SubjectID=vcf.SubjectID.str.replace('SubjectId=','')
	vcf.QueryID=vcf.QueryID.str.replace('QueryId=','')
	#vcf.Alignment=vcf.Alignment.str.replace('Alignment=\[','')
	#vcf.Alignment=vcf.Alignment.str.replace('\]','')
	vcf.QueryLength=vcf.QueryLength.str.replace('QueryLength=','')
	vcf.QueryStart=vcf.QueryStart.str.replace('QueryStart=','')
	vcf.QueryEnd=vcf.QueryEnd.str.replace('QueryEnd=','')
	vcf.SubjectLength=vcf.SubjectLength.str.replace('SubjectLength=','')
	vcf.SubjectStart=vcf.SubjectStart.str.replace('SubjectStart=','')
	vcf.SubjectEnd=vcf.SubjectEnd.str.replace('SubjectEnd=','')
	vcf.Type=vcf.Type.str.replace('Type:','')
	vcf.QPOS=vcf.QPOS.str.replace('QPOS:','')
	for i,row in known.iterrows():
		##If the protein is from Uniprot, following regex will work.
		m=re.search("(?<=OS=)(.*)(?= GN=)|(?<=OS=)(.*)(?= PE=)",row['Protein ID'])
		species=""
		if m:
			species=m.group(0)
		outfile.write(row['ORF Id']+","+row['Protein ID']+",known,0,"+species+"\n")
		        
	for i, row in kVar.iterrows():
		m=re.search("(?<=OS=)(.*)(?= GN=)|(?<=OS=)(.*)(?= PE=)",row['Protein ID'])
		species=""
		if m:
			species=m.group(0)
		#orfId=row['ORF Id']
		'''
		if not vcf.QueryID.str.contains(" "):
			## this means the ORFs ids as it was in the fasta file. nothing is removed after the first space.
			queryIDs=pd.DataFrame(row['ORF Id'].str.split(' ').tolist(),columns=['QueryID','GeneID','ORF','GeneID2','QueryID2','Type','Length','Strand','Location'])
		'''
		varCount=countVariation(vcf,row['ORF Id'])
		outfile.write(row['ORF Id']+","+row['Protein ID']+",known variation,"+str(varCount)+","+species+"\n")
    
	for i,row in iso.iterrows():
		##If the protein is from Uniprot, following regex will work.
		m=re.search("(?<=OS=)(.*)(?= GN=)|(?<=OS=)(.*)(?= PE=)",row['Protein ID'])
		species=""
		if m:
			species=m.group(0)
		outfile.write(row['ORF Id']+","+row['Protein ID']+","+row['Type']+",0,"+species+"\n")

	for i, row in isoVar.iterrows():
		m=re.search("(?<=OS=)(.*)(?= GN=)|(?<=OS=)(.*)(?= PE=)",row['Protein ID'])
		species=""
		if m:
			species=m.group(0)
		#orfId=row['ORF Id']
		'''
		if not vcf.QueryID.str.contains(" "):
			## this means the ORFs ids as it was in the fasta file. nothing is removed after the first space.
			queryIDs=pd.DataFrame(row['ORF Id'].str.split(' ').tolist(),columns=['QueryID','GeneID','ORF','GeneID2','QueryID2','Type','Length','Strand','Location'])
		'''
		varCount=countVariation(vcf,row['ORF Id'])
		outfile.write(row['ORF Id']+","+row['Protein ID']+","+row['Type']+","+str(varCount)+","+species+"\n")
	mapped=known['ORF Id'].tolist()+kVar['ORF Id'].tolist()+iso['ORF Id'].tolist()+isoVar['ORF Id'].tolist()
	novel=blast.loc[-blast['query_name'].isin(mapped)]
	novelList=novel['query_name'].tolist()
	print(len(novel['query_name'].tolist()))
	#print(novel['query_name'].tolist())
	for i in range(0,len(novelList)):
		outfile.write(novelList[i]+",-,novel,0,-\n")
