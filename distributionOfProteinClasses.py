##This Code readsTGE list for a dataset with multiple samples and count how many of them are known, novel, isoform, TrEMBL
## and so on.
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse
import glob
import pandas as pd
import numpy as np
import re

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;

def readFastaFile(fastaFile):
    ##This function reads a fasta file using SeqIO, convert entries into record object and put in dataframe.
    fastaHandle = open(fastaFile, "rU")
    records = list(SeqIO.parse(fastaHandle, "fasta"))
    #print("Records:"+str(len(records)))
    fastaDF=pd.DataFrame(columns=("ORF Id","Sequence"))
    for r in records:
        fastaDF=fastaDF.append({'ORF Id':re.sub("\s+"," ",r.description),'Sequence':str(r.seq)}, ignore_index=True)
        #print("seqList:"+str(len(seqList)))
    return fastaDF

def vcfVarReader(filename):
    vcf=readFile(filename, '\t')
    info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','QueryLength','QueryStart','QueryEnd','SubjectLength','SubjectStart','SubjectEnd','Type','QPOS','PeptideCount','UniquePeptideCount','Peptides','Score'])
    vcf=vcf.drop('INFO',1)
    vcf=vcf.join(info)
    vcf.SubjectID=vcf.SubjectID.str.replace('SubjectId=','')
    vcf.QueryID=vcf.QueryID.str.replace('QueryId=','')
    vcf.QueryLength=vcf.QueryLength.str.replace('QueryLength=','')
    vcf.QueryStart=vcf.QueryStart.str.replace('QueryStart=','')
    vcf.QueryEnd=vcf.QueryEnd.str.replace('QueryEnd=','')
    vcf.SubjectLength=vcf.SubjectLength.str.replace('SubjectLength=','')
    vcf.SubjectStart=vcf.SubjectStart.str.replace('SubjectStart=','')
    vcf.SubjectEnd=vcf.SubjectEnd.str.replace('SubjectEnd=','')
    vcf.Type=vcf.Type.str.replace('Type=','')
    vcf.QPOS=vcf.QPOS.str.replace('QPOS=','')
    vcf.PeptideCount=vcf.PeptideCount.str.replace('PeptideCount=','')
    vcf.UniquePeptideCount=vcf.UniquePeptideCount.str.replace('UniquePeptideCount=','')
    vcf.Peptides=vcf.Peptides.str.replace('Peptides=','')
    vcf.Score=vcf.Score.str.replace('Score=','')
    return vcf

def vcfIsoReader(filename):
    vcf=readFile(filename, '\t')
    info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','QueryLength','QueryStart','QueryEnd','SubjectLength','SubjectStart','SubjectEnd','Type','QPOS','PeptideCount','UniquePeptideCount','Peptides','Evidence','Score'])
    vcf=vcf.drop('INFO',1)
    vcf=vcf.join(info)
    vcf.SubjectID=vcf.SubjectID.str.replace('SubjectId=','')
    vcf.QueryID=vcf.QueryID.str.replace('QueryId=','')
    vcf.QueryLength=vcf.QueryLength.str.replace('QueryLength=','')
    vcf.QueryStart=vcf.QueryStart.str.replace('QueryStart=','')
    vcf.QueryEnd=vcf.QueryEnd.str.replace('QueryEnd=','')
    vcf.SubjectLength=vcf.SubjectLength.str.replace('SubjectLength=','')
    vcf.SubjectStart=vcf.SubjectStart.str.replace('SubjectStart=','')
    vcf.SubjectEnd=vcf.SubjectEnd.str.replace('SubjectEnd=','')
    vcf.Type=vcf.Type.str.replace('Type=','')
    vcf.QPOS=vcf.QPOS.str.replace('QPOS=','')
    vcf.PeptideCount=vcf.PeptideCount.str.replace('PeptideCount=','')
    vcf.UniquePeptideCount=vcf.UniquePeptideCount.str.replace('UniquePeptideCount=','')
    vcf.Peptides=vcf.Peptides.str.replace('Peptides=','')
    vcf.Evidence=vcf.Evidence.str.replace('Evidence=','')
    vcf.Score=vcf.Score.str.replace('Score=','')
    return vcf

parser = argparse.ArgumentParser(description='This python code counts unique number of proteins and peptides in a dataset')
parser.add_argument("-a", "--annot", nargs=1, required=True, help="full path of protein annotation folder", metavar="PATH")
parser.add_argument("-f", "--fasta", nargs=1, required=True, help="full path of fasta files folder", metavar="PATH")
parser.add_argument("-v", "--vcf", nargs=1, required=True, help="full path of polymorphism vcf folder", metavar="PATH")
parser.add_argument("-i", "--vcfiso", nargs=1, required=True, help="full path of isoform vcf folder", metavar="PATH")

args = parser.parse_args()
#print(args)
#onlyfiles = [ f for f in os.listdir(args.blast[0]) if os.path.isfile(os.path.join(args.blast[0],f)) ]
onlyfiles=glob.glob(args.annot[0]+"/*_details_annotation.csv")
print("All files")
print(onlyfiles)
#allProtein=readFile(args.fasta[0]+"TGEs.tsv", "\t")
#vcfDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PITDB/Variations-proVCF/"
proteinMat=pd.DataFrame(columns=('ORF Id', 'Protein ID', 'Class', 'Variation', 'Species','Protein Name', 'Gene Name', 'Protein description', 'Source', 'Sequence'))
for f in onlyfiles:
    #print("F:"+f)
    fBase=re.sub(re.escape("_details_annotation.csv"),"",f)
    #print("fBase"+fBase)
    sample=os.path.basename(fBase)#.split(".",maxsplit=1)[0]
    #sample="human_adeno"
    print("sample:"+str(sample))
    protAnnot=readFile(f, ',')
    protAnnot['Sample']=sample
    #vcf=readFile(args.vcf[0]+sample+".assemblies.fasta.transdecoder.pep_pepEvd.vcf", '\t')
    vcf=vcfVarReader(args.vcf[0]+sample+"_pepEvd.vcf")
    vcf['RefCount']=vcf['REF'].str.len()
    vcf['ALTCount']=vcf['ALT'].str.len()
    altVCF=vcf[(vcf['RefCount']>9) | (vcf['ALTCount']>9)]
    uniqIsoIds=altVCF['QueryID'].unique().tolist()
    protAnnot['mRNA']=protAnnot['ORF Id'].str.extract("([^ ]+)",expand=False)
    protAnnot['NewClass']=protAnnot['Class']
    protAnnot.loc[(protAnnot['mRNA'].isin(uniqIsoIds)) & (protAnnot['Class']=='known variation'),'NewClass']="ALT_SPLICE"
    
    
    ##Adding total of variations with peptide evidence in a TGE 
    pepVar=vcf[vcf['PeptideCount']!='0']['QueryID'].value_counts().to_frame()
    pepVar.columns=['VarPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    ##Adding total of variations with unique peptide evidence in a TGE 
    pepVar=vcf[vcf['UniquePeptideCount']!='0']['QueryID'].value_counts().to_frame()
    pepVar.columns=['UnqVarPeptide']
    pepVar['mRNA']=pepVar.index.values
    print("Prot Annot shape 1:")
    print(protAnnot.shape)
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    print("Prot Annot shape 2:")
    print(protAnnot.shape)
    ##Adding total of SAP variations with peptide evidence in a TGE
    pepVar=vcf[(vcf['Type']=='SAP') & (vcf['PeptideCount']!='0')]['QueryID'].value_counts().to_frame()
    pepVar.columns=['SAPPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    
    ##Adding total of SAP variations with unique peptide evidence in a TGE
    pepVar=vcf[(vcf['Type']=='SAP') & (vcf['UniquePeptideCount']!='0')]['QueryID'].value_counts().to_frame()
    pepVar.columns=['UnqSAPPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    
    ##Adding total of SSAP variations with peptide evidence in a TGE
    pepVar=vcf[(vcf['Type']=='SSAP') & (vcf['PeptideCount']!='0')]['QueryID'].value_counts().to_frame()
    pepVar.columns=['SSAPPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    
    ##Adding total of SSAP variations with unique peptide evidence in a TGE
    pepVar=vcf[(vcf['Type']=='SSAP') & (vcf['UniquePeptideCount']!='0')]['QueryID'].value_counts().to_frame()
    pepVar.columns=['UnqSSAPPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    
    ##Adding total of ALT variations with peptide evidence in a TGE
    pepVar=vcf[(vcf['Type']=='ALT') & (vcf['PeptideCount']!='0')]['QueryID'].value_counts().to_frame()
    pepVar.columns=['ALTPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    
    ##Adding total of SALT variations with unique peptide evidence in a TGE
    pepVar=vcf[(vcf['Type']=='ALT') & (vcf['UniquePeptideCount']!='0')]['QueryID'].value_counts().to_frame()
    pepVar.columns=['UnqALTPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    
    ##Adding total of SALT variations with peptide evidence in a TGE
    pepVar=vcf[(vcf['Type']=='SALT') & (vcf['PeptideCount']!='0')]['QueryID'].value_counts().to_frame()
    pepVar.columns=['SALTPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    
    ##Adding total of SALT variations with unique peptide evidence in a TGE
    pepVar=vcf[(vcf['Type']=='SALT') & (vcf['UniquePeptideCount']!='0')]['QueryID'].value_counts().to_frame()
    pepVar.columns=['UnqSALTPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    
    ##Adding total of INS variations with peptide evidence in a TGE
    pepVar=vcf[(vcf['Type']=='INS') & (vcf['PeptideCount']!='0')]['QueryID'].value_counts().to_frame()
    pepVar.columns=['INSPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    
    ##Adding total of INS variations with unique peptide evidence in a TGE
    pepVar=vcf[(vcf['Type']=='INS') & (vcf['UniquePeptideCount']!='0')]['QueryID'].value_counts().to_frame()
    pepVar.columns=['UnqINSPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    
    ##Adding total of DEL variations with peptide evidence in a TGE
    pepVar=vcf[(vcf['Type']=='DEL') & (vcf['PeptideCount']!='0')]['QueryID'].value_counts().to_frame()
    pepVar.columns=['DELPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    
    ##Adding total of DEL variations with unique peptide evidence in a TGE
    pepVar=vcf[(vcf['Type']=='DEL') & (vcf['UniquePeptideCount']!='0')]['QueryID'].value_counts().to_frame()
    pepVar.columns=['UnqDELPeptide']
    pepVar['mRNA']=pepVar.index.values
    protAnnot=pd.merge(protAnnot,pepVar,on='mRNA',how ='outer')
    
    isoVcf=vcfIsoReader(args.vcfiso[0]+sample+"_isoform_pepEvd.vcf")
    isoVcf['mRNA']=isoVcf['QueryID'].str.extract("([^ ]+)",expand=False)
    isoPepEvidence=isoVcf[isoVcf['PeptideCount']!='0']['mRNA'].unique()
    ##Creating a new column that will tell us which of the isoforms have isoform specific peptide evidence
    protAnnot['IsoPeptide']="No"
    protAnnot.loc[protAnnot['mRNA'].isin(isoPepEvidence),'IsoPeptide']="Yes"
    ##Creating a new column that will tell us which of the isoforms have isoform specific unique peptide evidence, i.e. the peptide did not match to any other TGE
    ##Though indication of this might have been wrong in the vcf file, when this TGE came from multiple transcript, it has been considered as multiple TGE, hence
    ##though the peptide essentially mapped to same TGE sequence it will come as non unique map.
    isoUnqPepEvidence=isoVcf[isoVcf['UniquePeptideCount']!='0']['mRNA'].unique()
    protAnnot['IsoUnqPeptide']="No"
    if len(isoUnqPepEvidence):
        protAnnot.loc[protAnnot['mRNA'].isin(isoUnqPepEvidence),'IsoUnqPeptide']="Yes"
    protDF=readFastaFile(args.fasta[0]+sample+".identified.fasta")
    print("protDF:")
    print(protDF.shape)
    #protDF=readFastaFile("/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PITDB/AminoAcids-or-ORFs-orTGEs/"+sample+".assemblies.fasta.transdecoder.pep.identified.fasta")
    print("protMerged 1:")
    if 'protMerged' in locals():
        print(protMerged.shape)
    protMerged=pd.merge(protAnnot,protDF, on='ORF Id', how='outer')
    print("protMerged 2:")
    print(protMerged.shape)
    print("proteinMat 1:")
    print(proteinMat.shape)
    proteinMat=proteinMat.append(protMerged)
    print("proteinMat 2:")
    print(proteinMat.shape)
'''
vcfFiles=glob.glob(vcfDir+"/*.assemblies.fasta.transdecoder.pep_pepEvd.vcf")
for f in onlyfiles:
    fBase=re.sub(re.escape(".assemblies.fasta.transdecoder.pep_details_annotation.csv"),"",f)
    #print("fBase"+fBase)
    sample=os.path.basename(fBase).split(".",maxsplit=1)[0]
    print("sample:"+str(sample))
    vcf=readFile(f, ',')
'''
proteinMat=proteinMat.fillna(0)
#unqPeptideEvdMa=proteinMat.drop_duplicates(['Protein ID', 'Class', 'Variation', 'Species','Protein Name', 'Gene Name', 'Protein description', 'Source', "NewClass",'Sequence'], keep='first')
#unqPeptideEvdMa.to_csv(args.annot[0]+"/uniqueProteins_pepEvd.csv")
#proteinMat.to_csv(args.annot[0]+"/AllProteins_pepEvd.csv")
#pepEvdUnique=proteinMat.drop_duplicates(['Protein ID', 'Class', 'Variation', 'Species','Protein Name', 'Gene Name', 'Protein description', 'Source', "NewClass",'Sequence','VarPeptide','UnqVarPeptide'], keep='first')
#pepEvdUnique.to_csv(args.annot[0]+"/VarPeptideProteins_pepEvd.csv")

colNames=list(proteinMat.columns.values)

#colNames.remove('ORF Id')
#colNames.remove('mRNA')
##no need to take unique, because I am reporting unique sequence count at the end and
##groupingby for peptide evidence counts.
unqPeptideEvdMat=proteinMat.drop_duplicates(colNames)
print("Sample 0")
print(unqPeptideEvdMat[unqPeptideEvdMat['Sample']==0].shape[0])
unqPeptideEvdMat.to_csv(args.annot[0]+"/TGEs_pepEvd.tsv", sep="\t", index=False)
aggregations={
    'ALTPeptide':'sum',
    'UnqALTPeptide':'sum',
    'SALTPeptide':'sum',
    'UnqSALTPeptide':'sum',
    'SAPPeptide':'sum',
    'UnqSAPPeptide':'sum',
    'SSAPPeptide':'sum',
    'UnqSSAPPeptide':'sum',
    'INSPeptide':'sum',
    'UnqINSPeptide':'sum',
    'DELPeptide':'sum',
    'UnqDELPeptide':'sum',
    'VarPeptide':'sum',
    'UnqVarPeptide':'sum',
    'IsoPeptide': lambda x: ",".join(x),
    'IsoUnqPeptide': lambda x: ",".join(x),
    'NewClass': lambda x: ",".join(x)
}
swissprot=unqPeptideEvdMat[unqPeptideEvdMat['Source']=='sp']
trEMBL=unqPeptideEvdMat[unqPeptideEvdMat['Source']=='tr']
pepEvdSwiss=swissprot.groupby('Sequence').agg(aggregations)
pepEvdTrEMBL=trEMBL.groupby('Sequence').agg(aggregations)

pepEvdSwiss.to_csv(args.annot[0]+"/swiss_pepEvd_grpby.tsv", sep="\t")
pepEvdTrEMBL.to_csv(args.annot[0]+"/trembl_pepEvd_grpby.tsv", sep="\t")
print("Total Protein:")
print(len(unqPeptideEvdMat['Sequence'].unique()))
print("Reviewed Protein:")
print(len(unqPeptideEvdMat.loc[(unqPeptideEvdMat['Class']=='known') & (unqPeptideEvdMat['Source']=='sp')]['Sequence'].unique()))
print("TrEMBL Protein:")
print(len(unqPeptideEvdMat.loc[(unqPeptideEvdMat['Class']=='known') & (unqPeptideEvdMat['Source']=='tr')]['Sequence'].unique()))
print("Canonical Protein:")
print(len(unqPeptideEvdMat.loc[(unqPeptideEvdMat['Class']=='known') & (unqPeptideEvdMat['Source']=='sp') & ~(unqPeptideEvdMat['Protein ID'].str.contains('-'))]['Sequence'].unique()))
print("Known Isoform Protein:")
print(len(unqPeptideEvdMat.loc[(unqPeptideEvdMat['Class']=='known') & (unqPeptideEvdMat['Source']=='sp') & (unqPeptideEvdMat['Protein ID'].str.contains('-'))]['Sequence'].unique()))


#print(unqPeptideEvdMat[['Class','NewClass']])
print("Novel Isoform Protein:")
print(len(unqPeptideEvdMat.loc[(unqPeptideEvdMat['NewClass']!='known') & (unqPeptideEvdMat['NewClass']!='novel') & (unqPeptideEvdMat['NewClass']!="known variation")]['Sequence'].unique()))

print("Novel Isoform Protein (Swissprot):")
print(len(unqPeptideEvdMat.loc[(unqPeptideEvdMat['NewClass']!='known') & (unqPeptideEvdMat['NewClass']!='novel') & (unqPeptideEvdMat['NewClass']!="known variation") & (unqPeptideEvdMat['Source']=='sp')]['Sequence'].unique()))

print("Novel Isoform Protein (TrEMBL):")
print(len(unqPeptideEvdMat.loc[(unqPeptideEvdMat['NewClass']!='known') & (unqPeptideEvdMat['NewClass']!='novel') & (unqPeptideEvdMat['NewClass']!="known variation") & (unqPeptideEvdMat['Source']=='tr')]['Sequence'].unique()))

print("Proteins with variations:")
print(len(unqPeptideEvdMat.loc[unqPeptideEvdMat['NewClass']=="known variation"]['Sequence'].unique()))
print("Proteins with variations(Swissprot):"+str(len(swissprot[swissprot['NewClass']=='known variation']['Sequence'].unique())))
print("Proteins with variations(TrEMBL):"+str(len(trEMBL[trEMBL['NewClass']=='known variation']['Sequence'].unique())))
print("Novel TGEs:")
print(len(unqPeptideEvdMat.loc[unqPeptideEvdMat['NewClass']=="novel"]['Sequence'].unique()))

##for the prptide evidence DF I dont need to take unique because this df is groupedby i.e. for unique sequence
#for following counts we have to separate between variations from 'known variation' and isoform classes.
print("Total number of TGEs with different types of variation/isoform specific peptide peptide evidence")
if pepEvdSwiss.shape[0]>0:
    print("Swissprot novel Isoform: "+str(pepEvdSwiss[pepEvdSwiss['IsoPeptide'].str.contains('Yes')].shape[0]))

if pepEvdTrEMBL.shape[0]>0:
    print("TrEMBL novel Isoform: "+str(pepEvdTrEMBL[pepEvdTrEMBL['IsoPeptide'].str.contains('Yes')].shape[0]))

if pepEvdSwiss.shape[0]>0:
    print("Swissprot Known protein with Variations"+str(pepEvdSwiss[(pepEvdSwiss['NewClass'].str.contains('known variation')) & (pepEvdSwiss['VarPeptide']>0)].shape[0]))

if pepEvdTrEMBL.shape[0]>0:
    print("TrEMBL Known protein with variations"+str(pepEvdTrEMBL[(pepEvdTrEMBL['NewClass'].str.contains('known variation')) & (pepEvdTrEMBL['VarPeptide']>0)].shape[0]))


print("Total number of TGEs with different types of variation/isoform specific unique peptide evidence")
if pepEvdSwiss.shape[0]>0:
    print("Swissprot novel Isoform"+str(pepEvdSwiss[(pepEvdSwiss['IsoUnqPeptide'].str.contains('Yes')) | ((pepEvdSwiss['NewClass'].str.contains('ALT_SPLICE')) & (pepEvdSwiss['UnqVarPeptide']>0))].shape[0]))

if pepEvdTrEMBL.shape[0]>0:
    print("TrEMBL novel Isoform"+str(pepEvdTrEMBL[(pepEvdTrEMBL['IsoUnqPeptide'].str.contains('Yes')) | ((pepEvdTrEMBL['NewClass'].str.contains('ALT_SPLICE')) & (pepEvdTrEMBL['UnqVarPeptide']>0))].shape[0]))

if pepEvdSwiss.shape[0]>0:
    print("Swissprot Known protein with Variations: "+str(pepEvdSwiss[(pepEvdSwiss['NewClass'].str.contains('known variation')) & (pepEvdSwiss['UnqVarPeptide']>0)].shape[0]))
if pepEvdTrEMBL.shape[0]>0:
    print("TrEMBL Known protein with variations: "+str(pepEvdTrEMBL[(pepEvdTrEMBL['NewClass'].str.contains('known variation')) & (pepEvdTrEMBL['UnqVarPeptide']>0)].shape[0]))

##Print isoform class/variation class wise peptide evidence
#print("Total variation ")
novelTges=unqPeptideEvdMat.loc[unqPeptideEvdMat['NewClass']=="novel"]

if len(onlyfiles)==1:
    prt=readFile("/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/PSMs-Peptides-ORFs/human_adeno+fdr+th+grouping+prt_filtered.csv","\t")
