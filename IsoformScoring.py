###Isoform scoring based on sample specific consequence score

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os
import re
import argparse
import pandas as pd
#from itertools import chain

def readFile(filename,sep):
    #reads delimited file and create a panda object 
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj

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
    vcf['mRNA']=vcf['QueryID'].str.extract("([^ ]+)", expand=False)
    return vcf

def readFastaFile(fastaFile):
    ##This function reads a fasta file using SeqIO, convert entries into record object and put in dataframe.
    fastaHandle = open(fastaFile, "r")
    records = list(SeqIO.parse(fastaHandle, "fasta"))
    print("Records:"+str(len(records)))
    fastaDF=pd.DataFrame(columns=("Id","Sequence"))
    for r in records:
        fastaDF=fastaDF.append({'Id':re.sub("\s+"," ",r.description),'Sequence':str(r.seq)}, ignore_index=True)
        #print("seqList:"+str(len(seqList)))
    return fastaDF

def addMissingPeptideToDigest(peptides, digestedDF, flag):
    #check All Identified Peptide Exists In Digested List
    #peptides is a dataframe containing MSGF+ peptide file
    #digestedDF is a dataframe with four columns, 'Protein','Peptide','Start', 'Stop'
    unqPep=pd.Series(peptides['Sequence'].unique())
    if len(unqPep[~unqPep.isin(digestedDF['Peptide'])])>0:
        ##not all of the identified peptides are there.
        missingPep=unqPep[~unqPep.isin(digestedDF['Peptide'])].tolist()
        print("Total Number of missing peptides:"+str(len(missingPep)))
        
        missingProts=peptides.loc[peptides['Sequence'].isin(missingPep)][['Sequence','proteinacc_start_stop_pre_post_;']].drop_duplicates()
        s = missingProts['proteinacc_start_stop_pre_post_;'].str.split(';').apply(pd.Series, 1).stack()
        s.index = s.index.droplevel(-1)
        s.name = 'proteinacc_start_stop_pre_post_;'
        del missingProts['proteinacc_start_stop_pre_post_;']
        missingProts=missingProts.join(s)
        missingProts=missingProts.drop_duplicates()
        missingProts=missingProts[~missingProts['proteinacc_start_stop_pre_post_;'].str.contains("^CONT|^XXX_",na=False)]
        protacc_start_stop_pre_post=missingProts['proteinacc_start_stop_pre_post_;'].str.extract("(.*)_(\d+)_(\d+)_([A-Z]|-)_([A-Z]|-)",expand=False)
        protacc_start_stop_pre_post.columns=["Protein","Start","Stop","Pre","Post"]
        ##remove 'sp|' or 'tr|' and '|PROTNAME' if this was a standard serach/uniprot search
        if flag==1:
            protacc_start_stop_pre_post['Protein']=protacc_start_stop_pre_post['Protein'].str.extract(".*\|([^\|]+)\|.*",expand=False)
        del missingProts['proteinacc_start_stop_pre_post_;']
        missingProts=missingProts.join(protacc_start_stop_pre_post)
        ##renaming Sequence column to peptide because thats how the digest columns are named.
        ##Also taking only those columns that I need
        missingProts=missingProts.rename(columns={'Sequence':'Peptide'})[['Peptide','Protein','Start','Stop']]
        missingProts['Identified']='Identified'
        merged=[digestedDF,missingProts[["Protein","Peptide","Start","Stop","Identified"]]]
        newDigestedDF=pd.concat(merged)
        return newDigestedDF
    else:
        print("all found")
        return digestedDF

def peptideWise(peptides):
    ##This function converts PSM wise dataframe into a peptide wise data frame by keeping one PSM per peptide where the selected PSM has lowest q-value.
    peptideWiseObj=peptides.groupby(['Sequence']).agg({'PSM-level q-value':min})
    peptideWiseObj['Peptide']=peptideWiseObj.index
    return peptideWiseObj

def trainTestModel(annot, orfs, refs, digestFile, isoDigestFile, conseqFile, isoConseqFile, peptideFile, rpeptideFile, sssFile):
    ##this function digest ORFs classified as known protein based on their blast map
    ##and calculates consequence score, sample specific score. Using these two values
    ##the R code trains a model and saves it. Isoform specific sequences are then
    ##digested in the same way, consequence values are computed and SampleSpecific
    ##scores are predicted using the mkodel.
    #change following filter for known canonical TGEs
    #knownOrfs=filterFastaDF(orfs, annot[(annot['Class']=="known") & (annot['Source']=="sp") & ~(annot['Protein ID'].str.contains("-"))]["ORF Id"])
    ##For Bat we include trEMBL proteins
    #knownOrfIds=annot[(annot['Class']=="known") & (annot['Source']=="sp") & ~(annot['Protein ID'].str.contains("-"))]["mRNA"]
    knownOrfAnnot=annot[(annot['Class']=="known") & ~(annot['Protein ID'].str.contains("-"))]
    print("annot:"+str(annot.shape[0]))
    print("knownOrfAnnot:"+str(knownOrfAnnot.shape[0]))
    print("ORF ID:"+knownOrfAnnot.iloc[0]['ORF Id'])
    knownOrfs=filterFastaDF(orfs, knownOrfAnnot['ORF Id'])
    print("knowOrfs:"+str(knownOrfs.shape[0]))
    #sTime=time.clock()
    knownDigested=knownOrfs.apply(digest,1)
    print("known Digested shape:"+str(knownDigested.shape[0]))
    digestedDF=pd.DataFrame()
    for i in range(knownDigested.shape[0]):
        #print(i)
        #print(knownDigested.iloc[i])
        digestedDF=digestedDF.append(pd.DataFrame(knownDigested.iloc[i]),ignore_index=True)
    #eTime=time.clock()
    #print(digestedDF.columns.values)
    digestedDF=digestedDF[['Protein','Peptide','Start','Stop']]
    
    ##Check if all the identified peptides from known orfs are there in the digested peptide dataframe
    #knownOrfIds=annot[annot['Class']=="known"]["mRNA"]
    
    knownOrfIdsStr="|".join([re.escape(mrna)+"_" for mrna in  knownOrfAnnot['mRNA']])
    peptide=readFile(peptideFile,',')
    rpeptide=readFile(rpeptideFile,',')
    knownPeptide=peptide[peptide['proteinacc_start_stop_pre_post_;'].str.contains(knownOrfIdsStr)]
    
    digestedDF=addMissingPeptideToDigest(knownPeptide, digestedDF,0)
    ##Write the digested ORFs to a file so that the consequence tool can read it.
    digestedDF.to_csv(digestFile,sep="\t", index=False)
    
    ##Run consequence.r
    conDir="/data/home/btw796/Code2/Proteomics/Consequence"
    os.chdir(conDir)
    conCommand="Rscript code/consequence.r "+digestFile+" "+conseqFile+" "+conDir+" rf"
    os.system(conCommand)
    ##Now put together peptides from rest of the TGEs.
    trainORFIds=annot[~annot['ORF Id'].isin(knownOrfAnnot['ORF Id'])]
    #trainORFIds=annot[annot['Class']=="known"]['ORF Id']
    restOrfs=filterFastaDF(orfs, trainORFIds["ORF Id"])
    #sTime=time.clock()
    restDigested=restOrfs.apply(digest,1)
    restDigestedDF=pd.DataFrame()
    for i in range(restDigested.shape[0]):
        restDigestedDF=restDigestedDF.append(pd.DataFrame(restDigested.iloc[i]),ignore_index=True)
    restDigestedDF['Identified']='Unidentified'
    restDigestedDF.loc[restDigestedDF['Peptide'].isin(peptide['Sequence']),'Identified']='Identified'
    restDigestedDF=restDigestedDF[['Protein','Peptide','Start','Stop','Identified']]
    ##Check if all the identified peptides from the rest of the orfs are there in the digested peptide dataframe
    #restOrfIds=annot[~(annot['ORF Id'].isin(trainORFIds)) & (annot['Class']!="novel")]["mRNA"]
    restOrfIds=trainORFIds["mRNA"]
    restOrfIdsStr="|".join([re.escape(mrna)+"_" for mrna in  restOrfIds])
    restPeptide=peptide[peptide['proteinacc_start_stop_pre_post_;'].str.contains(restOrfIdsStr)]
    restDigestedDF=addMissingPeptideToDigest(restPeptide, restDigestedDF,0)
    peptideWiseVar=peptideWise(peptide)
    restDigestedDF=restDigestedDF.merge(peptideWiseVar, how='left', on='Peptide')
    restDigestedDF.loc[restDigestedDF['PSM-level q-value'].isnull(),'PSM-level q-value']=1
    #Now digest reference proteins that has a map to the ORFs
    rIds=annot[~(annot['ORF Id'].isin(trainORFIds['ORF Id'])) & (annot['Class']!="novel")]["Protein ID"].unique()
    restOrfsRef=filterFastaDF(refs, rIds)
    restRefDigested=restOrfsRef.apply(digest,1)
    restRefDigestedDF=pd.DataFrame()
    for i in range(restRefDigested.shape[0]):
        restRefDigestedDF=restRefDigestedDF.append(pd.DataFrame(restRefDigested.iloc[i]),ignore_index=True)
    restRefDigestedDF['Identified']='Unidentified'
    restRefDigestedDF.loc[restRefDigestedDF['Peptide'].isin(rpeptide['Sequence']),'Identified']='Identified'
    restRefDigestedDF=restRefDigestedDF[['Protein','Peptide','Start','Stop','Identified']]
    ##Check if all the identified peptides from known orfs are there in the digested peptide dataframe
    #restOrfRefIds=annot[~(annot['ORF Id'].isin(trainORFIds)) & (annot['Class']!="novel")]["mRNA"]
    restOrfRefIdsStr="|".join([re.escape(mrna)+"_" for mrna in  rIds])
    ##not all of these would be found in standard search. Because there is PIT only identification. This 
    ##point has to be considered while doing the scoring.
    restRefPeptide=rpeptide[rpeptide['proteinacc_start_stop_pre_post_;'].str.contains(restOrfRefIdsStr)]
    restRefDigestedDF=addMissingPeptideToDigest(restRefPeptide, restRefDigestedDF,1)
    refPeptideWiseVar=peptideWise(rpeptide)
    restRefDigestedDF=restRefDigestedDF.merge(refPeptideWiseVar, how='left', on='Peptide')
    restRefDigestedDF.loc[restRefDigestedDF['PSM-level q-value'].isnull(),'PSM-level q-value']=1
    print("restRefDigestedDF shape")
    print(restRefDigestedDF.shape)
    #print(restRefDigestedDF.columns.values)
    #Merge rest ORFs peptides with reference peptides
    restDigestedDF=restDigestedDF.append(restRefDigestedDF)
    
    ##Write the digested ORFs to a file so that the consequence tool can read it.
    restDigestedDF.to_csv(isoDigestFile,sep="\t", index=False)
    
    ##Run consequence.r
    conDir="/data/home/btw796/Code2/Proteomics/Consequence"
    os.chdir(conDir)
    conCommand="Rscript code/consequence.r "+isoDigestFile+" "+isoConseqFile+" "+conDir+" rf"
    os.system(conCommand)
    ##Pass this matrix to train polinomial
    rCode="/data/home/btw796/Code2/Proteomics/R"
    os.chdir(rCode)
    sssCommand="Rscript SampleSpecificScoringandModelling.R "+conseqFile+" "+peptideFile+" "+isoConseqFile+" "+sssFile
    os.system(sssCommand)


def scoring(vcfFile, vcfIsoFile, annot, sssFile, restOut):
    ##This function computes score of an novel isoform or known protein with variation.
    ##Extract peptides of the part of reference or ORFs sequence that does not map to each other
    #trainORFIds=annot[(annot['Class']=="known") & (annot['Source']=="sp") & ~(annot['Protein ID'].str.contains("-"))]['ORF Id']

    ##We need to apply different scoring formula for non-overlappig part from TGE and reference sequence
    ##extra/alternative part from TGE is the part we consider to be 'peptides that should exist' and
    ##reference seq only part at 'peptides that should not exist'.
    #annot[(annot['Class']=="known") & (annot['Source']=="sp") & ~(annot['Protein ID'].str.contains("-"))]["ORF Id"]
    testORFIds=annot[(annot['Class']!="known") & (annot['Class']!="novel")]['ORF Id']
    restOrfsAnnotOrg=annot[(annot['ORF Id'].isin(testORFIds))]
    restOrfsAnnotOrg.loc[:,'Score']=0
    restOrfsAnnotOrg.loc[:,'RefScore']=0
    restOrfsAnnot=restOrfsAnnotOrg #restOrfsAnnotOrg.reset_index()

    testSSS=readFile(sssFile,'\t')
    #print(testSSS.loc[1:5])
    vcf=vcfVarReader(vcfFile)
    vcfIso=vcfIsoReader(vcfIsoFile)
    vcfIso['SubjectID']=vcfIso['SubjectID'].str.extract(".*\|([^\|]+)\|.*",expand=False)
    vcfIso['QueryID']=vcfIso['QueryID'].str.extract("([^ ]+)",expand=False)
    for i in restOrfsAnnot.index.values: #"asmbl_10881|m.135267""asmbl_6449|m.83120","asmbl_15021|m.190447","asmbl_4103|m.42327"].index.values:#restOrfsAnnot.index.values:
        ##For each TGE, find all variations, SAPs, ALT, INDELs, boundary variations.
        #orfId=restOrfsAnnot.iloc[i]['mrna']
        ##If this one has class "known variation" then check vcf file for variation location.
        #print("i:"+str(i))
        if restOrfsAnnot.loc[i]['Class']=="known variation":
            #We only need to check vcf
            #if this orf is not found in the vcf file, that means the variation happend
            #at the boundary of the ORF, but length of the variation was max of 2 AAs.
            orfVars=vcf[vcf['QueryID']==restOrfsAnnot.loc[i]['mRNA']]
            if orfVars.shape[0]>0:
                #internal variations, we can calculate the score here.
                varPeps=getPeptideMutation(orfVars,testSSS)
                ##varPeps is a dictionary with two keys, 'ref' and 'orf'
                ##calculate score for this orf and corresponding ref seq.
                #scores=isoformScore2(varPeps)
                #scores=isoformScore3(varPeps)
                scores=isoformScore4(varPeps)
                restOrfsAnnot.loc[i,'Score']=scores['orf']
                restOrfsAnnot.loc[i,'RefScore']=scores['ref']
                if scores['orf']==0 and scores['ref']==0:
                    print("For further Investigation (known variation):"+str(restOrfsAnnot.loc[i]['mRNA']))
            else:
                print("the variation happend at the boundary of the ORF, but length of the variation was max of 2 AAs.")
        elif restOrfsAnnot.loc[i]['Variation']==0 or restOrfsAnnot.loc[i]['Variation']=='0':
            ##These are novel-isoforms without any internal variations. Hence I should only check vcfIso.
            orfVars=vcfIso[vcfIso['QueryID']==restOrfsAnnot.loc[i]['mRNA']]
            isoVarPeps=getPeptideBoundary(orfVars,testSSS)
            ##isoVarPeps is a dictionary with two keys, 'ref' and 'orf'
            #scores=isoformScore2(isoVarPeps)
            #scores=isoformScore3(varPeps)
            scores=isoformScore4(isoVarPeps)
            restOrfsAnnot.loc[i,'Score']=scores['orf']
            restOrfsAnnot.loc[i,'RefScore']=scores['ref']
            if scores['orf']==0 and scores['ref']==0:
                print("For further Investigation (isoform):"+str(restOrfsAnnot.loc[i]['mRNA']))
        else:
            ##these are novvel-isoforms with variations.
            orfVars=vcf[vcf['QueryID']==restOrfsAnnot.loc[i]['mRNA']]
            varPeps=getPeptideMutation(orfVars,testSSS)
            print("varPeps:")
            #print(varPeps['ref'].shape)
            #print(varPeps)
            #scores1=isoformScore2(varPeps)
            #scores1=isoformScore3(varPeps)
            scores1=isoformScore4(varPeps)
            ##varPeps is a dictionary with two keys, 'ref' and 'orf'
            orfVars=vcfIso[vcfIso['QueryID']==restOrfsAnnot.loc[i]['mRNA']]
            isoVarPeps=getPeptideBoundary(orfVars,testSSS)
            print("isVarPeps:")
            #print(isoVarPeps.shape)
            #print(isoVarPeps)
            if restOrfsAnnot.loc[i]['mRNA']=='asmbl_4103|m.42327':
                print("Polymorpic peptides:")
                #print(varPeps)
                print("Iso peptides:")
                #print(isoVarPeps)
            ##isoVarPeps is a dictionary with two keys, 'ref' and 'orf'
            #scores2=isoformScore2(isoVarPeps)
            #scores2=isoformScore3(isoVarPeps)
            scores2=isoformScore4(isoVarPeps)
            restOrfsAnnot.loc[i,'Score']=scores1['orf']+scores2['orf']
            restOrfsAnnot.loc[i,'RefScore']=scores1['ref']+scores2['ref']
            if scores1['orf']==0 and scores1['ref']==0:
                print("For further Investigation (isoform variation):"+str(restOrfsAnnot.loc[i]['mRNA']))
            if scores2['orf']==0 and scores2['ref']==0:
                print("For further Investigation (Isoform variation):"+str(restOrfsAnnot.loc[i]['mRNA']))
    restOrfsAnnot.to_csv(restOut,sep='\t')


def refPeptide(SSS, prtId, refVarPeps, start, end):
    print("REFPEPTIDE function")
    print("refVarPeps object")
    #print(refVarPeps)
    if refVarPeps[refVarPeps['Identified']=='Identified'].shape[0]==0:
        ##No peptide was identified for standard/reference protein either. Hence the score entirely depends on SSS score.
        #refVarPepGrp=refVarPeps.groupby('Start')
        ##take non overlapping peptides to compute the score.
        print("NO identified peptide")
        seqRefVarPeps=nonOverlappingPeptide(refVarPeps, start, end)
        if isinstance(seqRefVarPeps, pd.DataFrame) and seqRefVarPeps.shape[1]==1:
            seqRefVarPeps=seqRefVarPeps.transpose()
        else:
            if isinstance(seqRefVarPeps, pd.Series):
                seqRefVarPepsDF=pd.DataFrame(seqRefVarPeps)
                if seqRefVarPepsDF.shape[1]==1:
                    seqRefVarPeps=seqRefVarPepsDF.transpose()
        print("seqRefVarPeps:"+str(seqRefVarPeps.shape[0]))
        return seqRefVarPeps
    else:
        ##For each identified peptide, I need to call nonoverlappingPeptide for start=end+1 of previous identified peptide and end=strt-1 of current 
        ##So we subset refVarPeps
        #print("Print in refPeptide, before sort")
        #seqRefVarPeps=pd.DataFrame()
        print("Identified peptide")
        refIdent=refVarPeps[refVarPeps['Identified']=='Identified']
        refIdentSorted=refIdent.sort_values(by='Start')
        #print(refIdentSorted)
        #print("Print in refPeptide, after sort")
        seqRefVarPeps=pd.DataFrame()
        for j in range(0,refIdentSorted.shape[0]):
            if int(refIdentSorted.iloc[j]['Start'])!=start:
                varPeps=extractPeptide(refVarPeps, prtId, start, int(refIdentSorted.iloc[j]['Start'])-1)
                if varPeps.shape[0]>0:
                    nonOverPep=nonOverlappingPeptide(varPeps, start, int(refIdentSorted.iloc[j]['Start'])-1)
                    if isinstance(nonOverPep,pd.DataFrame) and nonOverPep.shape[1]==1:
                        nonOverPep=nonOverPep.transpose()
                    else:
                        if isinstance(nonOverPep,pd.Series):
                            nonOverPepDF=pd.DataFrame(nonOverPep)
                            if nonOverPepDF.shape[1]==1:
                                nonOverPep=nonOverPepDF.transpose()
                    if seqRefVarPeps.shape[0]>0:
                        seqRefVarPeps=seqRefVarPeps.append(nonOverPep)
                    else:
                        seqRefVarPeps=nonOverPep #nonOverlappingPeptide(varPeps, start, int(refIdentSorted.iloc[j]['Start'])-1)
                    print("J:"+str(j)+ "seqRefVarPeps:")
                    #print(seqRefVarPeps)
            #seqRefVarPeps=seqRefVarPeps.append(refIdentSorted.iloc[j])
            start=int(refIdentSorted.iloc[j]['Stop'])+1
            print("Adding Identified peptide")
            identPep=refIdentSorted.iloc[j]
            if isinstance(identPep,pd.Series):
                identPep=pd.DataFrame(identPep).transpose()
            else:
                if isinstance(identPep,pd.Series) and identPep.shape[1]==1:
                    identPep=identPep.transpose()
            if seqRefVarPeps.shape[0]>0:
                if isinstance(seqRefVarPeps,pd.Series):
                    seqRefVarPeps=pd.DataFrame(seqRefVarPeps).transpose()
                else:
                    if isinstance(seqRefVarPeps,pd.DataFrame) and seqRefVarPeps.shape[1]==1:
                        seqRefVarPeps=seqRefVarPeps.transpose()
                seqRefVarPeps=seqRefVarPeps.append(identPep)
            else:
                seqRefVarPeps=identPep
            print("After adding identified peptide seqRefVarPeps")
            print(seqRefVarPeps)
            #print(seqRefVarPeps)
        ##following code is to add all the peptides from end of last identified peptide till end of refVarPeps list.
        #end=max(refVarPeps['Start'])
        if start<=end:
            varPeps=extractPeptide(refVarPeps, prtId, start, end)
            if varPeps.shape[0]>0:
                nonOverPep=nonOverlappingPeptide(varPeps, start, end)
                if isinstance(nonOverPep,pd.DataFrame) and nonOverPep.shape[1]==1:
                        nonOverPep=nonOverPep.transpose()
                else:
                    if isinstance(nonOverPep,pd.Series):
                        nonOverPepDF=pd.DataFrame(nonOverPep)
                        if nonOverPepDF.shape[1]==1:
                            nonOverPep=nonOverPepDF.transpose()
                if seqRefVarPeps.shape[0]>0:
                    seqRefVarPeps=seqRefVarPeps.append(nonOverPep)
                else:
                    seqRefVarPeps=nonOverPep #nonOverlappingPeptide(varPeps, start, end)
        else:
            print("Last identified covers till the end of the variation")
        print("In refPeptide seqRefVarPeps")
        #print(seqRefVarPeps)
        return seqRefVarPeps

def getPeptideBoundary(orfVars,SSS):
    ##orfVars is a dataframe of vcf file, containing variations of ORFS. One ORF can have multiple
    ##variations.
    ##For all variation boundary extract all peptides and remove unnecessary miscleaved peptides.
    print("GETPEPTIDEBOUNDARY")
    refVarNonOverPeps=pd.DataFrame()
    orfVarNonOverPeps=pd.DataFrame()    
    ##SSS=readFile("/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/AminoAcids-or-O                          RFs-orTGEs/human_adeno.assemblies.fasta.transdecoder.pep.identified.fasta_conseq_sss.tsv",sep='\t')
    for i in range(0,orfVars.shape[0]):
        type=orfVars.iloc[i]['Type']
        print("Type:"+type)
        prtId=orfVars.iloc[i]['SubjectID']
        ##extracting peptide from the reference protein from the location of the variation to check if that has any identified peptide.
        if '5prime' in type:
            if int(orfVars.iloc[i]['SubjectStart'])==1:
                ##do not have ref only portion. I need to the scoring part here. Otherwise I lose this information.
                rstart=-1
                rstop=-1
            else:
                rstart=1
                rstop=int(orfVars.iloc[i]['SubjectStart'])-1
            if int(orfVars.iloc[i]['QueryStart'])==1:
                ##do not have ref only portion. I need to the scoring part here. Otherwise I lose this information.
                ostart=-1
                ostop=-1
            else:
                ostart=1
                ostop=int(orfVars.iloc[i]['QueryStart'])-1
        elif '3prime' in type:
            if int(orfVars.iloc[i]['SubjectEnd'])==int(orfVars.iloc[i]['SubjectLength']):
                ##do not have ref only portion. I need to the scoring part here. Otherwise I lose this information.
                rstart=-1
                rstop=-1
            else:
                rstart=int(orfVars.iloc[i]['SubjectEnd'])+1
                rstop=int(orfVars.iloc[i]['SubjectLength'])
            if int(orfVars.iloc[i]['QueryEnd'])==int(orfVars.iloc[i]['QueryLength']):
                ##do not have ref only portion. I need to the scoring part here. Otherwise I lose this information.
                ostart=-1
                ostop=-1
            else:
                ostart=int(orfVars.iloc[i]['QueryEnd'])+1
                ostop=int(orfVars.iloc[i]['QueryLength'])
        else:
            print("ERROR: Iso variation type should contain either 5prime or 3prime")
        refVarPeps=extractPeptide(SSS, prtId, rstart, rstop)
        print("getPeptideBooundary:refVarpeps before refPeptide")
        #print(refVarPeps)
        if refVarPeps.shape[0]>0:
            if refVarNonOverPeps.shape[0]>0:
                print("Calling refPeptide with refVarNonOverPeps, non-zero refVarNonOverPeps")
                refVarNonOverPeps=refVarNonOverPeps.append(refPeptide(SSS, prtId, refVarPeps, rstart, rstop))
                #print(refVarNonOverPeps)
            else:
                print("Calling refPeptide with refVarNonOverPeps, zero refVarNonOverPeps")
                refVarNonOverPeps=refPeptide(SSS, prtId, refVarPeps, rstart, rstop)
                print("refVarNonOverPeps shape 0")
                #print(refVarNonOverPeps)
        ##TGE specific peptides, even though they have not been identified, we need to know these for scoring
        orfId=orfVars.iloc[i]['QueryID']
        orfVarPeps=extractPeptide(SSS, orfId, ostart, ostop)
        print("getPeptideBoundary: ostart n ostop:"+str(ostart)+","+str(ostop))
        print("orfVarPeps:")
        #print(orfVarPeps)
        if orfVarPeps.shape[0]>0:
            if orfVarNonOverPeps.shape[0]>0:
                print("Calling refPeptide with orfVarNonOverPeps, non-zero orfVarNonOverPeps")
                orfVarNonOverPeps=orfVarNonOverPeps.append(refPeptide(SSS, orfId, orfVarPeps,ostart, ostop))
                print("orfVarNonOverPeps is zero size")
                #print(orfVarNonOverPeps)
            else:
                print("Calling refPeptide with orfVarNonOverPeps, zero orfVarNonOverPeps")
                orfVarNonOverPeps=refPeptide(SSS, orfId, orfVarPeps,ostart, ostop)
                print("orfVarNonOverPeps shape 0")
                #print(orfVarNonOverPeps)
    return {'ref':refVarNonOverPeps, 'orf':orfVarNonOverPeps}

def getPeptideMutation(orfVars,SSS):
    ##orfVars is a dataframe of vcf file, containing variations of ORFS. One ORF can have multiple
    ##variations.
    ##For all variation boundary extract all peptides and remove unnecessary miscleaved peptides.
    print("GETPEPTIDEMUTATION")
    refVarNonOverPeps=pd.DataFrame()
    orfVarNonOverPeps=pd.DataFrame()
    ##SSS=readFile("/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/AminoAcids-or-ORFs-orTGEs/human_adeno.assemblies.fasta.transdecoder.pep.identified.fasta_conseq_sss.tsv",sep='\t')
    ##for consequitive variations, we will have repetations of peptides. I should take unique for the Dataframe
    for i in range(0,orfVars.shape[0]):
        ##check peptideCount value to quickly figure out if we have peptide evidence.
        ##Check ref peptides for the 'POS' location, ref seq gets extra points if 'POS' has peptide evidence.
        prtId=orfVars.iloc[i]['SubjectID']
		##extracting peptide from the reference protein from the location of the variation to check if that has any identified peptide.
		#start and stop location used here is suitable for internal variations, I need to change it for boundary variations. because 'ALT' and 'REF' is going to have '.'.
        refVarPeps=extractPeptide(SSS, prtId, int(orfVars.iloc[i]['POS']), len(orfVars.iloc[i]['REF'])+int(orfVars.iloc[i]['POS'])-1)
        print("refvarpeps:"+str(refVarPeps.shape[0]))
        #print(refVarPeps)
        if refVarPeps.shape[0]>0:
            ##this may happen due to too long peptide covering the area.
            if refVarNonOverPeps.shape[0]>0:
                refVarNonOverPeps=refVarNonOverPeps.append(refPeptide(SSS, prtId, refVarPeps, int(orfVars.iloc[i]['POS']), len(orfVars.iloc[i]['REF'])+int(orfVars.iloc[i]['POS'])-1))
            else:
                refVarNonOverPeps=refPeptide(SSS, prtId, refVarPeps, int(orfVars.iloc[i]['POS']), len(orfVars.iloc[i]['REF'])+int(orfVars.iloc[i]['POS'])-1)
        else:
            print("This reference protein does not have peptide covering the variation location:"+str(prtId))
        print("refVarNonOverPeps:"+str(refVarNonOverPeps.shape[0]))
        ##TGE specific peptides, even though they have not been identified, we need to know these for scoring
        orfId=orfVars.iloc[i]['QueryID']
        orfVarPeps=extractPeptide(SSS, orfId, int(orfVars.iloc[i]['QPOS']), len(orfVars.iloc[i]['ALT'])+int(orfVars.iloc[i]['QPOS'])-1)
        if orfVarPeps.shape[0]>0:
            ##This may happen due to the length restrictions on the peptide.
            if orfVarNonOverPeps.shape[0]>0:
                orfVarNonOverPeps=orfVarNonOverPeps.append(refPeptide(SSS, orfId, orfVarPeps, int(orfVars.iloc[i]['QPOS']), len(orfVars.iloc[i]['ALT'])+int(orfVars.iloc[i]['QPOS'])-1))
            else:
                orfVarNonOverPeps=refPeptide(SSS, orfId, orfVarPeps, int(orfVars.iloc[i]['QPOS']), len(orfVars.iloc[i]['ALT'])+int(orfVars.iloc[i]['QPOS'])-1)
        else:
            print("This ORF does not have peptide covering the variation location:"+str(orfId))
    if isinstance(refVarNonOverPeps,pd.DataFrame) and isinstance(orfVarNonOverPeps,pd.DataFrame):
        return {'ref':refVarNonOverPeps, 'orf':orfVarNonOverPeps}
    else:
        if isinstance(refVarNonOverPeps,pd.Series):
            refVarNonOverPeps=pd.DataFrame(refVarNonOverPeps).transpose()
        if isinstance(orfVarNonOverPeps,pd.Series):
            orfVarNonOverPeps=pd.DataFrame(orfVarNonOverPeps).transpose()
        return {'ref':refVarNonOverPeps, 'orf':orfVarNonOverPeps}

def nonOverlappingPeptide(nonOverlapPeps, start, end):
	##this function will take a dataframe of peptides and start and end and return nonoverlappiing/sequencial peptides sothat its easier to calculate the score.
	##Here start and end is start and end of the identified peptide.
	#nonOverlapPeps=peps[~ ((peps['Start']<=start) & (peps['Stop']>=start))|((peps['Start']<=end) & (peps['Stop']>=end))]
	#if there is an identified peptide, I need to start finding sequential from start and end of the identified peptide. 
	#nonOverlapPeps=nonOverlapPeps.sort(columns="Start", inplace=True)
	#we start from the peptide that has smallest start. next peptide is which has start=end of previous peptide+1
	#nonOverlapPepsGrp=nonOverlapPeps.groupby("Start")
	##find first and shortest peptide (smallest end), take end of that and find next peptide which has start=end of this peptide+1
    print("NONOVERLAPPINGPEPTIDE")
    #print(nonOverlapPeps)
    nonOverlapPeps['Start']=nonOverlapPeps['Start'].astype('int')
    nonOverlapPeps['Stop']=nonOverlapPeps['Stop'].astype('int')
    nonOverlapPeps['distFromStart']=start-nonOverlapPeps['Start']
    #print(nonOverlapPeps)
    if nonOverlapPeps[nonOverlapPeps['distFromStart']>=0].shape[0]>0:
        firstGrp=nonOverlapPeps[nonOverlapPeps['distFromStart']==min(nonOverlapPeps[nonOverlapPeps['distFromStart']>=0]['distFromStart'])]
    else:
        print("part of the variation does not have any detectable peptide")
        ##here we are taking the max because the distances are negetive, -1,-2 and so on. We want the closest one from the start.
        firstGrp=nonOverlapPeps[nonOverlapPeps['distFromStart']==max(nonOverlapPeps['distFromStart'])]
    currPep=firstGrp[firstGrp['Stop']==min(firstGrp['Stop'])]
    seqPeps=pd.DataFrame()
    end=int(nonOverlapPeps[nonOverlapPeps['Start']==max(nonOverlapPeps['Start'])]['Start'].unique())
    if currPep.shape[0]>0:
        seqPeps=currPep.iloc[0]
        if currPep.shape[0]>1:
            print("there are 2 peptides with same start and length, this should not happen")
            #print(currPep)

        while int(currPep.iloc[0]['Stop'])<=end:
            nextPeps=nonOverlapPeps[nonOverlapPeps['Start']==int(currPep.iloc[0]['Stop'])+1]
            if nextPeps.shape[0]==0:
                print("There is no peptides starting at stop+1 location:"+str(currPep.iloc[0]['Stop']))
                break
            currPep=nextPeps[nextPeps['Stop']==min(nextPeps['Stop'])]
            if currPep.shape[0]==1:
                if isinstance(seqPeps,pd.DataFrame) and seqPeps.shape[1]==1:
                    seqPeps=seqPeps.transpose()
                else:
                    if isinstance(seqPeps,pd.Series):
                        seqPeps=pd.DataFrame(seqPeps).transpose()
                    else:
                        print("can not determine type of seqPeps")
                seqPeps=seqPeps.append(currPep.iloc[0])
            else:
                print("This should not happen")
    else:
        print("currPep is empty")
    print("In nonOverlappingPeptide")
    #print(seqPeps)
    return seqPeps
    '''
    if isinstance(seqPeps, pd.DataFrame):
        return seqPeps
    else:
        print("seqPeps is not dataframe. convert it to DataFrame")
        if isinstance(seqPeps, pd.Series):
            print("seqPeps is Series")
            print(pd.DataFrame([seqPeps],columns=seqPeps.index.values).transpose().index.values)
            return pd.DataFrame([seqPeps]).transpose() #.to_frame()
        else:
            print("Check seqPeps further")
    '''

	
def digest(seqObj):
    #SeqObj is a row of fastaObj. This function digest the sequence and
    #returns a dataframe object with four columns, 'Protein','Peptide','Start','Stop'
    seq=seqObj['Sequence']
    #print("Seq:"+seq)
    seq=seq.rstrip('*')
    frags=list(filter(None,re.findall("(.*?(?:K|R|$))",seq)))
    misCleavage=10 ##This is what I have used while running Jun's digest
    #missedFrags=pd.DataFrame([], columns=["Protein","Peptide", "Start", "Stop"])#[[frags[i],"".join(frags[i:i+2]),"".join(frags[i:i+3]),"".join(frags[i:i+4]),"".join(frags[i:i+5]),"".join(frags[i:i+6]),"".join(frags[i:i+7]),"".join(frags[i:i+8]),"".join(frags[i:i+9]),"".join(frags[i:i+10]),"".join(frags[i:i+11])] for i in range(0,len(frags)-10)]
    peptides=[]
    #print(frags)
    m= re.search("^([^ ]+)",seqObj['Id'])
    protId=""
    if m:
        protId=m.group(1)
    else:
        print("ERROR")
    for i in range(0,len(frags)):
        count=0
        if i==0:
            start=1
        else:
            start=len("".join(frags[:i]))+1
        if (len(frags[i])>=8) & (len(frags[i])<=40):
            stop=start+len(frags[i])-1
            #missedFrags=missedFrags.append(pd.Series({"Protein":seqObj['Id'],"Peptide":frags[i],"Start":start,"Stop":stop}),ignore_index=True)##check the login behind 1
            if "X" in frags[i]:
                print("Following sequence has peptide containing X:"+seqObj['Id']+" : "+frags[i])
            else:
                peptides.append({"Protein":protId,"Peptide":frags[i],"Start":start,"Stop":stop})
            '''
        if i==len(frags)-1:
            print(frags[i])
            '''
        if i!=len(frags)-1:
            for j in range(i+1,len(frags)):
                if count<=misCleavage:
                    pep="".join(frags[i:j+1])
                    if (len(pep)>=8) & (len(pep)<=40):
                        stop=start+len(pep)-1
                        #missedFrags=missedFrags.append(pd.Series({"Protein":seqObj['Id'],"Peptide":pep,"Start":start,"Stop":stop}),ignore_index=True)
                        if "X" in pep:
                            print("Following sequence has peptide containing X:"+seqObj['Id']+" : "+pep)
                        else:
                            peptides.append({"Protein":protId,"Peptide":pep,"Start":start,"Stop":stop})
                        count=count+1
                    elif len(pep)>40:
                        break
                else:
                    break
    #missedFrags=pd.DataFrame(peptides)
    #return missedFrags
    ##remove peptides with 'X' character

    return peptides
#digest(seqObj.iloc[i])

def extractPeptide(digestedDF, orfId, start, stop):
    #finds all the peptides from digestedDF that overlaps with the start and stop location
    print("\n\nEXTRACTPEPTIDE")
    print("Start and Stop")
    print(start)
    print(stop)
    subDigstedDF = digestedDF[digestedDF['Protein']==orfId]
    if subDigstedDF.shape[0]==0:
        print("The reference protein was not identified:"+orfId)
        return subDigstedDF
    if stop-start==0:
        ##i.e. single point variation
        print("single point variation or truncation/extention")
        subDigstedDF = subDigstedDF[(subDigstedDF['Start']<=start) & (subDigstedDF['Stop']>=start)]
    else:
        ##added peptides that could cover the whole variation region for DEL, INS and ALT events.
        print("Variation")
        subDigstedDF = subDigstedDF[((subDigstedDF['Start']>=start) & (subDigstedDF['Start']<=stop)) | ((subDigstedDF['Stop']>=start) & (subDigstedDF['Stop']<=stop)) | ((subDigstedDF['Start']<=start) & (subDigstedDF['Stop']>=stop))]
    #print(subDigstedDF)
    return subDigstedDF

def isoformScore(pepDict):
    ##This function calculates score of an orf and its correspnding reference protein.
    ##pepDict is a dictionary with two keys, 'ref' and 'orf' containing all the variation/isoform
    ##specific peptides along with their identification status and sample specific score.
    orfPeps=pepDict['orf']
    refPeps=pepDict['ref']
    print("In isoformSCore")
    #print(orfPeps.columns.values)
    #print(len(orfPeps.columns.values))
    print("shape of orfPeps")
    print(orfPeps.shape)
    #print(orfPeps)
    if orfPeps.shape[0]>0:
        if refPeps.shape[0]>0:
            orfScore=(orfPeps[orfPeps['Identified']=='Identified'].shape[0]+((1-orfPeps[orfPeps['Identified']!='Identified']['sss']).sum()/2)+((0.5+refPeps[refPeps['Identified']!='Identified']['sss']).sum()/2))/(orfPeps.shape[0]+refPeps.shape[0])
            refScore=(refPeps[refPeps['Identified']=='Identified'].shape[0]+((1-refPeps[refPeps['Identified']!='Identified']['sss']).sum()/2)+((0.5+orfPeps[orfPeps['Identified']!='Identified']['sss']).sum()/2))/(orfPeps.shape[0]+refPeps.shape[0])
        else:
            print("RefPeps Empty")
            refScore=0
            orfScore=(orfPeps[orfPeps['Identified']=='Identified'].shape[0]+((1-orfPeps[orfPeps['Identified']!='Identified']['sss']).sum()/2))/(orfPeps.shape[0])
    else:
        if refPeps.shape[0]>0:
            print("OrfPeps Empty")
            orfScore=0
            refScore=(refPeps[refPeps['Identified']=='Identified'].shape[0]+((1-refPeps[refPeps['Identified']!='Identified']['sss']).sum()/2))/(refPeps.shape[0])
        else:
            print("RefPeps and OrfPeps both Empty, HOW????")
            refScore=0
            orfScore=0
    return {'orf':orfScore, 'ref':refScore}

def isoformScore2(pepDict):
    ##This function calculates score of an orf and its correspnding reference protein.
    ##pepDict is a dictionary with two keys, 'ref' and 'orf' containing all the variation/isoform
    ##specific peptides along with their identification status and sample specific score.
    orfPeps=pepDict['orf']
    refPeps=pepDict['ref']
    print("In isoformSCore")
    #print(orfPeps.columns.values)
    #print(len(orfPeps.columns.values))
    print("shape of orfPeps")
    print(orfPeps.shape)
    print(orfPeps.columns.values)
    if orfPeps.shape[0]>0:
        if refPeps.shape[0]>0:
            orfScore=(sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value'])-((0.5+orfPeps[orfPeps['Identified']!='Identified']['sss']).sum()/2)+((0.5+refPeps[refPeps['Identified']!='Identified']['sss']).sum()/2) - sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value']))/(orfPeps.shape[0]+refPeps.shape[0])
            refScore=(sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value'])-((0.5+refPeps[refPeps['Identified']!='Identified']['sss']).sum()/2)+((0.5+orfPeps[orfPeps['Identified']!='Identified']['sss']).sum()/2) - sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value']))/(orfPeps.shape[0]+refPeps.shape[0])
        else:
            print("RefPeps Empty")
            refScore=(((0.5+orfPeps[orfPeps['Identified']!='Identified']['sss']).sum()/2) - sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value']))/(orfPeps.shape[0])
            orfScore=(sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value'])-((0.5+orfPeps[orfPeps['Identified']!='Identified']['sss']).sum()/2))/(orfPeps.shape[0])
    else:
        if refPeps.shape[0]>0:
            print("OrfPeps Empty")
            orfScore=(((0.5+refPeps[refPeps['Identified']!='Identified']['sss']).sum()/2) - sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value']))/(refPeps.shape[0])
            refScore=(sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value'])-((0.5+refPeps[refPeps['Identified']!='Identified']['sss']).sum()/2))/(refPeps.shape[0])
        else:
            print("RefPeps and OrfPeps both Empty, HOW????")
            refScore=0
            orfScore=0
    return {'orf':orfScore, 'ref':refScore}

def isoformScore3(pepDict):
    ##This function calculates score of an orf and its correspnding reference protein.
    ##pepDict is a dictionary with two keys, 'ref' and 'orf' containing all the variation/isoform
    ##specific peptides along with their identification status and sample specific score.
    orfPeps=pepDict['orf']
    refPeps=pepDict['ref']
    print("In isoformSCore")
    #print(orfPeps.columns.values)
    #print(len(orfPeps.columns.values))
    print("shape of orfPeps")
    print(orfPeps.shape)
    print(orfPeps.columns.values)
    if orfPeps.shape[0]>0:
        if refPeps.shape[0]>0:
            orfScore=(sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value'])-(orfPeps[orfPeps['Identified']!='Identified']['sss'].sum()/2)+(refPeps[refPeps['Identified']!='Identified']['sss'].sum()/2) - sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value'])/2)/(orfPeps.shape[0]+refPeps.shape[0])
            refScore=(sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value'])-(refPeps[refPeps['Identified']!='Identified']['sss'].sum()/2)+(orfPeps[orfPeps['Identified']!='Identified']['sss'].sum()/2) - sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value'])/2)/(orfPeps.shape[0]+refPeps.shape[0])
        else:
            print("RefPeps Empty")
            refScore=((orfPeps[orfPeps['Identified']!='Identified']['sss'].sum()/2) - sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value'])/2)/(orfPeps.shape[0])
            orfScore=(sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value'])-(orfPeps[orfPeps['Identified']!='Identified']['sss'].sum()/2))/(orfPeps.shape[0])
    else:
        if refPeps.shape[0]>0:
            print("OrfPeps Empty")
            orfScore=((refPeps[refPeps['Identified']!='Identified']['sss'].sum()/2) - sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value'])/2)/(refPeps.shape[0])
            refScore=(sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value'])-(refPeps[refPeps['Identified']!='Identified']['sss'].sum()/2))/(refPeps.shape[0])
        else:
            print("RefPeps and OrfPeps both Empty, HOW????")
            refScore=0
            orfScore=0
    return {'orf':orfScore, 'ref':refScore}

def isoformScore4(pepDict):
    ##This function calculates score of an orf and its correspnding reference protein.
    ##pepDict is a dictionary with two keys, 'ref' and 'orf' containing all the variation/isoform
    ##specific peptides along with their identification status and sample specific score.
    orfPeps=pepDict['orf']
    refPeps=pepDict['ref']
    if isinstance(orfPeps,pd.Series) or (isinstance(orfPeps,pd.DataFrame) and orfPeps.shape[1]==1):
        print("in isoformScoring4, orfPeps is Series")
        orfPeps=pd.DataFrame(orfPeps).transpose()
    if orfPeps.shape[0]>1:
        orfPeps=orfPeps.drop_duplicates(['Protein','Peptide','Start','Stop'])
    if isinstance(orfPeps,pd.Series) or (isinstance(orfPeps,pd.DataFrame) and orfPeps.shape[1]==1):
        print("in isoformScoring4, orfPeps is Series")
        orfPeps=pd.DataFrame(orfPeps).transpose()
    if isinstance(refPeps,pd.Series) or (isinstance(refPeps,pd.DataFrame) and refPeps.shape[1]==1):
        print("in isoformScoring4, orfPeps is Series")
        refPeps=pd.DataFrame(refPeps).transpose()
    if refPeps.shape[0]>1:
        refPeps=refPeps.drop_duplicates(['Protein','Peptide','Start','Stop'])
    if isinstance(refPeps,pd.Series) or (isinstance(refPeps,pd.DataFrame) and refPeps.shape[1]==1):
        print("in isoformScoring4, orfPeps is Series")
        refPeps=pd.DataFrame(refPeps).transpose()
    print("In isoformSCore")
    #print(orfPeps.columns.values)
    #print(len(orfPeps.columns.values))
    print("orfPeps")
    #print(orfPeps.shape)
    #print(refPeps.shape)
    #print(orfPeps.columns.values)
    print(orfPeps)
    print("refPeps")
    print(refPeps)
    print("Duplicated index:\nORFPEPS")
    print(orfPeps[orfPeps.index.duplicated()])
    print("REFSPEPS")
    print(refPeps[refPeps.index.duplicated()])
    if orfPeps.shape[0]>0:
        if refPeps.shape[0]>0:
            print("orf q-val:"+str(sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value'])))
            print("orf unidentified:"+str((orfPeps[orfPeps['Identified']!='Identified']['sss'].sum()/8)))
            print("ref unidentified:"+str((refPeps[refPeps['Identified']!='Identified']['sss'].sum()/8)))
            print("ref q-val:"+str(sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value'])/4))
            orfScore=(sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value'])-(orfPeps[orfPeps['Identified']!='Identified']['sss'].sum()/8)+(refPeps[refPeps['Identified']!='Identified']['sss'].sum()/8) - sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value'])/4)/(orfPeps.shape[0]+refPeps.shape[0])
            refScore=(sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value'])-(refPeps[refPeps['Identified']!='Identified']['sss'].sum()/8)+(orfPeps[orfPeps['Identified']!='Identified']['sss'].sum()/8) - sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value'])/4)/(orfPeps.shape[0]+refPeps.shape[0])
        else:
            print("RefPeps Empty")
            refScore=((orfPeps[orfPeps['Identified']!='Identified']['sss'].sum()/8) - sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value'])/4)/(orfPeps.shape[0])
            orfScore=(sum(1-orfPeps[orfPeps['Identified']=='Identified']['PSM.level.q.value'])-(orfPeps[orfPeps['Identified']!='Identified']['sss'].sum()/8))/(orfPeps.shape[0])
    else:
        if refPeps.shape[0]>0:
            print("OrfPeps Empty")
            orfScore=((refPeps[refPeps['Identified']!='Identified']['sss'].sum()/8) - sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value'])/4)/(refPeps.shape[0])
            refScore=(sum(1-refPeps[refPeps['Identified']=='Identified']['PSM.level.q.value'])-(refPeps[refPeps['Identified']!='Identified']['sss'].sum()/8))/(refPeps.shape[0])
        else:
            print("RefPeps and OrfPeps both Empty, HOW????")
            refScore=0
            orfScore=0
    return {'orf':orfScore, 'ref':refScore}


def filterFastaDF(fastaDF, ids):
    ##this function filter fasta dataframe based on the ids passed
    print("Fasta ID:"+fastaDF.iloc[0]['Id'])
    return fastaDF[fastaDF['Id'].isin(ids)]

'''
def main(args):
    #vcf=vcfIsoReader(args.vcfIso)
    annot=readFile(args.annot,',')
    annot['mRNA']=annot['ORF Id'].str.extract("([^ ]+)", expand=False)
    
    orfs=readFastaFile(args.orfs)

    #proteins=readFile(args.prot,',')
    refs=readFastaFile(args.refs)
    refs['Id']=refs['Id'].str.extract(".*?\\|([^\\|]+)\\|.*?", expand=False)
    ##Trains a model with consequence and sample specific score of peptides from known proteins.
    trainTestModel(annot, orfs, refs, args.digesttrain, args.digesttest, args.conseqtrain, args.conseqtest, args.pep, args.rpep, args.ssstest)
    scoring(vcfFile, vcfIsoFile, annot, sssFile)
'''


'''
aFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Annotation/human_adeno_mydb_pasa.assemblies.fasta.transdecoder.pep_details_annotation.csv"

oFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/AminoAcids-or-ORFs-orTGEs/human_adeno.assemblies.fasta.transdecoder.pep.identified.fasta"
rFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/AminoAcids-or-ORFs-orTGEs/human_adeno_mydb_pasa.standard.identified.fasta"
digestFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.train.digest.tsv"
isoDigestFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.test.digest.tsv"
conseqFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.train.conseq.tsv"
#isoConseqFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.test.conseq.tsv"
isoConseqFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.test.conseq.scoreFunction2.tsv"
peptideFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/PSMs-Peptides-ORFs/human_adeno+fdr+th+grouping_filtered.csv"
rpeptideFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/PASA/human_adeno_mydb_pasa.standard+fdr+th+grouping_filtered.csv"
#sssFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.test.sss.tsv"
sssFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.test.sss.scoreFunction2.tsv"
vcfFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Variations-proVCF/human_adeno.assemblies.fasta.transdecoder.pep_pepEvd.vcf"
vcfIsoFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Annotation/human_adeno_mydb_pasa.assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf"
'''
##Canonical
'''
aFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/blast/human_adeno_mydb_pasa.assemblies.fasta.transdecoder.pep.canonical_details_annotation.csv"

oFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/AminoAcids-or-ORFs-orTGEs/human_adeno.assemblies.fasta.transdecoder.pep.identified.fasta"
rFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/AminoAcids-or-ORFs-orTGEs/human_adeno_mydb_pasa.standard.identified.fasta"
digestFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.canonical.train.digest.tsv"
isoDigestFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.canonical.test.digest.tsv"
conseqFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.canonical.train.conseq.tsv"
isoConseqFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.canonical.test.conseq.tsv"
peptideFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/PSMs-Peptides-ORFs/human_adeno+fdr+th+grouping_filtered.csv"
rpeptideFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/PASA/human_adeno_mydb_pasa.standard+fdr+th+grouping_filtered.csv"
'''
#sssFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.test.sss.tsv"
#sssFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.canonical.test.sss.tsv"
'''
sssFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/human_adeno.assemblies.fasta.transdecoder.pep.identified.test.sss.scoreFunction2.tsv"
vcfFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/blast/human_adeno_mydb_pasa.assemblies.fasta.transdecoder.pep.canonical_pepEvd.vcf"
vcfIsoFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/blast/human_adeno_mydb_pasa.assemblies.fasta.transdecoder.pep.canonical_details_annotation.csv_isoform_pep.vcf"
restOut="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/conseq/scoresvNewFunction4.tsv"
'''

###Bat
'''
aFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PITDB/Annotation/11_C5LH8ACXX_GGCTAC_L005.assemblies.fasta.transdecoder.pep_details_annotation.csv"
oFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PITDB/AminoAcids-or-ORFs-orTGEs/11_C5LH8ACXX_GGCTAC_L005.assemblies.fasta.transdecoder.pep.identified.fasta"
rFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/identification/batNelson/Palecto_nelsonBay.standard.identified.fasta"
digestFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/conseq/11_C5LH8ACXX_GGCTAC_L005.assemblies.fasta.transdecoder.pep.identified.train.digest.tsv"
isoDigestFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/conseq/11_C5LH8ACXX_GGCTAC_L005.fasta.transdecoder.pep.identified.test.digest.tsv"
conseqFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/conseq/11_C5LH8ACXX_GGCTAC_L005.fasta.transdecoder.pep.identified.train.conseq.tsv"
isoConseqFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/conseq/11_C5LH8ACXX_GGCTAC_L005.assemblies.fasta.transdecoder.pep.identified.test.conseq.scoreFunction4.tsv"
peptideFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PITDB/PSMs-Peptides-ORFs/11_C5LH8ACXX_GGCTAC_L005+fdr+th+grouping_filtered.csv"
rpeptideFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/identification/batNelson/batNelson.standard+fdr+th+grouping_filtered.csv"

sssFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/conseq/11_C5LH8ACXX_GGCTAC_L005.assemblies.fasta.transdecoder.pep.identified.test.sss.scoreFunction4.tsv"
vcfFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PITDB/Variations-proVCF/11_C5LH8ACXX_GGCTAC_L005.assemblies.fasta.transdecoder.pep_pepEvd.vcf"
vcfIsoFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PITDB/Annotation/11_C5LH8ACXX_GGCTAC_L005.assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf"
restOut="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/conseq/11_C5LH8ACXX_GGCTAC_L005.scoresvNewFunction4.tsv"


'''
parser = argparse.ArgumentParser(description='Isoform Scoring')
parser.add_argument("--annot", required=True, help="Protein annotation file name")
parser.add_argument("--vcf", required=True, help="Peptide csv out file name")
parser.add_argument("--vcfiso", required=True, help="Isoform peptide evidence vcf file name")
parser.add_argument("--pep", required=True, help="Peptide csv out file name")
####parser.add_argument("--prot", required=True, help="Protein csv out file name")
parser.add_argument("--rpep", required=True, help="Peptide csv out file name")
####parser.add_argument("--rprot", required=True, help="Protein csv out file name")
parser.add_argument("--digesttrain", required=True, help="Digest tsv out file name")
parser.add_argument("--digesttest", required=True, help="Digest score tsv out file name")
parser.add_argument("--conseqtrain", required=True, help="Consequence score tsv out file name")
parser.add_argument("--conseqtest", required=True, help="Consequence score tsv out file name")
parser.add_argument("--ssstest", required=True, help="Sample specific consequence score tsv out file name")
parser.add_argument("--orfs", required=True, help="Identified ORFs fasta file name")
parser.add_argument("--refs", required=True, help="Identified ORFs corresponding reference seq fasta file name")
parser.add_argument("--score", required=True, help="Score file name")



args = parser.parse_args()
aFile=args.annot
oFile=args.orfs
rFile=args.refs
'''
main(args)
#orfs=readFastaFile("/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/AminoAcids-or-ORFs-orTGEs/human_adeno.assemblies.fasta.transdecoder.pep.identified.fasta")
'''

annot=readFile(aFile,',')
annot['mRNA']=annot['ORF Id'].str.extract("([^ ]+)", expand=False)
#'''
orfs=readFastaFile(oFile)
refs=readFastaFile(rFile)
refs['Id']=refs['Id'].str.extract(".*?\\|([^\\|]+)\\|.*?", expand=False)
print("ORFS Size:"+str(orfs.shape[0]))
print("REFS Size:"+str(refs.shape[0]))
#'''
##Trains a model with consequence and sample specific score of peptides from known proteins.
#trainTestModel(annot, orfs, refs, args.digesttrain, args.digesttest, args.conseqtrain, args.conseqtest, args.pep, args.rpep, args.ssstest)
#trainTestModel(annot, orfs, refs, digestFile, isoDigestFile, conseqFile, isoConseqFile, peptideFile, rpeptideFile, sssFile)
#scoring(vcfFile, vcfIsoFile, annot, sssFile, restOut)
trainTestModel(annot, orfs, refs, args.digesttrain, args.digesttest, args.conseqtrain, args.conseqtest, args.pep, args.rpep, args.ssstest)
scoring(args.vcf, args.vcfiso, annot, args.ssstest, args.score)
