## this code finds peptide evidence of isoforms.

import pandas as pd
import argparse
import re
from AminoAcidVariation import AminoAcidVariation

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;

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
    print("1. prots dim:"+str(prots.shape))
    prots=prots[~prots['protein accession'].str.contains(protRevStr)]
    print("2. prots dim:"+str(prots.shape))
    prots=prots[~prots['protein accession'].str.contains(protContStr)]
    print("3. prots dim:"+str(prots.shape))
    prots=prots[prots['distinct peptide sequences']>pepTh]
    print("4. prots dim:"+str(prots.shape))
    return prots

def filterBlast(blastFile, queryList):
    blast=readFile(blastFile,",")
    blast['query_name'].isin(queryList)
    identIdx=blast['query_name'].isin(queryList)
    idx=identIdx[identIdx==True].index.tolist()
    blastIdentified = blast.iloc[idx,:]
    return blastIdentified

def createIsoObject(pos,qpos,isoClass, varCount, ref, alt, iso, count, fillStr):
    varCount=varCount+1
    qual=-1
    alignInfo={'queryStart':iso['q_st'],'queryEnd':iso['q_end'],'qLen':iso['query_length'],'subjectStart':iso['s_st'],'subjectEnd':iso['s_end'],'sLen':iso['hit_length'],'qSerialNo':count}
    iso=AminoAcidVariation(iso['hit_def'].replace(';',''),iso['query_name'],pos,qpos,ref,alt,isoClass,iso['Location'],varCount, qual, fillStr, alignInfo)
    return (varCount,count, iso)

def findORFType(header):
    orfType=''
    if 'type:complete' in header:
        orfType='_complete'
    elif 'type:3prime_partial' in header:
        orfType='_C-terminus_partial'
    elif 'type:5prime_partial' in header:
        orfType='_N-terminus_partial'
    elif 'type:internal' in header:
        orfType='_internal'
    else:
        print("Error: ORF type was not found. This should not happen if the ORFs were predicted using transdecoder.")
    return orfType
def printValidatedVariations(variation, newVcfFile, headerFlag, fileFlag):
    ##This is another vcf
    if headerFlag==1:
        newVcfFile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    ##SubjectID','QueryID','Alignment','Type','QPOS'
    if fileFlag==1:
        newVcfFile.write(str(variation['#CHROM'])+"\t"+str(variation['POS'])+"\t"+str(variation['ID'])+"\t"+variation['REF']+"\t"+variation['ALT']+"\t"+str(variation['QUAL'])+"\t"+variation['FILTER']+"\t"+"SubjectId="+variation['SubjectID']+";QueryId="+variation['QueryID'].replace(',','&')+";QueryLength="+variation['QueryLength']+";QueryStart="+variation['QueryStart']+";QueryEnd="+variation['QueryEnd']+";SubjectLength="+variation['SubjectLength']+";SubjectStart="+variation['SubjectStart']+";SubjectEnd="+variation['SubjectEnd']+";Type="+variation['Type']+";QPOS="+str(variation['QPOS'])+";PeptideCount="+str(variation['PeptideCount'])+";UniquePeptideCount="+str(variation['UniquePeptideCount'])+";Peptides="+variation['UniquePeptide']+";Evidence="+variation['Evidence']+";Score="+str(variation['ConfidenceScore'])+"\n")
        #newVcfFile.write(str(variation['#CHROM'])+"\t"+str(variation['POS'])+"\t"+str(variation['ID'])+"\t"+variation['REF']+"\t"+variation['ALT']+"\t"+str(variation['QUAL'])+"\t"+variation['FILTER']+"\t"+"SubjectId="+variation['SubjectID']+";QueryId="+variation['QueryID'].replace(',','&')+";QueryLength="+variation['QueryLength']+";QueryStart="+variation['QueryStart']+";QueryEnd="+variation['QueryEnd']+";SubjectLength="+variation['SubjectLength']+";SubjectStart="+variation['SubjectStart']+";SubjectEnd="+variation['SubjectEnd']+";Type="+variation['Type']+";QPOS="+str(variation['QPOS'])+";PeptideCount="+str(variation['PeptideCount'])+";UniquePeptideCount="+str(variation['UniquePeptideCount'])+";Peptides="+variation['UniquePeptide']+";Score="+str(variation['ConfidenceScore'])+"\n")
    else:
        newVcfFile.write(str(variation['#CHROM'])+"\t"+str(variation['POS'])+"\t"+str(variation['ID'])+"\t"+variation['REF']+"\t"+variation['ALT']+"\t"+str(variation['QUAL'])+"\t"+variation['FILTER']+"\t"+"SubjectId="+variation['SubjectID']+";QueryId="+variation['QueryID'].replace(',','&')+";Alignment=["+variation['Alignment']+"];Type="+variation['Type']+";QPOS="+str(variation['QPOS'])+"\n")

def IsoformVariation(identifiedBlast,count):
    ##we dont want to store the sequence, hence alt and ref are both '.'
    alt="."
    ref="."
    fillStr="TRANS"
    ##positions need redoing.
    isoList=[]
    for i,iso in identifiedBlast.iterrows():
        #print(iso['query_name'])
        #print(int(iso['q_st']))
        iso['q_st']=int(iso['q_st'])
        iso['s_st']=int(iso['s_st'])
        iso['q_end']=int(iso['q_end'])
        iso['s_end']=int(iso['s_end'])
        iso['query_length']=int(iso['query_length'])
        iso['hit_length']=int(iso['hit_length'])
        varCount=0
        count=count+1
        orfType=findORFType(iso['query_name'])
        if iso['q_st']==1:
            if iso['s_st']==1:
                if iso['q_end']==iso['query_length']:
                    if iso['s_end']==iso['hit_length']:
                        print("This is a case of known or known with variation protein. This ORF should not be in the list of isoforms.")
                    elif iso['s_end']<iso['hit_length']:
                        #The protein is partially mapped to the ORF. Because the q_end==query_length, we can say that the protein is longer than the ORF.
                        #The ORF is shorter than the protein at the 3' end.
                        isoClass='3prime_shortened'+orfType   
                        varCount, count, isoObj = createIsoObject(iso['s_end'],iso['q_end'],isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                elif iso['q_end']<iso['query_length']:
                    if iso['s_end']==iso['hit_length']:
                        ##The ORF is larger than the protein, because the full length of the protein has a map  to the ORF.
                        isoClass='3prime_extended'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end'],iso['q_end'],isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        ##Both the ORF and the protein has sequence at the 3' end that does not map to each other. Might be result of alternative
                        ## end exon.
                        isoClass='3prime_alternative'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end'],iso['q_end'],isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                else:
                    print('queary_end can not be bigger than the query length')
            elif iso['s_st']>1:
                ##protein pos is sent as 1 because the AminoAcid variation print function deals with iso['s_st']>1 ans add the value as prefix.
                if iso['q_end']==iso['query_length']:
                    if iso['s_end']==iso['hit_length']:
                        ##protein's 5' side has extra amino acids. i.e. this is a deletion event at 5'
                        isoClass='5prime_shortened'+orfType
                        varCount, count, isoObj=createIsoObject(1, iso['q_st'], isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        ##part of the protein match with the ORF. But the whole ORF is mapping to the protein. So two events.
                        isoClass='5prime_shortened'+orfType   
                        varCount, count, isoObj=createIsoObject(1,iso['q_st'],isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                        isoClass='3prime_shortened'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end']-iso['s_st']+1, iso['q_end'], isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                elif iso['q_end']<iso['query_length']:
                    #q_end<quety_length, i.e. the whole ORF did not match.
                    if iso['s_end']==iso['hit_length']:
                        ##this is an event of INS at the 3' end. And also an event of deletion at the 5' end. because s_st>1
                        isoClass='5prime_shortened'+orfType   
                        varCount, count, isoObj=createIsoObject(1,iso['q_st'],isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                        isoClass='3prime_extended'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end']-iso['s_st']+1,iso['q_end'],isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        ##Both the ORF and the protein has sequence at the 3' end that does not map to each other. Might be result of alternative
                        ## end exon.
                        isoClass='5prime_shortened'+orfType   
                        varCount, count, isoObj=createIsoObject(1,iso['q_st'],isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                        orfType=findORFType(iso['query_name'])
                        isoClass='3prime_alternative'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end']-iso['s_st']+1,iso['q_end'],isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                else:
                    print('queary_end can not be bigger than the query length')
        elif iso['q_st']>1:
            if iso['s_st']==1:
                ##Insertion
                if iso['q_end']==iso['query_length']:
                    #3' side ORF complete match
                    if iso['s_end']==iso['hit_length']:
                        #3' protein full match
                        isoClass='5prime_extended'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_st'],1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        #3' protein short match
                        #so event of 5'INS and 3'del
                        isoClass='5prime_extended'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_st'],1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                        isoClass='3prime_shortened'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end'],iso['q_end']-iso['q_st']+1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                elif iso['q_end']<iso['query_length']:
                    #3' side ORF short match
                    if iso['s_end']==iso['hit_length']:
                        #3' protein full match
                        #as q_st>1 and q_end<q_length, we have insertion at both the end.
                        isoClass='5prime_extended'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_st'],1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                        isoClass='3prime_extended'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end'],iso['q_end']-iso['q_st']+1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        #3' protein short match
                        #so event of 5'INS and 3' ALT
                        isoClass='5prime_extended'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_st'],1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                        isoClass='3prime_alternative'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end'],iso['q_end']-iso['q_st']+1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                else:
                    print("Error: query_end can not be bigger than the query length")
            elif iso['s_st']>1:
                ##check Updated code. Here both query and subject 5 prime entries shouls have position value 1.
                ##Alteration at 5 prime.
                if iso['q_end']==iso['query_length']:
                    #3' side ORF complete match
                    if iso['s_end']==iso['hit_length']:
                        #3' protein full match
                        isoClass='5prime_alternative'+orfType   
                        varCount, count, isoObj=createIsoObject(1,1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        #3' protein short match
                        #so event of 5'INS and 3'del
                        isoClass='5prime_alternative'+orfType   
                        varCount, count, isoObj=createIsoObject(1,1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                        isoClass='3prime_shortened'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end']-iso['s_st']+1,iso['q_end']-iso['q_st']+1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                elif iso['q_end']<iso['query_length']:
                    #3' side ORF short match
                    if iso['s_end']==iso['hit_length']:
                        #3' protein full match
                        #as q_st>1 and q_end<q_length, we have insertion at both the end.
                        isoClass='5prime_alternative'+orfType   
                        varCount, count, isoObj=createIsoObject(1,1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                        isoClass='3prime_extended'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end']-iso['s_st']+1,iso['q_end']-iso['q_st']+1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        #3' protein short match
                        #so event of 5'INS and 3' ALT
                        isoClass='5prime_alternative'+orfType   
                        varCount, count, isoObj=createIsoObject(1,1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                        isoClass='3prime_alternative'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end']-iso['s_st']+1,iso['q_end']-iso['q_st']+1,isoClass, varCount, ref, alt, iso, count, fillStr)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                else:
                    print("Error: query_end can not be bigger than the query length")
                
    return (count, isoList)


def checkIsoformPeptideCoverage(ORFId, PSMsOfORF, variation):
    prevPeptide=''
    prevFound=0
    PSMEvidence=pd.DataFrame();
    print("variation id:"+str(variation['ID']))
    for i in range(0,len(PSMsOfORF)):
        protAcc=PSMsOfORF.iloc[i]['proteinacc_start_stop_pre_post_;']
        print("prot accession:"+protAcc)
        pepSeq=PSMsOfORF.iloc[i]['Sequence']
        if prevPeptide!=pepSeq:
            ## Check whether this PSM covers the location of the variation.
            ## If yes, compare the peptide sequence and the variation sequence to make sure the peptide does
            ## support the variation.
            ## Split the protId, get start and end and compare with the loc. Check the sequence, compare with altSeq.
            prevPeptide=pepSeq;
            ## prevFound flag is reset here. Because prevPeptide is not same as current one, therefore the PSM has not been
            ## added to the PSMRvidence yet.
            prevFound=0
            protList=protAcc.split(';')
            proacc_start_stop_pre_post=[s for s in protList if ORFId in s] #this should always have one element.[Wrong, one peptide might map to multiple location of an ORF]
            if len(proacc_start_stop_pre_post)>0: ##changed this so that multiple location of the peptide is considered.
                for j in range(0, len(proacc_start_stop_pre_post)):                
                    ##this should always be the case.
                    ## 'start' and 'stop' are both inclusive
                    startEndPrePost=proacc_start_stop_pre_post[j].replace(ORFId,'')
                    print("start end pre post:"+startEndPrePost)
                    pattern=re.compile('_(\d+)_(\d+)_([A-Z]|-)_([A-Z]|-)')
                    match=pattern.finditer(startEndPrePost)
                    count=0
                    ##check following code carefully and make necessary changes to make it work for terminal INDELs and ALT. Also figure out a way to represent splice junction peptides.
                    for m in match:
                        ##this loop should run only once.
                        proacc_start_stop_pre_postList=m
                        if count>0:
                            ##ERROR
                            print("ERROR: Match object should have single value."+protAcc)
                            print("ORFId:"+ORFId)
                            print("proacc_start_stop_pre_postList:"+proacc_start_stop_pre_postList.groups())
                            count=count+1;
                            break;
                        count=count+1;
                    ## this is what should happen
                    #print("Count:"+str(count))
                    if count==1:
                        ## Check whether location of the ALT/INDEL is covered by this peptide.
                        ## For SAPs its straight forward. 
                        ## It might happen that for ALT/INDELs the peptide covers the event partially.
                        #print("Count:"+str(count))
                        if '5prime' in variation['Type'] and 'shortened' not in variation['Type']:
                            altStart=1
                            #print("ORFId:"+ORFId)
                            #print("variation:"+variation.to_csv(None))
                            #print("alt:"+str(variation['ALT']))
                            #Minus 1 because the blast start and end is inclusive.
                            altEnd=int(variation['QPOS'])-1
                            ## Check whether this identified peptide overlaps with this alteration event.
                            ## The overlap can either be full or partial overlap.
                            if int(proacc_start_stop_pre_postList.group(1))>=altStart and altEnd>=int(proacc_start_stop_pre_postList.group(2)) and int(proacc_start_stop_pre_postList.group(1))<=altEnd:
                                ##this means the alteration event is fully covered by the peptide.
                                ##pythonic indexes are 0 based.
                                position=int(variation['QPOS'])-int(proacc_start_stop_pre_postList.group(1))
                                print("2. Overlap")
                                prevFound=1
                                PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'internal'})),ignore_index=True)
                                print("2.pep found")
                                
                            elif int(proacc_start_stop_pre_postList.group(1))>=altStart and altEnd>=int(proacc_start_stop_pre_postList.group(1)) and altEnd<=int(proacc_start_stop_pre_postList.group(2)):
                                ## juntion peptide
                                position=0
                                positionEnd=altEnd-int(proacc_start_stop_pre_postList.group(1))+1
                                altIdx=altStart-int(proacc_start_stop_pre_postList.group(1))
                                print("3. Overlap")
                                prevFound=1
                                PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'junction'})),ignore_index=True)
                                print("3.pep found")
                            elif int(proacc_start_stop_pre_postList.group(1))<=altStart:
                                ##this is not possible for 5prime casses
                                print("ERROR: 5prime variation events can not have a peptide which starts before the alt start.")
                            else:
                                print(str(variation['ID'])+": This variation is not supported by "+PSMsOfORF.iloc[i]['Sequence'])
                        elif "3prime" in variation['Type'] and 'shortened' not in variation['Type']:
                            altStart=int(variation['QPOS'])+1
                            altEnd=int(variation['QueryLength'])
                            ## Check whether this identified peptide overlaps with this insertion event.
                            ## The overlap can either be full or partial overlap.
                            if int(proacc_start_stop_pre_postList.group(1))<=altStart and altEnd>=int(proacc_start_stop_pre_postList.group(2)) and int(proacc_start_stop_pre_postList.group(2))>=altStart:
                                ##junction peptide
                                ##pythonic indexes are 0 based.
                                position=int(variation['QPOS'])-int(proacc_start_stop_pre_postList.group(1))
                                altLength=len(variation['ALT'])
                                print("5. Overlap")
                                prevFound=1
                                PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'junction'})),ignore_index=True)
                            elif int(proacc_start_stop_pre_postList.group(1))>=altStart and altEnd>=int(proacc_start_stop_pre_postList.group(1)) and altEnd>=int(proacc_start_stop_pre_postList.group(2)):
                                ## internal coverage
                                position=0
                                positionEnd=altEnd-int(proacc_start_stop_pre_postList.group(1))+1
                                ## part of the variation['ALT'] will be missing
                                altIdx=altStart-int(proacc_start_stop_pre_postList.group(1))
                                print("6. Overlap")
                                prevFound=1
                                PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'internal'})),ignore_index=True)
                                print("6.pep found")
                            else:
                                print(str(variation['ID'])+": This variation is not supported by "+PSMsOfORF.iloc[i]['Sequence'])
                    elif count==0:
                        print("Count is 0, i.e. protacc_start_end_prev_next does not match:"+proacc_start_stop_pre_post[j])
                        print("ORFId:"+ORFId)
                    else:
                        print("ERROR: count is not 1:"+protAcc)
                        print("ORFId:"+ORFId)
                        print("proacc_start_stop_pre_postList:"+proacc_start_stop_pre_postList.groups())
            else:
                ##ERROR
                print(proacc_start_stop_pre_post+" ERROR: proacc_start_stop_pre_post list should have at least one entry")
        else:
            ## This suggest multiple Spectra match to a peptide. Peptide sequence might be the same, but the PSM
            ## is not.
            if prevFound==1:
                print("11. Overlap")
                print("10.pep found")
                print("Length of PSMEvidence"+str(len(PSMEvidence)))
                ##this means prevPeptide was counted as an evidence of the variaton. Hence this PSM should also
                ##be counted.
                print(PSMsOfORF.iloc[i].to_csv(None))
                PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':PSMEvidence.iloc[len(PSMEvidence)-1]['Evidence']})),ignore_index=True)
            #else:
            #    print(str(variation['ID'])+": This variation is not supported by "+PSMsOfORF.iloc[i]['Sequence'])
    print("PSMEvd:"+PSMEvidence.to_csv())
    return PSMEvidence

def groupPeptide(psms, ORFId, pepColumn, rev):
    ##not considering CONT because CONT entries are not going to have ORFId.
    regex="(?<!"+rev+")"+re.escape(ORFId)+"_" ##"_" is added to make sure that id asmbl_27582|m.351368 is not picked when ORFSID is asmbl_27582|m.35136
    filteredPSMs= psms[psms[pepColumn].str.contains(regex)]
    return filteredPSMs

def groupVariation(variationDF, ORFId, rev):
    ##not considering CONT because CONT entries are not going to have ORFId.
    regex="(?<!"+rev+")"+re.escape(ORFId)+" "
    filteredVar= variationDF[variationDF['QueryID'].str.contains(regex)]
    return filteredVar

def printVariationListtoFile(varList, fileWrt):
    print("writing to file")
    for l in varList:
        fileWrt.write(l.toString()+"\n");
def calculateQualityScore(PSMsOfVariation,featureWeights):
    ## For time being this function returns 1, but in future I need to figure out a formula to represent the confidence of peptide evidences supporting a variation.
    return 1

def addPSMInfoToVariation(variation, PSMsOfVariation):
    ##Add PSM count, Peptide count, possibly an average score indicating PSM/peptide quality.
    ##For time being quality score is 1.
    score=calculateQualityScore(PSMsOfVariation,1)
    ## Number of PSMs
    if len(PSMsOfVariation)>0:
        PSMCount=len(PSMsOfVariation['PSM_ID'].unique())
        ## Number of peptide sequences
        peptideCount=len(PSMsOfVariation['Sequence'].unique())
        ## Count of peptide sequences that uniquely map to this ORF/Protein
        unqPeptideCount=len(PSMsOfVariation[(~PSMsOfVariation['proteinacc_start_stop_pre_post_;'].str.contains(';'))]['Sequence'].unique()) ##~ is doing the job of negation, i.e. does not contain. This means the peptide is not shared with other proteins
        evdType=",".join(PSMsOfVariation['Evidence'].unique().tolist())
        ## Unique peptides
        print("Unique Peptide:")
        print(PSMsOfVariation[(~PSMsOfVariation['proteinacc_start_stop_pre_post_;'].str.contains(';'))]['Sequence'].unique())
        unqPeptides=",".join(PSMsOfVariation[(~PSMsOfVariation['proteinacc_start_stop_pre_post_;'].str.contains(';'))]['Sequence'].unique())
        if not unqPeptides:
            unqPeptides="-"
        PSMInfo={'PSMCount':PSMCount,'PeptideCount':peptideCount,'UniquePeptideCount':unqPeptideCount,'UniquePeptide':unqPeptides,'Evidence':evdType,'ConfidenceScore':score}
    else:
        PSMInfo={'PSMCount':0,'PeptideCount':0,'UniquePeptideCount':0,'UniquePeptide':'-','Evidence':'-','ConfidenceScore':0}
    writeVar=variation.append(pd.Series(PSMInfo))
    print("writeVar")
    print(writeVar.to_csv(None))
    return writeVar

def printPSMsOfVariation(writeDF, newPSMFile, headerFlag):
    ##this is PSM wise view
    if headerFlag==1:
        ## Write column names
        writeDF.to_csv(newPSMFile, index=False)
    else:
        writeDF.to_csv(newPSMFile, header=False, index=False)

def addVariationInfoToPSM(PSMsOfVariations, variation):
    ## Add ORFId and variation information.
    varInf={'ORFID':variation['QueryID'],'SubjectID':variation['SubjectID'],'variationID':variation['ID'],'variationType':variation['Type'],'Location':variation['QPOS']}
    lenPSMVar=len(PSMsOfVariations)
    writeDF=pd.concat([PSMsOfVariations,pd.DataFrame([varInf]*lenPSMVar,index=range(lenPSMVar))],axis=1)
    return writeDF
##isoList is a list of AminoAcidVariation class.
def peptidePerIsoform(peptides, variationDF, isoIds, pepColumn, rev, newVcfFile,):
    print("Peptide Count:")
    ##Known proteins
    psms=filterPeptide(peptides, isoIds, pepColumn)
    #psms.to_csv(,index=False)
    headerFlag=1
    print("Peptides:"+str(len(psms['Sequence'].unique())))
    for i in range(0,len(isoIds)):
        PSMsOfORF=groupPeptide(psms, isoIds[i], pepColumn, rev)
        variations=groupVariation(variationDF, isoIds[i], rev)
        for j in range(0,variations.shape[0]):
            variation=variations.iloc[j]
            PSMsOfVariations=checkIsoformPeptideCoverage(isoIds[i], PSMsOfORF, variation)
            if len(PSMsOfVariations)>0:
                print("Variation "+str(variation['ID'])+" has peptide evidence")
                PSMsOfVariations=addVariationInfoToPSM(PSMsOfVariations, variation)
                variation=addPSMInfoToVariation(variation, PSMsOfVariations)
                variation['QUAL']=PSMsOfVariations['PSM-level q-value'].mean()
                variation['FILTER']='PASS'
                #printPSMsOfVariation(PSMsOfVariations, newPSMFile,headerFlag)
                printValidatedVariations(variation, newVcfFile,headerFlag,1)
                headerFlag=0
            else:
                print(str(variation['ID'])+" does not have peptide evidence")
                ##still write this variation in the vcffile.
                variation=addPSMInfoToVariation(variation, PSMsOfVariations)
                printValidatedVariations(variation, newVcfFile,headerFlag,1)
                headerFlag=0
    
protRevStr="XXX_"
protContStr="CONT_"
pepTh=1
pepColumn='proteinacc_start_stop_pre_post_;'
parser = argparse.ArgumentParser(description='This script read output of UniProteinLocation.py, SplitAnnotationFile.py and finds peptide evidence for protein isoforms')
parser.add_argument("-b", "--blast", nargs=1, required=True, help="full path of blast csv file", metavar="PATH")
parser.add_argument("-a", "--annotation", nargs=1, required=True, help="full path to the annotation file", metavar="PATH")
#parser.add_argument("-j", "--isovar", nargs=1, required=True, help="full path to the isoforms with variation file", metavar="PATH")
parser.add_argument("-p", "--pep", nargs=1, required=True, help="full path to the peptide identification file", metavar="PATH")
parser.add_argument("-o", "--output", nargs=1, required=True, help="full path to the output csv file", metavar="PATH")
args = parser.parse_args()
newVcfFile=args.output[0]
##Read annotation file
annotations=readFile(args.annotation[0],',')
print("Total number of TGE obs:"+str(annotations.shape[0]))

isoforms=annotations.loc[annotations['Class'].str.contains("prime_")]
isoIds=isoforms['ORF Id'].str.extract("([^ ]*)").tolist()
print("Total number of isoform TGE obs:"+str(isoforms.shape[0]))
### SEparates blast results for the isoform list
identifiedIsoBlast=filterBlast(args.blast[0], isoforms['ORF Id'])

##Identified PSMs
peptides=readFile(args.pep[0],',')
count=0
##isoList is a list of AminoAcid  class
count, isoList=IsoformVariation(identifiedIsoBlast,count)
with open(newVcfFile, 'w') as vcfFile:
    vcfFile.write(AminoAcidVariation.printHeader()+"\n");
    printVariationListtoFile(isoList,vcfFile)

vcf=readFile(newVcfFile,'\t')

info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','QueryLength','QueryStart','QueryEnd','SubjectLength','SubjectStart','SubjectEnd','Type','QPOS'])
#SubjectId=sp|Q92466-5|DDB2_HUMAN;QueryId=asmbl_10851|m.134573;QueryLength=244;QueryStart=1;QueryEnd=237;SubjectLength=244;SubjectStart=1;SubjectEnd=237;Type=3prime_alternative_C-terminus_partial;QPOS=237

vcf=vcf.drop('INFO',1)
vcf=vcf.join(info)
vcf.SubjectID=vcf.SubjectID.str.replace('SubjectId=','')
vcf.QueryID=vcf.QueryID.str.replace('QueryId=','')
#print(vcf.QueryID[0:5])
vcf.QueryLength=vcf.QueryLength.str.replace('QueryLength=','')
vcf.QueryStart=vcf.QueryStart.str.replace('QueryStart=','')
vcf.QueryEnd=vcf.QueryEnd.str.replace('QueryEnd=','')
vcf.SubjectLength=vcf.SubjectLength.str.replace('SubjectLength=','')
vcf.SubjectStart=vcf.SubjectStart.str.replace('SubjectStart=','')
vcf.SubjectEnd=vcf.SubjectEnd.str.replace('SubjectEnd=','')
vcf.Type=vcf.Type.str.replace('Type=','')
vcf.QPOS=vcf.QPOS.str.replace('QPOS=','')    

print(vcf.columns.values)
newVcfFile2=args.annotation[0]+"_isoform_pep.vcf"
print("new file:"+newVcfFile2)
with open(newVcfFile2,'w') as isoPep:
    peptidePerIsoform(peptides, vcf, isoIds, pepColumn, protRevStr, isoPep)


print("\n\n\n TEST")
#for i in range(0,len(variationList)):
#    print(variationList[i].toString())

