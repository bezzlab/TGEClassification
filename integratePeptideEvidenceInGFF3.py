##This code integrates GFF3 file and protein identification file for genoverse.
##Input is *.assemblies.fasta.transdecoder.genome.gff3 and *+fdr+th+grouping+prt_filtered.csv
import pandas as pd
import re
import os
import argparse

def readFile(filename, sep, headerFlag):
    if headerFlag==0:
        fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    elif headerFlag==1:
        fileDFObj = pd.read_table(filename, sep=sep, header=None, keep_default_na=False, na_values=[''])
    else:
        print("Unrecognized Header Flag")
    return fileDFObj;

def nucleotideCount(gffObj):
    ##this function counts total nucleotide in a feature and add that as an extra column.
    gffObjCount=pd.concat([gffObj,pd.Series([0]*gffObj.shape[0],name='nCount',index=range(gffObj.shape[0]))],axis=1)
    gffObjCount['nCount']=gffObjCount['end']-gffObjCount['start']+1
    return gffObjCount

def cdsAggregate(prtGffCDS):
    print(prtGffCDS.shape)
    print(prtGffCDS.iloc[0]['strand'])
    aggreNCount=pd.Series([0]*prtGffCDS.shape[0],name='aggreNCount',index=prtGffCDS.index)
    if prtGffCDS.iloc[0]['strand']=="+":
        ##The gene/mRNA is on positive strand, hence nucleotide aggregation of CDS is done from the first CDS in the dataframe
        idx=prtGffCDS.index.tolist()
        for i in range(0,len(idx)):
            if i==0:
                aggreNCount.loc[idx[i]]=prtGffCDS.loc[idx[i]]['nCount']
                print(aggreNCount.loc[idx[i]])
            else:
                aggreNCount.loc[idx[i]]=prtGffCDS.loc[idx[i]]['nCount']+aggreNCount.loc[idx[i-1]]
    else:
        print("negetive strand")
        ##This mRNA is on negetive strand. Hence, nucleotide aggregation starts from last CDS.
        idx=prtGffCDS.index.tolist()
        for i in range(len(idx)-1,-1,-1):
            #print(i)
            if i==len(idx)-1:
                aggreNCount.loc[idx[i]]=prtGffCDS.loc[idx[i]]['nCount']
                print(aggreNCount.loc[idx[i]])
            else:
                aggreNCount.loc[idx[i]]=prtGffCDS.loc[idx[i]]['nCount']+aggreNCount.loc[idx[i+1]]
    return aggreNCount

def peptidePosition(subGffObj, prtAcc, start_stop_pepqval_mod_charge_psmcount):
    ##this function take a subset of gffObj [containing only CDS of a mRNA with aggregated nucleotide sum], protein/orf id and peptide location and calculate genomic location for the peptide.
    ##I subtract 2 because if the peptide match from first AA, then start is 1st nucleotide. in the same way, if start is 5,  nucleotide
    ##start is 4*3+1 i.e. 1st nucleotide of 5th AA.
    ## This function is called for in a loop for all peptides mapping to a gene.
    pepGffEntries=pd.DataFrame(columns=subGffObj.columns[:-2])
    for i in range(0,len(start_stop_pepqval_mod_charge_psmcount)):
        nucStart=int(start_stop_pepqval_mod_charge_psmcount[i][0])*3-2
        nucEnd=int(start_stop_pepqval_mod_charge_psmcount[i][1])*3
        print("nucStart - nunEnd")
        print(nucStart)
        print(nucEnd)
        overlappingCDS1=subGffObj[subGffObj['aggreNCount']>=nucStart]
        overlappingCDS2=subGffObj[subGffObj['aggreNCount']>=nucEnd]
        print("dim:"+str(overlappingCDS1.shape))
        #print(overlappingCDS1)
        #print(overlappingCDS2)
        if overlappingCDS1.iloc[0]['strand']=="+":
            #positive strand
            print("Positive Strand")
            print("CDS1-index")
            print(overlappingCDS1.index)
            print("CDS1-aggreNCount")
            print(overlappingCDS1['aggreNCount'])
            print("CDS2-index")
            print(overlappingCDS2.index)
            print("CDS2-aggreNCount")
            print(overlappingCDS2['aggreNCount'])
            #CDS that contains start of the peptide
            idx1=min(overlappingCDS1.index)
            idx2=min(overlappingCDS2.index)
            pepStartCDS=overlappingCDS1.loc[min(overlappingCDS1.index)]
            pepStopCDS=overlappingCDS2.loc[min(overlappingCDS2.index)]
            print(str(pepStartCDS['start'])+","+str(pepStopCDS['start']))
            if min(overlappingCDS1.index) == min(overlappingCDS2.index):
                print("single CDS peptide")
                pepStart=pepStartCDS['end'] - (pepStartCDS['aggreNCount']-nucStart)
                pepStop=pepStartCDS['end'] - (pepStartCDS['aggreNCount']-nucEnd)
                temp=pd.DataFrame([subGffObj.loc[idx1][['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase','attributes']]])
                temp['type']='peptide'
                temp['start']=pepStart
                temp['end']=pepStop
                cdsid=temp['attributes'].str.extract('ID=([^;]+)').tolist()[0]
                #temp['attributes']='ID=pep'+str(i)+'.'+cdsid+";Parent="+cdsid #+";Gap=M"+str(nucEnd-nucStart+1)
                mod='None'
                if start_stop_pepqval_mod_charge_psmcount[i][3]!='':
                	mod=start_stop_pepqval_mod_charge_psmcount[i][3]
                
                temp['attributes']='ID=pep'+str(i)+'.'+cdsid+";Parent="+cdsid+";Q-Value="+str(start_stop_pepqval_mod_charge_psmcount[i][2])+";Modifications="+mod+";Charge="+str(start_stop_pepqval_mod_charge_psmcount[i][4])+";PSM-Count="+str(start_stop_pepqval_mod_charge_psmcount[i][5])
                
                pepGffEntries=pepGffEntries.append(temp)
                print("Peptide Start-End:"+str(pepStart)+"-"+str(pepStop))
            else:
                print("Multiple CDS peptide:"+str(i))
                sIdx=subGffObj.index.tolist()
                cdsIdx=[s for s in sIdx if min(overlappingCDS2.index)>=s>=min(overlappingCDS1.index)]
                for j in range(0,len(cdsIdx)):
                    pepStart=None
                    pepStop=None
                    if j==0:
                        pepStart=subGffObj.loc[cdsIdx[j]]['end'] - (subGffObj.loc[cdsIdx[j]]['aggreNCount']-nucStart)
                        pepStop=subGffObj.loc[cdsIdx[j]]['end']
                    elif j==len(cdsIdx)-1:
                        pepStart=subGffObj.loc[cdsIdx[j]]['start']
                        pepStop=subGffObj.loc[cdsIdx[j]]['end'] - (subGffObj.loc[cdsIdx[j]]['aggreNCount']-nucEnd)
                    else:
                        pepStart=subGffObj.loc[cdsIdx[j]]['start']
                        pepStop=subGffObj.loc[cdsIdx[j]]['end']
                    temp=pd.DataFrame([subGffObj.loc[idx1][['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase','attributes']]])
                    temp['type']='peptide'
                    temp['start']=pepStart
                    temp['end']=pepStop
                    print("Peptide start-stop:"+str(pepStart)+"-"+str(pepStop))
                    cdsid=temp['attributes'].str.extract('ID=([^;]+)').tolist()[0]
                    #temp['attributes']='ID=pep'+str(i)+'_'+str(j)+'.'+cdsid+";Parent="+cdsid #+";Gap=M"+str(nucEnd-nucStart+1)
                    mod='None'
                    if start_stop_pepqval_mod_charge_psmcount[i][3]!='':
                    	mod=start_stop_pepqval_mod_charge_psmcount[i][3]
                    temp['attributes']='ID=pep'+str(i)+'_'+str(j)+'.'+cdsid+";Parent="+cdsid+";Q-Value="+str(start_stop_pepqval_mod_charge_psmcount[i][2])+";Modifications="+mod+";Charge="+str(start_stop_pepqval_mod_charge_psmcount[i][4])+";PSM-Count="+str(start_stop_pepqval_mod_charge_psmcount[i][5])
                    pepGffEntries=pepGffEntries.append(temp)
        else:
            #negetive strand
            print("Negetive Strand")
            print("CDS1-index")
            print(overlappingCDS1.index)
            print("CDS1-aggreNCount")
            print(overlappingCDS1['aggreNCount'])
            print("CDS2-index")
            print(overlappingCDS2.index)
            print("CDS2-aggreNCount")
            print(overlappingCDS2['aggreNCount'])
            #CDS that contains start of the peptide
            idx1=max(overlappingCDS1.index)
            idx2=max(overlappingCDS2.index)
            pepStartCDS=overlappingCDS1.loc[idx1]
            pepStopCDS=overlappingCDS2.loc[idx2]
            print(str(pepStartCDS['start'])+","+str(pepStopCDS['end']))
            if idx1 == idx2:
                print("single CDS peptide")
                pepStart=pepStartCDS['end'] - (pepStartCDS['aggreNCount']-nucStart)
                pepStop=pepStartCDS['end'] - (pepStartCDS['aggreNCount']-nucEnd)
                temp=pd.DataFrame([subGffObj.loc[idx1][['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase','attributes']]])
                temp['type']='peptide'
                temp['start']=pepStart
                temp['end']=pepStop
                cdsid=temp['attributes'].str.extract('ID=([^;]+)').tolist()[0]
                #temp['attributes']='ID=pep'+str(i)+'.'+cdsid+";Parent="+cdsid #+";Gap=M"+str(nucEnd-nucStart+1)
                mod='None'
                if start_stop_pepqval_mod_charge_psmcount[i][3]!='':
                	mod=start_stop_pepqval_mod_charge_psmcount[i][3]
                
                temp['attributes']='ID=pep'+str(i)+'.'+cdsid+";Parent="+cdsid+";Q-Value="+str(start_stop_pepqval_mod_charge_psmcount[i][2])+";Modifications="+mod+";Charge="+str(start_stop_pepqval_mod_charge_psmcount[i][4])+";PSM-Count="+str(start_stop_pepqval_mod_charge_psmcount[i][5])
                pepGffEntries=pepGffEntries.append(temp)
                print("Peptide Start-End:"+str(pepStart)+"-"+str(pepStop))
            else:
                print("Multiple CDS peptide"+str(i))
                sIdx=subGffObj.index.tolist()
                cdsIdx=[s for s in sIdx if max(overlappingCDS2.index)<=s<=max(overlappingCDS1.index)]
                for j in range(len(cdsIdx)-1,-1,-1):
                    pepStart=None
                    pepStop=None
                    if j==0:
                        pepStart=subGffObj.loc[cdsIdx[j]]['start']
                        pepStop=subGffObj.loc[cdsIdx[j]]['end'] - (subGffObj.loc[cdsIdx[j]]['aggreNCount']-nucEnd)
                    elif j==len(cdsIdx)-1:
                        pepStart=subGffObj.loc[cdsIdx[j]]['end'] - (subGffObj.loc[cdsIdx[j]]['aggreNCount']-nucStart)
                        pepStop=subGffObj.loc[cdsIdx[j]]['end']
                    else:
                        pepStart=subGffObj.loc[cdsIdx[j]]['start']
                        pepStop=subGffObj.loc[cdsIdx[j]]['end']
                    temp=pd.DataFrame([subGffObj.loc[idx1][['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase','attributes']]])
                    temp['type']='peptide'
                    temp['start']=pepStart
                    temp['end']=pepStop
                    print("Peptide start-stop:"+str(pepStart)+"-"+str(pepStop))
                    cdsid=temp['attributes'].str.extract('ID=([^;]+)').tolist()[0]
                    #temp['attributes']='ID=pep'+str(i)+'_'+str(j)+'.'+cdsid+";Parent="+cdsid #+";Gap=M"+str(nucEnd-nucStart+1)
                    mod='None'
                    if start_stop_pepqval_mod_charge_psmcount[i][3]!='':
                        mod=start_stop_pepqval_mod_charge_psmcount[i][3]
                    temp['attributes']='ID=pep'+str(i)+'_'+str(j)+'.'+cdsid+";Parent="+cdsid+";Q-Value="+str(start_stop_pepqval_mod_charge_psmcount[i][2])+";Modifications="+mod+";Charge="+str(start_stop_pepqval_mod_charge_psmcount[i][4])+";PSM-Count="+str(start_stop_pepqval_mod_charge_psmcount[i][5])
                    pepGffEntries=pepGffEntries.append(temp)
    return pepGffEntries
    

def findPeptidesStartStop(prt, pepObj):
    ##This function finds all the peptides for a protein id, format the 'prtacc_start_stop_prev_post;' string and extract the start and stop AAs.
    
    #Keep only ['Sequence' and 'prtacc_start_stop_prev_post;' columns]
    #pepObjUnique=pepObj[pepObj['proteinacc_start_stop_pre_post_;'].str.contains(re.escape(prt)+"_")][['Sequence','proteinacc_start_stop_pre_post_;']].drop_duplicates()
	aggregations = {
    ##It used to be MS-GF:PepQVal, which is missing if tda=0 is used for MSGF+ search.
	'PSM-level q-value': lambda x: x.mean(),
	'Modifications': lambda x: ','.join(x.dropna().unique()),
	'Charge': lambda x: ','.join([str(s) for s in x.dropna().unique()]),
	'Spectrum ID': lambda x: x.shape[0]
	}
	pepObjUnique = pepObj[pepObj['proteinacc_start_stop_pre_post_;'].str.contains(re.escape(prt)+"_")].groupby(['Sequence','proteinacc_start_stop_pre_post_;']).agg(aggregations).reset_index()
	#extract 'prt''s start, stop, pre, post. Extractall will also extracts multiple match of the pattern, i.e. if the peptid is mapping to multiples locations of the protein, this will return all of that mapping locations.
	start_stop_pre_post= pepObjUnique['proteinacc_start_stop_pre_post_;'].str.extractall(re.escape(prt)+"_"+"(?P<start>\d+)_(?P<stop>\d+)_(?P<pre>[A-Z]|-)_(?P<post>[A-Z]|-|\?)(?:;|$)")
	start_stop_pre_post[['PepQValue','Modifications','Charge','PSM Count']] = pepObjUnique.reindex(start_stop_pre_post.index,level=0)[['PSM-level q-value','Modifications','Charge','Spectrum ID']]
	#The dataframe is flattened to list of lists for easy access of all the peptide locations.
	start_stop_pepqval_mod_charge_psmcount=start_stop_pre_post[['start','stop','PepQValue','Modifications','Charge','PSM Count']].values.tolist()
	return start_stop_pepqval_mod_charge_psmcount
    
#'asmbl_28157|m.401505_294_315_R_R;asmbl_28158|m.401599_291_312_R_R;asmbl_28157|m.401505_291_312_R_R'
def findProteinAnnotation(prt, gene, gffObjCount):
    ##this function take protein accession [mRNA id if protein identification was MSGF+ and PASA was used] and find annotation from the gffObj.
    prtGffAnn=gffObjCount[gffObjCount['attributes'].str.contains(re.escape("ID="+gene+";")+"|"+re.escape(prt)+"(?:;|\.)")]
    if prtGffAnn.shape[0]>0: 
        prtGffCDS=prtGffAnn[prtGffAnn['type']=='CDS']
        aggreNCount=cdsAggregate(prtGffCDS)
        prtGffCDSAggre=pd.concat([prtGffCDS,aggreNCount],axis=1)
        print("CDS Aggre :"+str(prtGffCDSAggre))
        return prtGffCDSAggre
    else:
        print("WARNING:"+prt+" does not have a genomic location")
        return prtGffAnn

def main(prtObj,pepObj,gffObj,outFile):
    #this function takes three input, prtObj=MSGF+ protein wise result, pepObj=MSGF+ peptide results and gffObj=gff3 file produced by Transdecoder and filtered by PIT-DB_dataprocessing code.
    #Adds nucletide counts to this dataframe
    gffObjCount=nucleotideCount(gffObj)
    pIdx=prtObj.index.tolist()
    pepGff=pd.DataFrame(columns=gffObj.columns[:-2])
    for i in range(0,len(pIdx)):
        prt=prtObj.loc[pIdx[i]]['protein accession']
        gene=prtObj.loc[pIdx[i]]['description'].split(' ')[1]
        print("Protein:"+prt)
        prtGffCDSAggre=findProteinAnnotation(prt, gene, gffObjCount)
        if prtGffCDSAggre.shape[0]>0:
            #start_stop=findPeptidesStartStop(prt, pepObj)
            start_stop_pepqval_mod_charge_psmcount=findPeptidesStartStop(prt, pepObj)
            #pepGff=pepGff.append(peptidePosition(prtGffCDSAggre, prt, start_stop))
            pepGff=pepGff.append(peptidePosition(prtGffCDSAggre, prt, start_stop_pepqval_mod_charge_psmcount))
        else:
            print("WARNING:"+prt+" no genomic location or no CDS")
    print(pepGff.columns)
    pepGff.to_csv(outFile,sep="\t",columns=gffObj.columns.tolist(),index=False)

parser = argparse.ArgumentParser(description='Representing identified peptides in GFF3 file format')
parser.add_argument("--psm", required=True, help="mzIdentML-lib PSM export file name")
parser.add_argument("--gff3", required=True, help="Annotation GFF3 file name")
parser.add_argument("--proteins", required=True, help="mzIdentML-lib protein export file")

args = parser.parse_args()

PSMFile=args.psm
prtFile=args.proteins
gffFile=args.gff3
outFile=os.path.splitext(gffFile)[0]+"_peptide.gff3"
print ("PSM file:"+PSMFile)
print ("Protein file:"+prtFile)
print ("GFF3 file:"+gffFile)

pepObj=readFile(PSMFile,',',0)
prtObj=readFile(prtFile,',',0)
gffObj=readFile(gffFile,'\t',0)

main(prtObj,pepObj,gffObj,outFile)

#psmFile='/data/SBCS-BessantLab/shyama/Data/Oliver/Identification/G102/G102+fdr+th+grouping_filtered.csv'
#psm=readFile(psmFile,',',0)


