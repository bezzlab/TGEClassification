###### this code reads the the output of contigStat.py which is essentially BLAST result in csv format. Second input to this an excel file downloaded
###### from the Uniprot website that tells us the chomosome location of the proteins. In future this step needs to be done programatically. Third and
###### the last parameter to this script is the output path. Output of this script is similar to that of the first input file, with an extra column for
###### the location of the protein.

import os
import csv
import re
import argparse

#findProtein('P62258',protList)
#Consider chromosome X and Y.
def findProtein(proteinId, proteinList):
    chromosome=[element[1] for element in proteinList if element[0]==proteinId]## should always have one element in this list. unless this proteinId does not exist in the proteinList, which is possible.
    if len(chromosome)==0:
        return 0;
    else:
        if len(chromosome)==1:
            #this is what expected.
            ##split the element and extract only the chromosome number, e.g. outputs ['UP000005640: Chromosome 17, UP000005640: Unplaced'] or ['UP000005640: Chromosome 17']
            print(chromosome)
            chrX=re.split(',',chromosome[0]);
            #now remove 'UP000005640: Chromosome' from it.
            chrStr="";
            for ch in chrX:
                #print(ch)
                c=re.search('(\d+$)|(X$)|(Y$)|(M$)|(Mitochondrion$)',ch)
                #total number of pattern in above statement decides length of c.groups.
                #print(ch[0])
                if c!=None:
                    print(c.groups())
                    cStr=[ele for ele in c.groups()[0:] if ele!=None]
                    if len(cStr)==1:
                        cStr[0]=re.sub("Mitochondrion","MT",cStr[0])
                        if chrStr=="":
                            chrStr=cStr[0];
                        else:
                            chrStr=chrStr+","+cStr[0]
                    else:
                        print("\n\nERROR: more than chromosome given:"+ch)
                else:
                    print("Inadequate pattern to capture the location info or this protein is unplaced:"+ch)
            if chrStr=="":
                print(proteinId+": chromosome not found");
                chrStr="N/A"
            print(chrStr)
            return chrStr;
        else:
            return -1;
    

def readUniFile(filename,nameCol,locationCol):
    proteins=[]
    with open(filename, newline='') as protFile:
        next(protFile)
        reader = csv.reader(protFile, delimiter='\t');
        for r in reader:
            proteins.append([r[nameCol],r[locationCol]])
    protFile.close()
    return proteins;

#read('D:/data/blast/blastCSV/human_adeno_mydb_pasa.assemblies_ORFs.csv','D:/data/Data/uniprot-Homo+sapiens960613_04_2015.xls','D:/data/blast/blastCSV/human_adeno_mydb_pasa.assemblies_ORFs_with_Location.csv')
def read(blastFileName, uniProteinFileName, newBlastFileName, index, nameCol, locationCol):
    subjectProtCol=4;
    with open(blastFileName, 'r', newline='') as blastFile, open(newBlastFileName,'w',newline='') as newBlastFile:
        reader = csv.reader(blastFile, delimiter=',') #,quoting=csv.QUOTE_NONNUMERIC
        fileWriter = csv.writer(newBlastFile, delimiter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        count=0
        protList=readUniFile(uniProteinFileName, nameCol, locationCol)
        #print(len(protList))
        for line in reader:
            
            if count>0:
                pId=re.search('(?:\|)(.+?)(?:\|)',line[index]) # here group(0) will return text including '|'s, and group(1) returns only (.+?). this regex is very uniprot specific. in future, if we get reference proteome from non uniprot source, this should be revisited.
                if pId!=None:
                    #print(line)
                    #print(pId.groups())
                    pIdCan=pId.group(1).split('-')
                    print(pIdCan)
                    if len(pIdCan)>0:
                        line.append(findProtein(pIdCan[0],protList))
                    else:
                        print("This should not happen about the protein Id.")
                        sys.exit();
                else:
                    line.append("")
                if line[0]=="Dataset_A_asmbl_43340_ORF4_Frame_2_26-313":
                    print(line);
                fileWriter.writerow(line)
            else:
                line.append("Location")
                fileWriter.writerow(line)
                count=count+1
        #fileWriter.close()

parser = argparse.ArgumentParser(description='This script read output of contigStat.py and uniprot protein location file and merge the information')
parser.add_argument("-b", "--blast", nargs=1, required=True, help="full path of blast csv file", metavar="PATH")
parser.add_argument("-u", "--uniprot", nargs=1, required=True, help="full path of uniprot file", metavar="PATH")
parser.add_argument("-o", "--output", nargs=1, required=True, help="full path of output file", metavar="PATH")
parser.add_argument("-n", "--name", nargs=1, default=[0], help="full path of output file", metavar="PATH")
parser.add_argument("-l", "--location", nargs=1, default=[8], help="full path of output file", metavar="PATH")
args = parser.parse_args()
print(args)
idx=4
read(args.blast[0], args.uniprot[0], args.output[0],idx, args.name[0], args.location[0])
#read('D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs.csv','D:/data/Data/uniprot-Homo+sapiens960613_04_2015.xls','D:/data/blast/blastCSV/human_adeno_mydb_pasa.assemblies_ORFs_with_LocationV3.csv')
#D:\Code\Proteomics\Python\IsoformsSAP>python UniProteinLocation.py -b D:\data\Bristol\HumanAdeno\PASATransdecoder\pasa_transdecoder_nonstar_identified.csv -u D:/data/Data/uniprot-Homo+sapiens960613_04_2015.xls -o D:\data\Bristol\HumanAdeno\PASATransdecoder\pasa_transdecoder_nonstar_identified_Location.csv
