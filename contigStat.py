#!/usr/bin/python

from __future__ import division
from os.path import expanduser
import sys
#sys.path.append('/nfs/ma/home/shyama/installed_soft/BiopythonCluster/lib64/python')
import Bio
import os
import glob
import math

from Bio.Blast import NCBIXML

#path = '/Volumes/ma-home/SYBARIS/WGS-AF293/Shyama/test/xml_output/'
#path = '/Volumes/ma-home-1/shyama/outputs/SYBARIS/BLAST/OG/xml_output_default_gap/'
#path = '/Volumes/ma-home-1/shyama/outputs/SYBARIS/BLAST/AF293/xml_output_gap/'
#path = '/Volumes/ma-home/shyama/outputs/SYBARIS/BLAST/OG/xml_output_gap_coverage/'
#print sys.argv[1:];
#path = '/nfs/ma/home/shyama/DATA/SYBARIS/data/blast/' #sys.argv[1] #
#files = glob.glob('contigsAllReferenceGuidedBlast/*.csv')
#files = glob.glob('whole*/Trinity_blast_transcript.xml')
#files = glob.glob('E:\Data\HUMAN\database\Trinity\*.xml');
#files = glob.glob(r"D:\New folder\blast\uni*.xml");
#infile=path
print("In:"+sys.argv[1])
#files=glob.glob(sys.argv[1]+"*.xml");
files=glob.glob(sys.argv[1]);
print(files)
home = expanduser("~");

path=sys.argv[2]
print("out:"+path)

if not os.path.exists(path):
    os.makedirs(path)

#files=["blast/B07BNABXX_2_4.out","blast/C023MABXX_3_8.out","blast/D0ACKACXX_1_12.out","blast/B07BNABXX_1_9.out","blast/C023MABXX_5_2.out"]

for infile in files: #glob.glob('blast/*.out'):
    if os.path.getsize(infile)>0:
        file_handle = open(infile)
        
        out_handle = open(os.path.join(path,os.path.basename(infile)+'2.csv'),'w')
        #out_handle.write("file,query_name,query_length,match,hit_def,hit_length,e-value,p-value,identitie,align_length,gap,good_match,long_match,score,more,s_st,s_end\n") #seq,sseq,
        out_handle.write("query_name,query_length,match,hit_count,hit_def,hit_length,e-value,p-value,identitie,align_length,gap,good_match,long_match,score,more,q_st,q_end,s_st,s_end,seq,sseq,match_seq\n")
        blast_record_itr = Bio.Blast.NCBIXML.parse(file_handle)
        for blast_record in blast_record_itr:
            #out_handle.write(os.path.basename(infile) + ",")
            out_handle.write(blast_record.query.replace(",",";") + ",")
            out_handle.write(str(blast_record.query_letters)+",")
            if(len(blast_record.alignments)>0):
                out_handle.write("yes,")
                hit_count=len(blast_record.alignments)
                out_handle.write(str(hit_count))
                out_handle.write(",")
                for alignment in blast_record.alignments:
                    count=len(alignment.hsps)
                    out_handle.write(alignment.hit_def.replace(",",";"))
                    out_handle.write(",")
                    out_handle.write(str(alignment.length))
                    out_handle.write(",")
                    for hsp in alignment.hsps:
                        out_handle.write(str(hsp.expect))
                        out_handle.write(",")
                        out_handle.write(str((1-math.exp(-(hsp.expect)))))
                        out_handle.write(",")
                        out_handle.write(str(hsp.identities))
                        out_handle.write(",")
                        out_handle.write(str(hsp.align_length))
                        out_handle.write(",")
                        out_handle.write(str(hsp.gaps))
                        out_handle.write(",")
                        out_handle.write(str(round(hsp.identities/hsp.align_length,4)))
                        out_handle.write(",")
                        out_handle.write(str(round(hsp.identities/alignment.length,4)))
                        out_handle.write(",")
                        out_handle.write(str(hsp.score))
                        out_handle.write(",")
                        out_handle.write(str(count))
                        out_handle.write(",")
                        out_handle.write(str(hsp.query_start))
                        out_handle.write(",")
                        out_handle.write(str(hsp.query_end))
                        out_handle.write(",")
                        out_handle.write(str(hsp.sbjct_start))
                        out_handle.write(",")
                        out_handle.write(str(hsp.sbjct_end))
                        out_handle.write(",")
                        out_handle.write(hsp.query)
                        out_handle.write(",")
                        out_handle.write(hsp.sbjct)
                        out_handle.write(",")
                        out_handle.write(hsp.match+"\n")
                        break
                    break
            else:
                out_handle.write("no, , , , , , , , , , , , , , , , , , \n")
    else:
        print(infile+" is empty")


