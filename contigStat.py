#!/usr/bin/python

from __future__ import division
from os.path import expanduser
import sys
import Bio
import os
import glob
import math

from Bio.Blast import NCBIXML

print("In:"+sys.argv[1])
#files=glob.glob(sys.argv[1]+"*.xml");
files=glob.glob(sys.argv[1]);
print(files)
home = expanduser("~");

path=sys.argv[2]
print("out:"+path)

if not os.path.exists(path):
    os.makedirs(path)

for infile in files: #glob.glob('blast/*.out'):
    if os.path.getsize(infile)>0:
        file_handle = open(infile)
        
        out_handle = open(os.path.join(path,os.path.basename(infile)+'.csv'),'w')
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


