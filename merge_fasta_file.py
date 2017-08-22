from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

print("sys")
print(sys.argv[1])
with open(sys.argv[len(sys.argv)-1],"w") as wrt:
	for i in range(1,len(sys.argv)-1):
		with open(sys.argv[i], "rU") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				wrt.write(">"+record.description+"\n")
				wrt.write(str(record.seq)+"\n")
				handle.close()
wrt.close()
