#!/bin/sh
#$ -cwd              # Set the working directory for the job to the current directory
#$ -V
#$ -l h_rt=120:0:0    # Request 240 hour runtime
#$ -l h_vmem=16G      # Request 4GB RAM

###This code create BLAST database from Uniprot proteome fasta file and BLAST the identified ORFs (passed as argument) agaist it.

usage()
{
	echo "usage: $0 [-u proteomefasta] [-p orffasta] [-o outfolder] [-m BLOSUM_MATRIX]"
}


while getopts u:p:o:m: opt
do
    case $opt in
      (u)  uni="$OPTARG";;
      (p)  prot="$OPTARG";;
      (o)  out="$OPTARG";;
      (m)  matrix="$OPTARG";;
      (\?)
      	  usage
	  exit;;
    esac
done
shift `expr $OPTIND - 1`

if [ -z "$uni" ]; then
	usage
	exit
fi

if [ -z "$prot" ]; then
        usage
        exit
fi

if [ -z "$out" ]; then
        usage
        exit
fi

if [ -z "$matrix" ]; then
        matrix="BLOSUM80"
fi

#echo $uni
#echo $prot
#echo $out

proteome=$(basename $uni | sed -e "s/.fasta//g") 
echo $proteome

fastaPath=$(dirname $uni)
echo $fastaPath

if [ -f $fastaPath/$proteome.aa.pin ]; then
	echo "BLAST DB exist"
else
	echo "ERROR: The BLASTDb does not exist"
	echo "makeblastdb -in $uni -dbtype prot -out $fastaPath/$proteome.aa"
	#makeblastdb -in $uni -dbtype prot -out $fastaPath/$proteome.aa
fi

if [ -f $prot ]; then
    sample=$(basename $prot | sed -e "s/.assemblies.fasta.transdecoder.pep.identified.fasta//g")
	#cd $DIR$sample
	echo $sample
	cp $prot $prot.nostar.fasta
	mv $prot.nostar.fasta $out
	sed -i 's/\*//g' $out$sample.assemblies.fasta.transdecoder.pep.identified.fasta.nostar.fasta
	if [ ! -d $out ]; then
		mkdir -p $out
	fi
	echo "blastp -db $fastaPath/$proteome.aa -query $prot.nostar.fasta -gapopen 6 -gapextend 2 -out $out/$sample.assemblies.fasta.transdecoder.pep.identified.canonical.xml -outfmt 5 -evalue 1 -matrix BLOSUM80"
	#blastp -db $fastaPath/$proteome.aa -query $prot.nostar.fasta -gapopen 6 -gapextend 2 -out $out/$sample.assemblies.fasta.transdecoder.pep.identified.canonical.xml -outfmt 5 -evalue 1 -matrix BLOSUM80
	blastp -db $fastaPath/$proteome.aa -query $out$sample.assemblies.fasta.transdecoder.pep.identified.fasta.nostar.fasta -gapopen 6 -gapextend 2 -out $out/$sample.assemblies.fasta.transdecoder.pep.identified.xml -outfmt 5 -evalue 1 -matrix $matrix
	cd /data/home/btw796/Code2/SYBARIS/Python
	#python3.4 contigStat.py $out/$sample.assemblies.fasta.transdecoder.pep.identified.canonical.xml $out
	#python contigStat.py $out/$sample.assemblies.fasta.transdecoder.pep.identified.xml $out
fi



