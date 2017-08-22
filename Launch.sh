#!/bin/sh
###This code run all the components of the TGE Classification.
usage()
{
	echo "usage: $0 [-u proteomefasta] [-p orffasta] [-r transcriptfasta] [-o outfolder] [-s mgffile] [-m modificationsFile] [-t tolerance deafault 10ppm] [-i instrument, deafault 1, options 1. Orbitrap/FTICR, 2. TOF, 3. Q-Exactive] [-f fragment method, default 1. Possoble values 1.CID, 2. ETD, 3. HCD] [-d decoysearch 0|1 default 1] [-c contaminantfile default crap.fasta] [-l minlengthpeptide interger default 8] [-v TSV file conatining reference protein location] [-g transcript gff3 file]"
}

while getopts u:p:r:o:s:m:t:i:f:d:c:l:v:g opt
do
	case $opt in
	  (u)  uni="$OPTARG";;
	  (p)  prot="$OPTARG";;
	  (r)  transcript="$OPTARG";;
	  (o)  out="$OPTARG";;
	  (s)  ms="$OPTARG";;
	  (m)  mod="$OPTARG";;
	  (t)  tolerance="$OPTARG";;
	  (i)  inst="$OPTARG";;
	  (f)  frag="$OPTARG";;
	  (d)  decoy="$OPTARG";;
	  (c)  cont="$OPTARG";;
	  (l)  peplen="$OPTARG";;
	  (v)  uniLoc="$OPTARG";;
	  (g)  gff="$OPTARG";;
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

if [ -z "$transcript" ]; then
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

if [ -z "$ms" ]; then
	usage
	exit
fi

if [ -z "$mod" ]; then
		usage
		exit
fi

if [ -z "$inst" ]; then
	inst='1'
fi

if [ -z "$frag" ]; then
	frag='1'
fi

if [ -z "$decoy" ]; then
	decoy='1'
fi

if [ -z "$tolerance" ]; then
	tolerance='10ppm'
fi

if [ -z "$peplen" ]; then
	peplen='8'
fi

if [ -z "$decoy" ]; then
	decoy='1'
fi

if [ -z "$cont" ]; then
	echo "cont"
	cont="./data/crap.fasta"
	echo $cont
fi

PIT_DBData_processing=PIT-DBData_processing.py
uniloc=UniProteinLocation.py
isoSAP=IdentifyProteinIsoformSAP.py
pepEvd=peptideEvidence.py
annotCode=annotationMatrix.py
preliAnnot=PreliminaryProteinAnnotationForPITDBV2.py
sample=$(basename $prot | sed -e "s/.assemblies.fasta.transdecoder.pep//g")

if [ ! -d "$out/AminoAcids-or-ORFs-orTGEs/" ]; then
	echo "AminoAcids-or-ORFs-orTGEs directory should be created"
  	mkdir $out/AminoAcids-or-ORFs-orTGEs
else
	echo "$out/AminoAcids-or-ORFs-orTGEs exist"
fi

if [ ! -d $out/transcripts ]; then
  mkdir $out/transcripts
fi

if [ ! -d $out/PSMs-Peptides-ORFs ]; then
  mkdir $out/PSMs-Peptides-ORFs
fi

if [ ! -d $out/Annotation ]; then
  mkdir $out/Annotation
fi

if [ ! -d $out/Score ]; then
  mkdir $out/Score
fi

if [ ! -d $out/Summary ];then
  mkdir $out/Summary
fi

if [ ! -d $out/Variations-proVCF ];then
  mkdir $out/Variations-proVCF
fi

if [ ! -d $out/GFF3 ];then
  mkdir $out/GFF3
fi

#PIT
python runProteinIdentificationAndPostProcessing_cluster.py -i $ms -o $out/PSMs-Peptides-ORFs/$sample.mzid -d $prot -m $mod -t $tolerance --m_msgf $frag \
  --inst $inst --minLength $peplen --tda $decoy --contaminants $cont > $out/run_$sample.cmds.sh
sh $out/run_$sample.cmds.sh

#Standard
python runProteinIdentificationAndPostProcessing_cluster.py -i $ms -o $out/PSMs-Peptides-ORFs/$sample.standard.mzid -d $uni -m $mod -t $tolerance --m_msgf $frag \
  --inst $inst --minLength $peplen --tda $decoy --contaminants $cont > $out/run_$sample.standard.cmds.sh
sh $out/run_$sample.standard.cmds.sh

# #Post-processing
python $PIT_DBData_processing --ORFs $prot --transcripts $transcript --proteins $out/$sample+fdr+th+grouping+prt.csv --peptides $out/$sample+fdr+th+grouping.csv --ORFsOutFile $out/AminoAcids-or-ORFs-orTGEs/$sample.identified.fasta --transcriptsOutFile $out/transcripts/$sample.identified.fasta --proteinOutFile $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping+prt_filtered.csv --peptideOutFile $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping_filtered.csv
python standardSearchResultProcessing.py --protein $out/$sample.standard+fdr+th+grouping+prt.csv --protOut $out/$sample.standard+fdr+th+grouping+prt_filtered.csv --peptide $out/$sample.standard+fdr+th+grouping.csv --pepOut $out/$sample.standard+fdr+th+grouping_filtered.csv --fasta $uni --fastaOut $out/$sample.standard.identified.fasta
sh runBlastBristol.sh -u $uni -p $out/AminoAcids-or-ORFs-orTGEs/$sample.identified.fasta -o $out
python contigStat.py $out/$sample.assemblies.fasta.transdecoder.pep.identified.xml $out
python $uniloc -b $out/$sample.identified.xml.csv -u $uniLoc -o $out/$sample.identified.loc.csv	
python $isoSAP -b $out/$sample.identified.loc.csv -k $out/$sample_known.csv -s $out/$sample_knownVar.csv -v $out/$sample.vcf -i $out/$sample_iso.csv -j $out/$sample_isoVar.csv > $out/$sample.log
python $preliAnnot -b $out/$sample.identified.loc.csv -o $out/Annotation/
python SplitAnnotationFile.py -a $out/Annotation/ -o $out/Annotation/$sample.assemblies.fasta.transdecoder.pep_details_annotation.csv
# if [ -z "$gff" ]; then
# 	echo "No GFF3 was paased, hence TGEs will be unavailable for visualization."
# else
# 	python $annotCode -p $out/$sample+fdr+th+grouping+prt_filtered.csv -g $out/$sample.genome.gff3
# fi
# python $pepEvd -p $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping_filtered.csv -s $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping_filtered_variation.csv -v $out/$sample.vcf -o $out/$sample_pepEvd.vcf
# python peptideEvidenceIsoforms.py -b $out/$sample.identified.loc.csv -a -p -o
# python IsoformScoring.py --annot "$pitdb"Annotation/"$sample".assemblies.fasta.transdecoder.pep_details_annotation.csv \
# 	--vcf "$pitdb"Variations-proVCF/"$sample".assemblies.fasta.transdecoder.pep_pepEvd.vcf \
# 	--vcfiso "$pitdb"Variations-proVCF/"$sample".assemblies.fasta.transdecoder.pep_annotation.csv_isoform_pep.vcf \
# 	--pep "$pitdb"PSMs-Peptides-ORFs/"$sample"+fdr+th+grouping_filtered.csv \
# 	--rpep "$ident$sampleBase"/$sample.standard+fdr+th+grouping_filtered.csv \
# 	--digesttrain "$out$sample.train.digest.tsv" \
# 	--digesttest "$out$sample.test.digest.tsv" \
# 	--conseqtrain "$out$sample.train.conseq.tsv" \
# 	--conseqtest "$out$sample.test.conseq.tsv" \
# 	--ssstest "$out$sample.test.sss.tsv" \
# 	--orfs "$pitdb"AminoAcids-or-ORFs-orTGEs/"$sample".assemblies.fasta.transdecoder.pep.identified.fasta \
# 	--refs $ident$sampleBase/$sample.standard.identified.fasta \
# 	--score "$out$sample.scoreFunction4.tsv"