#!/bin/sh
###This code run all the components of the TGE Classification.
usage()
{
	echo "usage: $0 [-u proteomefasta] [-p orffasta] [-r transcriptfasta] [-o outfolder] [-s mgffile] [-m modificationsFile] [-t tolerance deafault 10ppm] [-i instrument, deafault 1, options 1. Orbitrap/FTICR, 2. TOF, 3. Q-Exactive] [-f fragment method, default 1. Possoble values 1.CID, 2. ETD, 3. HCD] [-d decoysearch 0|1 default 1] [-c contaminantfile default crap.fasta] [-l minlengthpeptide interger default 8] [-v TSV file conatining reference protein location] [-g transdecoder generated sample.genome.gff3 file]"
}

while getopts u:p:r:o:s:m:t:i:f:d:c:l:v:g: opt
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

if [ ! -d $out/blast ];then
  mkdir $out/blast
fi

#PIT
python runProteinIdentificationAndPostProcessing_cluster.py -i $ms -o $out/PSMs-Peptides-ORFs/$sample.mzid -d $prot -m $mod -t $tolerance --m_msgf $frag \
  --inst $inst --minLength $peplen --tda $decoy --contaminants $cont > $out/run_$sample.cmds.sh
sh $out/run_$sample.cmds.sh #| qsub -cwd -l h_vmem=40G -l h_rt=140:0:0

#Standard
python runProteinIdentificationAndPostProcessing_cluster.py -i $ms -o $out/PSMs-Peptides-ORFs/$sample.standard.mzid -d $uni -m $mod -t $tolerance --m_msgf $frag \
  --inst $inst --minLength $peplen --tda $decoy --contaminants $cont > $out/run_$sample.standard.cmds.sh
sh $out/run_$sample.standard.cmds.sh


#Post-processing
python $PIT_DBData_processing --ORFs $prot --transcripts $transcript --proteins $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping+prt.csv --peptides $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping.csv --ORFsOutFile $out/AminoAcids-or-ORFs-orTGEs/$sample.identified.fasta --transcriptsOutFile $out/transcripts/$sample.identified.fasta --proteinOutFile $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping+prt_filtered.csv --peptideOutFile $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping_filtered.csv
python standardSearchResultProcessing.py --protein $out/PSMs-Peptides-ORFs/$sample.standard+fdr+th+grouping+prt.csv --protOut $out/PSMs-Peptides-ORFs/$sample.standard+fdr+th+grouping+prt_filtered.csv --peptide $out/PSMs-Peptides-ORFs/$sample.standard+fdr+th+grouping.csv --pepOut $out/PSMs-Peptides-ORFs/$sample.standard+fdr+th+grouping_filtered.csv --fasta $uni --fastaOut $out/$sample.standard.identified.fasta
sh runBlastBristol.sh -u $uni -p $out/AminoAcids-or-ORFs-orTGEs/$sample.identified.fasta -o $out/blast/
python contigStat.py $out/blast/$sample.identified.xml $out/blast/
python $uniloc -b $out/blast/$sample.identified.xml.csv -u $uniLoc -o $out/blast/$sample.identified.loc.csv 	
python $isoSAP -b $out/blast/$sample.identified.loc.csv -k $out/blast/$sample"_known.csv" -s $out/blast/$sample"_knownVar.csv" -v $out/blast/$sample.vcf -i $out/blast/$sample"_iso.csv" -j $out/blast/$sample"_isoVar.csv" > $out/$sample.log
python $preliAnnot -b $out/blast/$sample.identified.loc.csv -o $out/Annotation/ -k $out/blast/$sample"_known.csv" -s $out/blast/$sample"_knownVar.csv" -i $out/blast/$sample"_iso.csv" -j $out/blast/$sample"_isoVar.csv" -v $out/blast/$sample.vcf
python SplitAnnotationFile.py -a $out/Annotation/$sample"_annotation.csv" -o $out/Annotation/$sample"_details_annotation.csv"
if [ -z "$gff" ]; then
	echo "No GFF3 was paased, hence TGEs will be unavailable for visualization."
else
	python $annotCode -p $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping+prt_filtered.csv -g $transcript.transdecoder.genome.gff3 -o $out/GFF3/$sample"_identified.gff3"
fi
python $pepEvd -p $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping_filtered.csv -s $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping_filtered_variation.csv -v $out/blast/$sample.vcf -o $out/Variations-proVCF/$sample"_pepEvd.vcf"
python peptideEvidenceIsoforms.py -b $out/blast/$sample.identified.loc.csv -a $out/Annotation/$sample"_details_annotation.csv" -p $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping_filtered.csv -o $out/Variations-proVCF/$sample"_isoform.vcf"
python IsoformScoring.py --annot $out/Annotation/$sample"_details_annotation.csv" \
	--vcf $out/Variations-proVCF/$sample"_pepEvd.vcf" \
	--vcfiso $out/Variations-proVCF/$sample"_isoform_pepEvd.vcf" \
	--pep $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping_filtered.csv \
	--rpep $out/PSMs-Peptides-ORFs/$sample".standard+fdr+th+grouping_filtered.csv" \
	--digesttrain "$out/Score/$sample.train.digest.tsv" \
	--digesttest "$out/Score/$sample.test.digest.tsv" \
	--conseqtrain "$out/Score/$sample.train.conseq.tsv" \
	--conseqtest "$out/Score/$sample.test.conseq.tsv" \
	--ssstest "$out/Score/$sample.test.sss.tsv" \
	--orfs $out/AminoAcids-or-ORFs-orTGEs/$sample".identified.fasta" \
	--refs $out/AminoAcids-or-ORFs-orTGEs/$sample".standard.identified.fasta" \
	--score $out/Score/$sample".scoreFunction4.tsv"

python integratePeptideEvidenceInGFF3.py --psm $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping_filtered.csv --gff3 $out/GFF3/$sample"_identified.gff3" --proteins $out/PSMs-Peptides-ORFs/$sample+fdr+th+grouping+prt_filtered.csv
