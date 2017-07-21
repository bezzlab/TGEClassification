
args <- commandArgs(trailingOnly = TRUE)
args.length <- length(args)
if(args.length < 4){
	print ("Error: Not sufficient parameters")
    print ("Rscript SampleSpecificScoringandModelling.R coseq.train.tsv peptides.csv conseq.test.tsv conseq.test.sss.tsv")
}
conseq.train.file=args[1]
peptide.file=args[2]
conseq.test.file=args[3]
conseq.test.outfile=args[4]
#conseq.file="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/AminoAcids-or-ORFs-orTGEs/human_adeno.assemblies.fasta.transdecoder.pep.identified.fasta_conseq.tsv"
#peptide.file="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/PSMs-Peptides-ORFs/human_adeno+fdr+th+grouping_filtered.csv"
conseq.matrix.train=read.table(conseq.train.file, header=TRUE, sep='\t')
peptides=read.table(peptide.file, header=TRUE, sep=',')
peptide.train.conseq=conseq.matrix.train[which(conseq.matrix.train[,'Peptide'] %in% peptides[,'Sequence']),]
conseq.matrix.test=read.table(conseq.test.file, header=TRUE, sep='\t')
##Train model.

data.points=seq(0.01,max(conseq.matrix.train$detectability.predicted.by.Random.Forest),0.01)
xy=matrix(0,nrow=length(data.points),ncol=4)
rownames(xy)=data.points
colnames(xy)=c('c2','sss','psm','specEval')
xy[,'c2']=data.points

sss=apply(xy,1,function(xy,x,y){
  windowStart=xy['c2']-0.05
  windowStop=xy['c2']+0.05
  ident=length(unique(y[which(windowStart<=y$detectability.predicted.by.Random.Forest & y$detectability.predicted.by.Random.Forest<=windowStop),'Peptide']))
  total=length(unique(x[which(windowStart<=x$detectability.predicted.by.Random.Forest & x$detectability.predicted.by.Random.Forest<=windowStop),'Peptide']))
  ident/total
}, conseq.matrix.train, peptide.train.conseq)

xy[,'sss']=sss
x=xy[,'c2']
y=xy[,'sss']
model <- lm(y ~ poly(x,3))

predicted.intervals.test=predict(model,data.frame(x=conseq.matrix.test[,'detectability.predicted.by.Random.Forest']),interval='confidence',level=0.99)
conseq.matrix.test[,'sss']=predicted.intervals.test[,1]
write.table(conseq.matrix.test, file=conseq.test.outfile, sep="\t", row.names=FALSE)