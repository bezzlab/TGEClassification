## Main AA script 
# Main script to run Consequence 2
#
# 1) Load in necessary files, models etc 
# 2) Input peptide sequence and extract AA index 
# 3) Run predictor 
#
#get the current path where the script locates
#therefore the filenames are relative to this folder
curr.dir = getwd()

args <- commandArgs(trailingOnly = TRUE)
args.length <- length(args)
if(args.length < 4){
	print ("Error: Not sufficient parameters")
	usage()
}

in.peptide.file <- args[1]
out.peptide.file <- args[2]
work.dir <- args[3]

models <- args[4:args.length]
for (predictor.interest in models){
	print (predictor.interest)
	#check the selected predictor is defined
	if(! (predictor.interest == "rf" | predictor.interest == "nn" | predictor.interest == "ada")){
		error.msg <- paste("Error: Wrong value (",predictor.interest,sep="")
		error.msg <- paste(error.msg,") for the model",sep="")
		print (error.msg)
		usage();
	}
}

## load libraries 
library(ada)
library(nnet)
library(randomForest)
library(seqinr)
#library(caret)

# Directories 
map.dir <- paste(work.dir,"/maps",sep="")
code.dir <- paste(work.dir,"/code",sep="")
predictor.dir <- paste(work.dir,"/predictors",sep="")

#load the library which has fixed location
setwd(code.dir)
source("main_functions.r")

#read in the input file which is relative to this script
setwd(curr.dir)
peptide.matrix <- table.reader(in.peptide.file)

#load the feature map which has fixed location
setwd(map.dir)
AA.interest <- readLines("feature_map.txt")

peptide.aa <- as.data.frame(returnAA(peptide.matrix[, "Peptide"], AA.interest))
colnames(peptide.aa) <- paste("f", 1:length(AA.interest), sep = "")

for (predictor.interest in models){
	setwd(predictor.dir)
	if(predictor.interest == "rf"){
		load("rf_predictor.RData")
		predictor.title <- "Random Forest"
		pred.run <- rfFit
	}
	if(predictor.interest == "nn"){
		load("nn_predictor.RData")
		predictor.title <- "Neural Network"
		pred.run <- nnFit
	}
	if(predictor.interest == "ada"){
		load("ada_predictor.RData")
		predictor.title <- "Adaptive booster"
		pred.run <- adaFit
	}

	# produce the AA indices 
	prob.matrix <- round(predict(pred.run, peptide.aa, type = "prob"), 2)

	detectability.matrix <- as.matrix(prob.matrix[,1])
	colnames(detectability.matrix) <- c(paste("detectability predicted by",predictor.title,sep=" "))
	peptide.matrix <- cbind(peptide.matrix, detectability.matrix)
}
#write to the output file which is relative to this script
setwd(curr.dir)
write.table(peptide.matrix, out.peptide.file, row.names = FALSE, quote = FALSE, sep = "\t")


















