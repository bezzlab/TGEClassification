## main functions 
usage <- function(){
print ("Usage: Rscript consequence.r <input.file> <output.file> <model>")
print ("Model should be one of nn(neural network), rf(random forest) and ada(Adaptive booster)")
q()	
}

AAIndex.Matrix.Generator <- function(){
library(seqinr)
data(aaindex)
aaindex.Matrix <- matrix(nrow = length(aaindex), ncol = 20)
colnames(aaindex.Matrix) <- c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val")
rownames(aaindex.Matrix) <- c(names(aaindex))

for(i in 1:length(aaindex)){
aaindex.Matrix[i, ] <- aaindex[[i]][["I"]]
}
rm(i)

# Change from triple to single code
colnames(aaindex.Matrix) <- a(colnames(aaindex.Matrix))
return(aaindex.Matrix)
}

returnAA <- function(In.Peptides, In.Optimal){
Positive.Set <- lapply(as.list(In.Peptides), function(x){
Initial.Split <- strsplit(x, split = "", fixed = TRUE)[[1]];
Initial.Split <- Initial.Split[which(Initial.Split != "X")];
Out.Split <- Initial.Split[which(Initial.Split != "U")]
})

AA.Matrix <- AAIndex.Matrix.Generator()
AA.Optimal <- AA.Matrix[In.Optimal, ]

positive.matrix <- matrix(nrow = length(In.Peptides), ncol = length(In.Optimal))
rownames(positive.matrix) <- In.Peptides
colnames(positive.matrix) <- In.Optimal
for(i in 1:length(Positive.Set)){
positive.matrix[i, ] <- t(as.matrix(round(apply(AA.Optimal[, Positive.Set[[i]]], 1, function(y){
mean(y[which(is.na(y) == FALSE)])
}), 2)))
}
rm(i) 
return(positive.matrix)
}


table.reader <- function(in.file){
# in.file <- "peptides_variability.tsv"
file.in <- readLines(in.file)
file.headers <- strsplit(file.in[1], split = "\t", fixed = TRUE)[[1]]
max.columns <- length(file.headers)

tmp.matrix <- t(apply(as.matrix(file.in[2:length(file.in)]), 1, function(x){
tmp.split <- strsplit(x, split = "\t", fixed = TRUE)[[1]]; 
if(length(tmp.split) != max.columns){
add.set <- max.columns - length(tmp.split) 
tmp.split <- c(tmp.split, rep("", add.set))
}
return(tmp.split)
}))

colnames(tmp.matrix) <- file.headers
return(tmp.matrix)
}


