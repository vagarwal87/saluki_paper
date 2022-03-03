library(glmnet)

#mouse model tested w/ human data
load("../Robj/BC3MS_CV-Lasso.Robj")

#half lives, kmer kreqs, & codon freqs features
a=readme("seqFeatWithKmerFreqs.txt.gz",gzip=T)
a=a[,c(1:10,grep("Codon",colnames(a)),grep("3UTR",colnames(a)))]

#augment matrix w/ mouse miRs -- set to 0 since they don't exist in human
b=readme("../CWCS.txt.gz",gzip=T)
mirnames = paste0("MIR.",colnames(b)[2:ncol(b)])
df <- data.frame(matrix(ncol = length(mirnames), nrow = nrow(a)))
df[is.na(df)] = 0
colnames(df) <- mirnames
a=cbind(a, df)
say("with miRNA features")

seqpreds <- function(tool){
  mouseorf=readme(paste0(tool,"_predictions_mm10_90/ORF_avg.txt.gz"), gzip=T)
  mouse3p=readme(paste0(tool,"_predictions_mm10_90/3pUTR_avg.txt.gz"), gzip=T)
  mouse5p=readme(paste0(tool,"_predictions_mm10_90/5pUTR_avg.txt.gz"), gzip=T)

  colnames(mouseorf)[2:ncol(mouseorf)]=paste0(tool,colnames(mouseorf)[2:ncol(mouseorf)],".ORF")
  colnames(mouse3p)[2:ncol(mouse3p)]=paste0(tool,colnames(mouse3p)[2:ncol(mouse3p)],".3UTR")
  colnames(mouse5p)[2:ncol(mouse5p)]=paste0(tool,colnames(mouse5p)[2:ncol(mouse5p)],".5UTR")
  a=merge(mouseorf,mouse3p,by=1,all=T)
  a=merge(a,mouse5p,by=1,all=T)
  a[is.na(a)]=0
  a
}

say("with SeqWeaver features")
a=merge(a,seqpreds("SeqWeaver"),by=1)

say("With Ensembl IDs + Features: ", nrow(a))

say("Dimensions: ", dim(a))
genes = a[,1]
a[,1]=NULL

# x=scale(a[,2:ncol(a)])
x=scale(a[,rownames(cvfit$glmnet.fit$beta)]) #reorder columns by the coefficients stored in mouse model
x[is.na(x)]=0
y=scale(a[,1])

binnedgenes = read.delim("binnedgenes.txt")
yhat <- rep(NA,length(y))

for (i in 1:10) {
  predRows <- as.numeric(rownames(binnedgenes)[binnedgenes$BIN == i])
  yhat[predRows] = predict(cvfit,newx=x[predRows,],type="response",s=cvfit$lambda.min)
  say(cor(yhat[predRows], y[predRows]))
}
