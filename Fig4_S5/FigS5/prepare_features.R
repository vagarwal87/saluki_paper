library(imputeR)

#half lives, simple seq features, & codon freqs
a=readme("seqFeatWithKmerFreqs.txt.gz", gzip=T)
a=a[,c(1,grep("Codon",colnames(a)),grep("3UTR",colnames(a)))]
dim(a)

b=readme("CWCS.txt.gz",gzip=T)
colnames(b)[2:ncol(b)]=paste0("MIR.",colnames(b)[2:ncol(b)])
b[2:ncol(b)] = -1*b[2:ncol(b)]
a=merge(a,b,by=1,all.x=1)
a[is.na(a)]=0
say("with miRNA features")
dim(a)

seqpreds <- function(tool){
  humanorf=readme(paste0(tool,"_predictions_v83/ORF_avg.txt.gz"), gzip=T)
  human3p=readme(paste0(tool,"_predictions_v83/3pUTR_avg.txt.gz"), gzip=T)
  human5p=readme(paste0(tool,"_predictions_v83/5pUTR_avg.txt.gz"), gzip=T)

  colnames(humanorf)[2:ncol(humanorf)]=paste0(tool,colnames(humanorf)[2:ncol(humanorf)],".ORF")
  colnames(human3p)[2:ncol(human3p)]=paste0(tool,colnames(human3p)[2:ncol(human3p)],".3UTR")
  colnames(human5p)[2:ncol(human5p)]=paste0(tool,colnames(human5p)[2:ncol(human5p)],".5UTR")
  a=merge(humanorf,human3p,by=1,all=T)
  a=merge(a,human5p,by=1,all=T)
  a[is.na(a)]=0
  a
}

say("with SeqWeaver features")
b=merge(a,seqpreds("SeqWeaver"),by=1)
colnames(b)=paste0("GENETIC.",colnames(b))
dim(b)

load("alldata.RData")
#BEeM model
a=a[,c(1:10,grep("^ENCORI|^M6ACLIP|^eCLIP",colnames(a)))]
colnames(a)[3:10]=paste0("BASIC.",colnames(a)[3:10])
colnames(a)[11:ncol(a)]=paste0("BIOCHEMICAL.",colnames(a)[11:ncol(a)])

dim(a)
a=merge(a,b,by=1)
dim(a)

save(a,file="alldataGeneticBiochemical.RData")
