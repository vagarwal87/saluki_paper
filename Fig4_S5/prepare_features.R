library(imputeR)

#half lives, simple seq features, & codon freqs
a=readme("seqFeatWithKmerFreqs.txt.gz",gzip=T)
a=a[,1:10]
dim(a)

say("RIP/m6aCLIP")
b=readme("allRIPCLIPdata.txt.gz",gzip=T)
a=merge(a, b, by=1, all.x=T)
dim(a)

say("ENCORI")
load("peakMatENCORI.RData")
colnames(b)[2:ncol(b)]=paste0("ENCORI.",colnames(b)[2:ncol(b)])
a=merge(a, b, by=1, all.x=T)
dim(a)

say("PARCLIP")
b=readme("PARCLIP-master-table.txt")
colnames(b)[2:ncol(b)]=paste0("PARCLIP.",colnames(b)[2:ncol(b)])
b=aggregate(b[2:ncol(b)],by=list(b[,1]),max)
a=merge(a,b,by=1,all.x=T)
dim(a)

say("RIP")
b=readme("RIP-master-table.txt")
colnames(b)[2:ncol(b)]=paste0("RIP.",toupper(colnames(b)[2:ncol(b)]))
b=merge(a,b,by=1, all.x=T)
dim(b)

say("eCLIP")
load("eCLIP.RData")
a=merge(b, a, by=1, all.x=T)
dim(a)

apply(a,2,function(x) sum(is.na(x)))

#those without peaks have 0 peaks
clip=a[, grepl("CLIP|ENCORI",colnames(a))]
colnames(clip)
clip[is.na(clip)]=0
clip=log10(clip+1)

#impute continuous (i.e., non-peak) features
noclip=a[, !grepl("CLIP|ENCORI",colnames(a))]
colnames(noclip)
noclip[,3:ncol(noclip)]=impute(noclip[,3:ncol(noclip)],lmFun='plsR')$imp #impute continuous features

a=cbind(noclip, clip)

save(a,file="alldata.RData")
