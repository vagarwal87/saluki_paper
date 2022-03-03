library(biomaRt)
library(tidyr)
library(ggplot2)

ensembl = useMart(host='ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl')
annotations = getBM(attributes=c('ensembl_gene_id', 'external_gene_name','refseq_mrna', "chromosome_name"), mart = ensembl)
head(annotations,10)
validchr = c(1:100,'X','Y','MT')
annotations = annotations[annotations$chromosome_name %in% validchr,]

#gene symbols
a=read.delim("mouse/Sharova/dsn030-dsn030_Table2.csv")
a=a[,1:5]
a[,1]=gsub("\\*", "", a[,1])
a[a==0]=NA
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
colnames(a)=c("gene",paste0("Ko_ActD_mESC_",colnames(a)[2:ncol(a)]))
b=a

a=readme("mouse/Darnell/Ke.txt.gz",gzip=T)
a=a[,c(1,3)]
a[,2]=log10(a[,2]+0.1)
colnames(a)=c("gene","Darnell_ActD_mESC_1")
b=merge(b,a,by=1,all=T)

a=read.delim("mouse/Dolken/41586_2019_1369_MOESM3_ESM.csv")
a=a[!is.infinite(a[,2]),]
a[,2]=log10(a[,2]+0.1)
a=aggregate(a[,2],by=list(a[,1]),mean)
colnames(a)=c("gene","Dolken_4sU_M2-10B4_1")
b=merge(b,a,by=1,all=T)

a=read.delim("mouse/Hanna/degradation_rate.non_m6a.txt")
z=read.delim("mouse/Hanna/degradation_rate.m6a.txt")
a=rbind(a,z)
a=a[,c(2,3,5)]
a=a[a[,2]>0 & a[,3]>0, ]
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
colnames(a)=c("gene","Hanna_ActD_EB_1","Hanna_ActD_mESC_1")
b=merge(b,a,by=1,all=T)

a=read.delim("mouse/Herzog/Herzog_mESC_half_life.txt",skip=2)
a=a[,c(4,7)]
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
a=aggregate(a[,2],by=list(a[,1]),mean)
colnames(a)=c("gene","Ameres_4sU_mESC_1")
b=merge(b,a,by=1,all=T)

a=read.delim("mouse/Koszinowski/half_life_by_gene.txt")
a[,2]=log10(a[,2]+1)
a=aggregate(a[,2],by=list(a[,1]),mean)
colnames(a)=c("gene","Koszinowski_4sU_3T3_1")
b=merge(b,a,by=1,all=T)

a=read.delim("mouse/Rissland_DRUID/GSE99517_NIH3T3.4SU.HL.one.csv",sep=',')
a=a[,1:2]
a[,2]=log10(a[,2]+0.1)
colnames(a)=c("gene","Rissland2_4sU_3T3_1")
b=merge(b,a,by=1,all=T)

a=read.delim("mouse/Rissland_DRUID/GSE99517_NIH3T3.4SU.HL.two.csv",sep=',')
a=a[,1:2]
a[,2]=log10(a[,2]+0.1)
colnames(a)=c("gene","Rissland2_4sU_3T3_2")
b=merge(b,a,by=1,all=T)

a=read.delim("mouse/Schoefield/MEF_41592_2018_BFnmeth4582_MOESM5_ESM.csv")
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
a=a[,c(1,4,5)]
colnames(a)=c("gene","Simon_4sU_MEF_1","Simon_4sU_MEF_2")
b=merge(b,a,by=1,all=T)

a=read.delim("mouse/Wilusz/Dataset_S1.csv",skip=2)
a=a[,c(2,4,7,10)]
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
a=aggregate(a[,2:ncol(a)],by=list(a[,1]),mean)
colnames(a)=c("gene","Wilusz_ActD_C2C12_1","Wilusz_ActD_C2C12_2","Wilusz_ActD_C2C12_3")
b=merge(b,a,by=1,all=T)

a=read.delim("mouse/Zimmer/half_life_by_gene.txt")
a[,2]=log10(a[,2]+1)
a=aggregate(a[,2],by=list(a[,1]),mean)
colnames(a)=c("gene","Zimmer_ActD_3T3_1")
b=merge(b,a,by=1,all=T)

b=merge(b,annotations[,1:2],by.x=1,by.y=2,all=F)
b$gene=NULL
nrow(b)
head(b)
b=aggregate(b[,1:(ncol(b)-1)],by=list(b$ensembl_gene_id),function(x){ mean(x, na.rm=T) } )
colnames(b)[1]=c("gene")
nrow(b)

#Refseq
a=read.delim("mouse/Mattick/half_life_by_NMID2.txt",skip=1)
a[,2]=log10(a[,2]+0.1)
a=aggregate(a[,2],by=list(a[,1]),mean)
colnames(a)=c("gene","Mattick_ActD_Neuro2a_1")
c=a

a=read.delim("mouse/Eisen/mmc2.csv")
a=a[,c(1,3,4)]
colnames(a)=c("gene","Bartel2_5EU_3T3_1","Bartel2_5EU_3T3_2")
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
c=merge(c,a,by=1,all=T)

a=read.delim("mouse/Rabani/dg_rate_predictions_2011.txt")
a=a[,c(1,3)]
a[,2]=log10(a[,2]+1)
colnames(a)=c("gene","Regev_4sU_Dendritic_1")
c=merge(c,a,by=1,all=T)

a=read.delim("mouse/Rabani/rates_gene.dg.txt")
a=a[,c(1,3)]
a[,2]=log10(a[,2]+1)
colnames(a)=c("gene","Regev2_4sU_Dendritic_1")
c=merge(c,a,by=1,all=T)

a=read.delim("mouse/Schwanhausser/half_lives.txt",F)
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
colnames(a)=c("gene","Selbach_4sU_3T3_1","Selbach_4sU_3T3_2")
c=merge(c,a,by=1,all=T)

c=merge(c,annotations[,c(1,3)],by.x=1,by.y=2,all=F)
c$gene=NULL
nrow(c)
c=aggregate(c[,1:(ncol(c)-1)],by=list(c$ensembl_gene_id),function(x){ mean(x, na.rm=T) } )
colnames(c)[1]="gene"
nrow(c)


#Ensembl
a=read.delim("mouse/Spies/spies_half_lives.txt")
a=a[,c(2,1)]
a=a[!duplicated(a),]
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
a=aggregate(a[,2],by=list(a[,1]),mean)
colnames(a)=c("gene","Bartel_ActD_3T3_1")
c=merge(c,a,by=1,all=T)

a=read.delim("mouse/Zhang/1-s2.0-S2213671116302132-mmc2.csv",skip=2)
a=a[,2:3]
a[,2]=log10(a[,2]+0.1)
colnames(a)=c("gene","Hu_ActD_mESC_1")
c=merge(c,a,by=1,all=T)

nrow(c)
c=merge(c,b,by=1,all=T)
c=c[!duplicated(c), ]
c=merge(annotations[,1:2],c,by=1)
c=c[!duplicated(c), ]
nrow(c)

#in rare cases where gene maps to multiple ENSIDs, take the one with less missing info
hls = do.call(rbind, lapply(unique(c$external_gene_name), function(x){
  tmp = c[c$external_gene_name==x, ]
  if (nrow(tmp) > 1) {
    missingvals = apply(tmp,1,function(x) sum(is.na(x)) )
    tmp = tmp[order(missingvals),]
  }
  tmp[1,]
}))

writefile(hls,"all_HLs_mouse.txt")
