library(biomaRt)
library(tidyr)
library(ggplot2)

ensembl = useMart(host='uswest.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl') #ensembl75 feb2014.archive.
annotations = getBM(attributes=c('ensembl_gene_id', 'external_gene_name','refseq_mrna', "chromosome_name"), mart = ensembl)
head(annotations,10)
validchr = c(1:100,'X','Y','MT')
annotations = annotations[annotations$chromosome_name %in% validchr,]

#GENE SYMBOLS
a=read.delim("human/Arango/mmc5.csv")
a=a[,1:2]
colnames(a)=c("gene","Oberdoerffer_BrU_HeLa_1") #,"Arango_HeLa_Nat10_KO"
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
b=a

a=read.delim("human/Dieterich/GSE49831_HEK293_halflives.txt", sep=' ', skip=1, header=F)
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+1)
colnames(a)=c("gene","Dieterich_4sU_HEK293_1","Dieterich_4sU_HEK293_2")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Dieterich/GSE49831_MCF7_halflives.txt", sep=' ', skip=1, header=F)
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+1)
colnames(a)=c("gene","Dieterich_4sU_MCF7_1","Dieterich_4sU_MCF7_2")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Wang/GSE49339_D-mRNA_lifetime-rep1.csv", skip=1, header=T)
a=a[,c(1,28)]
a[a[,2]<0, 2] = 0
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+1)
a=aggregate(a[,2],by=list(a[,1]),mean)
colnames(a)=c("gene","He_ActD_HeLa_1")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Wang/GSE49339_E-mRNA_Lifetime-rep2.csv", skip=1, header=T)
a=a[,c(1,28)]
a[a[,2]<0, 2] = 0
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+1)
a=aggregate(a[,2],by=list(a[,1]),mean)
colnames(a)=c("gene","He_ActD_HeLa_2")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Marks/HeLa.tsv",skip=2)
a=a[,1:2]
colnames(a)=c("gene","Marks_ActD_HeLa_1")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Marks/HepG2.tsv",skip=2) #provides data from Darnell study
a=a[,c(3,4)]
a=a[!duplicated(a), ]
a=aggregate(a[,2],by=list(a[,1]),mean)
colnames(a)=c("gene","Darnell_ActD_HepG2_1")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Mortazavi/Supplementary_Table_2.1.txt",sep=' ')
a=a[,c(1,5)]
a[,2]=log10(a[,2]+1)
a=aggregate(a[,2],by=list(a[,1]),mean)
colnames(a)=c("gene","Mortazavi_4sU_GM12878_1")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Rissland_DRUID/GSE99517_ActD.HEK.HL.one.csv",sep=',')
a=a[,1:2]
a[,2]=log10(a[,2]+0.1)
colnames(a)=c("gene","Rissland2_ActD_HEK293_1")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Rissland_DRUID/GSE99517_ActD.HEK.HL.two.csv",sep=',')
a=a[,1:2]
a[,2]=log10(a[,2]+0.1)
colnames(a)=c("gene","Rissland2_ActD_HEK293_2")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Rissland_DRUID/GSE99517_Aman.HEK.HL.one.csv",sep=',')
a=a[,1:2]
a[,2]=log10(a[,2]+0.1)
colnames(a)=c("gene","Rissland2_Aman_HEK293_1")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Rissland_DRUID/GSE99517_Aman.HEK.HL.two.csv",sep=',')
a=a[,1:2]
a[,2]=log10(a[,2]+0.1)
colnames(a)=c("gene","Rissland2_Aman_HEK293_2")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Rissland_DRUID/GSE99517_HEK.4SU.HL.one.csv",sep=',')
a=a[,1:2]
a[,2]=log10(a[,2]+0.1)
colnames(a)=c("gene","Rissland2_4sU_HEK293_1")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Rissland_DRUID/GSE99517_HEK.4SU.HL.two.csv",sep=',')
a=a[,1:2]
a[,2]=log10(a[,2]+0.1)
colnames(a)=c("gene","Rissland2_4sU_HEK293_2")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Zimmer/half_life_by_gene.txt")
a[,2]=log10(a[,2]+1)
a=aggregate(a[,2],by=list(a[,1]),mean)
colnames(a)=c("gene","Zimmer_ActD_Bcell_1")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Duan/LCL_HFs.csv",skip=1,sep=',')
a=a[,c(2,grep("GM", colnames(a)))]
header=colnames(a)
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
a=aggregate(a[,2:ncol(a)],by=list(a[,1]),mean)
colnames(a)=header
colnames(a)=c("gene",paste0("Gejman_4sU_",colnames(a)[2:ncol(a)]))
b=merge(b,a,by=1,all=T)

a=read.delim("human/Schoefield/Schofield_K562_half_lives.txt")
a=a[,c(1,4,5)]
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
colnames(a)=c("gene","Simon_4sU_K562_1","Simon_4sU_K562_2")
b=merge(b,a,by=1,all=T)

a=read.delim("human/Rissland_ORFeome/Supplemental_Table_S1.csv",skip=9,sep=',')
a=a[,1:(ncol(a)-1)]
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
colnames(a)=c("gene","Rissland_4sU_HEK293_1","Rissland_4sU_HEK293_2","Rissland_4sU_HEK293_3","Rissland_4sU_HEK293_4")
b=merge(b,a,by=1,all=T)

b=merge(b,annotations[,1:2],by.x=1,by.y=2,all=F)
b$gene=NULL
b=aggregate(b[,1:(ncol(b)-1)],by=list(b$ensembl_gene_id),function(x){ mean(x, na.rm=T) } )
colnames(b)[1]=c("gene")
nrow(b)

#ENSEMBL GENE ID
a=read.delim("human/Bazzini_ORFome/293T-endog.csv")
a=a[,c(3,5)]
colnames(a)=c("gene","Bazzini_ActD_HEK293_1")
c=a

a=read.delim("human/Bazzini_ORFome/HeLa-endog.csv")
a=a[,3:4]
colnames(a)=c("gene","Bazzini_ActD_HeLa_1")
c=merge(c,a,by=1,all=T)

a=read.delim("human/Bazzini_ORFome/RPE-endog.csv")
a=a[,3:4]
colnames(a)=c("gene","Bazzini_ActD_RPE_1")
c=merge(c,a,by=1,all=T)

a=read.delim("human/Bazzini_ORFome/K562-SLAMseq.csv")
a=a[,c(3,5)]
a[,2]=log10(1/-a[,2]+1)
colnames(a)=c("gene","Bazzini_4sU_K562_1")
c=merge(c,a,by=1,all=T)

#NM IDs
a=read.delim("human/Tani/Tani_Supp_Tables_revised2.lineSep.tsv",skip=2)
a=merge(annotations,a,by.x=3,by.y=1)
a=a[,c(2,7)]
a=a[!duplicated(a), ]
a[,2]=log10(a[,2]+0.1)
a=aggregate(a[,2],by=list(a[,1]),mean)
colnames(a)=c("gene","Akimitsu_BrU_HeLa_1")
c=merge(c,a,by=1,all=T)

#ENSEMBL GENE ID
a=read.delim("human/Rinn/Supplemental_Table_S10_Rinn.txt")
a=a[,c(1,grep("halflife", colnames(a)))]
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
a[,1]=substring(a[,1],1,15)
colnames(a)=c("gene","Rinn_ActD_K562_1","Rinn_ActD_K562_2","Rinn_ActD_K562_3","Rinn_ActD_H1ESC_1","Rinn_ActD_H1ESC_2","Rinn_ActD_H1ESC_3")
c=merge(c,a,by=1,all=T)

a=read.delim("human/Akimitsu/media-3.csv",F,skip=3)
a=a[,c(1,4)]
a[,2]=log10(a[,2]+0.00001)
colnames(a)=c("gene","Akimitsu2_BrU4sU_HeLa_1")
c=merge(c,a,by=1,all=T)

a=read.delim("human/Jaffrey/hek293_halflife.txt")
a=a[,c(1,3)]
a[,2]=log10(a[,2]+0.1)
colnames(a)=c("gene","Jaffrey_ActD_HEK293_1")
c=merge(c,a,by=1,all=T)

a=read.delim("human/Cramer/elife-45056-supp3-v2")
a=a[,c(1,5)]
a[,2]=log10(a[,2]+1)
a[,1]=substring(a[,1],1,15)
colnames(a)=c("gene","Cramer_4sU_K562_1")
c=merge(c,a,by=1,all=T)

a=read.delim("human/Shendure/41587_2020_480_MOESM3_ESM.csv",skip=1)
a=a[,c(1,4)]
a[,2]=log10(a[,2]+0.1)
a[,1]=substring(a[,1],1,15)
colnames(a)=c("gene","Shendure_4sU_A549_1")
c=merge(c,a,by=1,all=T)

a=read.delim("human/Cramer_TTseq/reciprocal_best_noAmbig.txt")
a[,2]=log10(a[,2]+1)
colnames(a)=c("gene","Cramer2_4sU_K562_1")
c=merge(c,a,by=1,all=T)

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

writefile(hls,"all_HLs_human.txt")
