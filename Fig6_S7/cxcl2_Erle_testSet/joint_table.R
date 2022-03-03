a=readme("cxcl2_deltaA.wig",skip=2,header=F, sep="\t")
c=readme("cxcl2_deltaC.wig",skip=2,header=F)
g=readme("cxcl2_deltaG.wig",skip=2,header=F)
u=readme("cxcl2_deltaU.wig",skip=2,header=F)

all=cbind(a,c[,2],g[,2],u[,2])
colnames(all)=c("pos","A","C","G","T")

all
writefile(all,"mutagenesis_MPRA.txt")
