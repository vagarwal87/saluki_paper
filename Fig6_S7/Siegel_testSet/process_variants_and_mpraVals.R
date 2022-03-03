library(gdata)

prepare = function(x){
  a = read.csv(paste0("AU-Rich-Elements/sequence_level_data_",x,".csv"))
  b = a[a$iscontrol!=1 & complete.cases(a[,c("effect_size_T4T0","effect_size_T4T0_GC_resid")]),
                          c("region","seq","effect_size_T4T0_GC_resid")]
  a = a[a$iscontrol==1,c("region","seq")]
  c = data.frame(seq=a$seq[match(b$region,a$region)])
  b$region = NULL
  c$effect_size_T4T0_GC_resid = 0
  b = interleave(b, c)
  colnames(b)[2]=paste0("DeltaStability_",x)
  b
}

a = prepare("Jurkat")
b = prepare("Beas2B")
# writefile(a, "fastUTR_mpra_Jurkat.txt", col.names=F)
# writefile(b, "fastUTR_mpra_Beas2B.txt", col.names=F)

c = merge(a,b,by=1,all=T)
nrow(c)
c = c[complete.cases(c),]
nrow(c)
cor(c[,2],c[,3],use='complete.obs')
cor(c[,2],c[,3],method='spearman',use='complete.obs')
# writefile(c, "fastUTR_mpra.txt", col.names=F)
