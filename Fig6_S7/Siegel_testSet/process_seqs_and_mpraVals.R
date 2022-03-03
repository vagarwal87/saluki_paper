a = read.csv("AU-Rich-Elements/sequence_level_data_Jurkat.csv")
b = read.csv("AU-Rich-Elements/sequence_level_data_Beas2B.csv")

nrow(a)
nrow(b)
a = a[complete.cases(a$ratios_T4T0_GC_resid),c("seq","ratios_T4T0_GC_resid")]
b = b[complete.cases(b$ratios_T4T0_GC_resid),c("seq","ratios_T4T0_GC_resid")]
colnames(a)[2]="Stability_Jurkat"
colnames(b)[2]="Stability_Beas2B"
nrow(a)
nrow(b)

c = merge(a,b,by=1,all=T)
cor(c[,2],c[,3],use='complete.obs')
cor(c[,2],c[,3],method='spearman',use='complete.obs')
writefile(c, "fastUTR_mpra.txt", col.names=F)
