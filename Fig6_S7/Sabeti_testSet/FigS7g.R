library(psych)

a=readme("fastUTR_mpra_with_predictions.txt")
idx = grep("Ref",colnames(a))
a=a[,idx]
head(a)
a$Ref_GC = NULL
colnames(a) = gsub("log2FoldChange_Ref_","", colnames(a))
head(a)

png("all_spearman_Ref.png",width=800,height=800,type="cairo") #Figure S7g
pairs.panels(a, method = "spearman", hist.col = "#00AFBB", density = T, scale=T, lm=T, cex.cor = 2, cex.axis = 1.6, las = 2)
dev.off()
