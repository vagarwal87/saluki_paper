library(psych)

#FOR Ref (pred vs obs)
a=readme("fastUTR_mpra_with_predictions.txt")
idx = grep("Ref",colnames(a))
a=a[,idx] #change all Skew to Ref or Alt
# idx = grep("log2FoldChange_Ref_",colnames(a))
# a$log2FoldChange_Ref_Mean=apply(a[,idx],1,mean) #change all Skew to Ref or Alt
head(a)
a$Ref_GC = NULL
colnames(a) = gsub("log2FoldChange_Ref_","", colnames(a))
head(a)
#
# lapply(3:ncol(a), function(x){
#   print(colnames(a)[x])
#   print(summary(lm(a[,x]~Ref_pred,data=a)))
#   print(summary(lm(a[,x]~Ref_GC,data=a)))
#   print(summary(lm(a[,x]~Ref_GC+Ref_pred,data=a)))
# })

# png("all_pearson_Ref.png",width=800,height=800,type="cairo")
# pairs.panels(a, method = "pearson", hist.col = "#00AFBB", density = T, scale=T, lm=T, cex.cor = 2, cex.axis = 1.6, las = 2)
# dev.off()
png("all_spearman_Ref.png",width=800,height=800,type="cairo")
pairs.panels(a, method = "spearman", hist.col = "#00AFBB", density = T, scale=T, lm=T, cex.cor = 2, cex.axis = 1.6, las = 2)
dev.off()
#
# #FOR Alt (pred vs obs)
# a=readme("fastUTR_mpra_with_predictions.txt")
# idx = grep("Alt",colnames(a))
# a=a[,idx] #change all Skew to Ref or Alt
# idx = grep("log2FoldChange_Alt_",colnames(a))
# a$log2FoldChange_Alt_Mean=apply(a[,idx],1,mean) #change all Skew to Ref or Alt
# colnames(a) = gsub("log2FoldChange_Alt_","", colnames(a))
# head(a)
#
# png("all_pearson_Alt.png",width=800,height=800,type="cairo")
# pairs.panels(a, method = "pearson", hist.col = "#00AFBB", density = T, scale=T, lm=T, cex.cor = 2, cex.axis = 1.6, las = 2)
# dev.off()
# png("all_spearman_Alt.png",width=800,height=800,type="cairo")
# pairs.panels(a, method = "spearman", hist.col = "#00AFBB", density = T, scale=T, lm=T, cex.cor = 2, cex.axis = 1.6, las = 2)
# dev.off()
#
# #FOR SKEW (pred vs obs)
# a=readme("fastUTR_mpra_with_predictions.txt")
# a$Skew_pred = a$Ref_pred - a$Alt_pred
# a$Skew_GC = a$Ref_GC - a$Alt_GC
# idx = grep("Skew",colnames(a))
# a=a[,idx] #change all Skew to Ref or Alt
# idx = grep("log2FoldChange_Skew_",colnames(a))
# a$log2FoldChange_Skew_Mean=apply(a[,idx],1,mean) #change all Skew to Ref or Alt
# colnames(a) = gsub("log2FoldChange_Skew_","", colnames(a))
# head(a)
#
# png("all_pearson_Skew.png",width=800,height=800,type="cairo")
# pairs.panels(a, method = "pearson", hist.col = "#00AFBB", density = T, scale=T, lm=T, cex.cor = 2, cex.axis = 1.6, las = 2)
# dev.off()
# png("all_spearman_Skew.png",width=800,height=800,type="cairo")
# pairs.panels(a, method = "spearman", hist.col = "#00AFBB", density = T, scale=T, lm=T, cex.cor = 2, cex.axis = 1.6, las = 2)
# dev.off()
