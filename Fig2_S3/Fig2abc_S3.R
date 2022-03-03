library(missMDA)
library(RColorBrewer)
library(plyr)
library(preprocessCore)
library(ggplot2)

c=readme("all_HLs_human.txt.gz",gzip=T)

c$Bazzini_ActD_HeLa_1 = -1*c$Bazzini_ActD_HeLa_1
c$Bazzini_ActD_RPE_1 = -1*c$Bazzini_ActD_RPE_1
c$Akimitsu2_BrU4sU_HeLa_1 = -1*c$Akimitsu2_BrU4sU_HeLa_1
c$Marks_ActD_HeLa_1 = -1*c$Marks_ActD_HeLa_1
c$Darnell_ActD_HepG2_1 = -1*c$Darnell_ActD_HepG2_1

numassays = apply(c[,3:ncol(c)],1,function(x) length(x)-sum(is.na(x)) )
nrow(c)
c=c[numassays >= 10,]
"Genes with >=10 measurements:"
nrow(c)
b=c[,1:2]
c=c[,3:ncol(c)]
c = scale(c)

ncomp <- estim_ncpPCA(c, ncp.min = 0, ncp.max = 20) #takes long time to compute
"Ncomp to use:"
ncomp$ncp

res.imp <- imputePCA(c, ncp = ncomp$ncp) #10 is output of above & stored in ncomp$ncp
c=res.imp$completeObs
header = colnames(c)
c=normalize.quantiles(c)
colnames(c) = header

writefile(cbind(b,c), gzfile("all_HLs_human_imputed.txt.gz"))

c=readme("all_HLs_human_imputed.txt.gz",gzip=T)
b=c[,1:2]
c=c[,3:ncol(c)]

c=c[,-grep("Gejman",colnames(c))]

pc = stats::prcomp(as.matrix(c))
t(t(apply(c, 2, function(x) cor(x, -1*pc$x[, 1]))))
a = as.data.frame(cbind(b,-1*pc$x[, 1]))
colnames(a)[3]="human_PC1"

writefile(a,"all_HLs_human_PC1.txt")


pc = stats::prcomp(as.matrix(t(c)))

code=as.factor(unlist(lapply(strsplit(colnames(c),"\\_"),function(x) x[1])))
type=as.factor(unlist(lapply(strsplit(colnames(c),"\\_"),function(x) x[2])))
cell=as.factor(unlist(lapply(strsplit(colnames(c),"\\_"),function(x) x[3])))
rep=as.factor(unlist(lapply(strsplit(colnames(c),"\\_"),function(x) x[4])))

samecelltype=c()
diffcelltype=c()
for (i in 1:(ncol(c)-1)){
  for (j in (i+1):ncol(c)){
    if (code[i] == code[j]) next
    if (cell[i] == cell[j]) samecelltype=c(cor(c[,i], c[,j], method='pearson'), samecelltype)
    else diffcelltype=c(cor(c[,i], c[,j], method='pearson'), diffcelltype)
  }
}

t.test(samecelltype, diffcelltype, alternative="greater")
wilcox.test(samecelltype, diffcelltype, alternative="greater")

data = data.frame(
  value=c(samecelltype,diffcelltype),
  variable=c(
    rep("same", length(samecelltype)),
    rep("diff", length(diffcelltype))
  )
)

pdf("Fig2c.pdf")
boxplot(value~variable, data=data, col="red", ylim=c(0,1))
dev.off()

## Experimented w/ weighted avg (each column weight = inverse of # of samples from study) as alternative to PC
## results were practically identical
# counts = table(code)
# weights = 1/counts[match(code, names(counts))]
# weights = weights/sum(weights)
# c[1:10,1:5]
# sapply(1:ncol(c),function(x){ c[,x]=c[,x]*weights[x] })
# weightedavg = rowSums(c)
# t(t(apply(c, 2, function(x) cor(x, weightedavg))))
# a = as.data.frame(cbind(b,weightedavg))
# colnames(a)[3]="human_avg"
# writefile(a,"all_HLs_human_weightedavg.txt")

cell
revalue(cell, c("GM07019" = "LCL")) -> cell
revalue(cell, c("GM07029" = "LCL")) -> cell
revalue(cell, c("GM10835" = "LCL")) -> cell
revalue(cell, c("GM12813" = "LCL")) -> cell
revalue(cell, c("GM12812" = "LCL")) -> cell
revalue(cell, c("GM12814" = "LCL")) -> cell
revalue(cell, c("GM12815" = "LCL")) -> cell
revalue(cell, c("GM12878" = "LCL")) -> cell
cell

jColors <- data.frame(samples = levels(cell), cols=brewer.pal(nlevels(cell), name = "Paired"))
jColors2 <- data.frame(samples = levels(type), cols=c("black","purple","blue","darkorange","darkred"))

(all = data.frame(PC1=pc$x[, 1], PC2=pc$x[, 2], type = type, cell = cell, code = code))
all = aggregate(all[,1:2], by=list(type, cell, code), mean)
all$measurement = ifelse(all$Group.1=='BrU' | all$Group.1=='BrU4sU' | all$Group.1=='4sU', "Pulse labeling", "Transcriptional shutoff")
all

#also made FigS3 with this, except including all samples (commented out line 40 which filters out Gejman samples)
pdf("Fig2ab.pdf")
par(mar=c(5,7,4,3), xpd=TRUE)
plot(pc$x[, 1], pc$x[, 2], bty="n", col = jColors2$cols[match(type, jColors2$samples)], pch=19, las=1,
xlab = sprintf("PC1 (%.2f%%)", 100*summary(pc)$importance[2,1]),
ylab = sprintf("PC2 (%.2f%%)", 100*summary(pc)$importance[2,2]))
text(pc$x[, 1]+10, pc$x[, 2], labels=code, col = jColors$cols[match(cell, jColors$samples)])
legend("topleft", inset=c(-0.3,0), cex=1, bg="white", bty="n", legend = jColors$samples, text.col=jColors$cols)
legend("bottomleft", inset=c(-0.3,0), cex=1, bg="white", bty="n", legend = jColors2$samples, pch=19, col=jColors2$cols, text.col=jColors2$cols)

boxplot(PC1~measurement, data=all, main="PC1 vs Treatment")
boxplot(PC2~measurement, data=all, main="PC2 vs Treatment")
wilcox.test(PC2~measurement, data=all)
dev.off()
