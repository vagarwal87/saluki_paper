library(missMDA)
library(RColorBrewer)
library(plyr)
library(preprocessCore)
library(ggplot2)

c=readme("all_HLs_mouse.txt.gz",gzip=T)

numassays = apply(c[,3:ncol(c)],1,function(x) length(x)-sum(is.na(x)) )
nrow(c)
c=c[numassays >= 5,]
"Genes with >=5 measurements:"
nrow(c)
b=c[,1:2]
c=c[,3:ncol(c)]
c = scale(c)

ncomp <- estim_ncpPCA(c, ncp.min = 0, ncp.max = 10) #takes long time to compute
"Ncomp to use:"
ncomp$ncp

res.imp <- imputePCA(c, ncp = ncomp$ncp) #4 is output of above & stored in ncomp$ncp
c=res.imp$completeObs
header = colnames(c)
c=normalize.quantiles(c)
colnames(c) = header

writefile(cbind(b,c), gzfile("all_HLs_mouse_imputed.txt.gz"))

pc = stats::prcomp(as.matrix(c))
t(t(apply(c, 2, function(x) cor(x, -1*pc$x[, 1]))))

a = as.data.frame(cbind(b,-1*pc$x[, 1]))
colnames(a)[3]="mouse_PC1"
writefile(a,"all_HLs_mouse_PC1.txt")


c=readme("all_HLs_mouse_imputed.txt.gz",gzip=T)
b=c[,1:2]
c=c[,3:ncol(c)]

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

pdf("Fig2e.pdf")
boxplot(value~variable, data=data, col="red")
dev.off()

cell
revalue(cell, c("E14ESC" = "mESC")) -> cell
cell

jColors <- data.frame(samples = levels(cell), cols=brewer.pal(nlevels(cell), name = "Paired"))
jColors2 <- data.frame(samples = levels(type), cols=c("black","purple","blue"))

pdf("Fig2d.pdf")
par(mar=c(5,7,4,3), xpd=TRUE)
plot(pc$x[, 1], pc$x[, 2], bty="n", col = jColors2$cols[match(type, jColors2$samples)], pch=19, las=1,
xlab = sprintf("PC1 (%.2f%%)", 100*summary(pc)$importance[2,1]),
ylab = sprintf("PC2 (%.2f%%)", 100*summary(pc)$importance[2,2])) #, col="blue" , axis=F, ann=F, xlim=c(-20,50), ylim=c(-40,20)
text(pc$x[, 1]+10, pc$x[, 2], labels=code, col = jColors$cols[match(cell, jColors$samples)])
legend("topleft", inset=c(-0.3,0), cex=1, bg="white", bty="n", legend = jColors$samples, text.col=jColors$cols)
legend("bottomleft", inset=c(-0.3,0), cex=1, bg="white", bty="n", legend = jColors2$samples, pch=19, col=jColors2$cols, text.col=jColors2$cols)
dev.off()
