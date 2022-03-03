library(ggplot2)

getvals <- function(dir, pattern, filter){
    files = list.files(path=dir, pattern=paste0('^',pattern), full.names=T, recursive=T)
    files = files[grep(filter, files)]
    say(files)
    table <- as.data.frame(do.call("rbind", lapply(files, FUN=function(file){
      vals=readme(file,skip=95,header=F,gzip=T)
      fold = strsplit(file, "/")[[1]][2]
      foldrun = strsplit(fold, "_")[[1]]
      tmp = cbind(data.frame(fold=foldrun[2], run=foldrun[3]), vals)
      tmp
    })))
    table=aggregate(table[6:ncol(table)], list(table$fold, table$V1, table$V3), mean)
    table
}

a = getvals("spikeinmotif", "3utr.*txt", "")
b = getvals("spikeinmotif", "splice3utr.*txt", "")
a = rbind(a, b)
a[,5:54]=a[,5:54]-a[,4]
b=aggregate(a[,5:54], list(a$Group.1,a$Group.3), mean)
utr3=aggregate(b[,3:ncol(b)], list(b$Group.2), mean)

a = getvals("spikeinmotif", "orf.*txt", "")
b = getvals("spikeinmotif", "spliceorf.*txt", "")
a = rbind(a, b)
a[,5:54]=a[,5:54]-a[,4]
b=aggregate(a[,5:ncol(a)], list(a$Group.1,a$Group.3), mean)
orf=aggregate(b[,3:ncol(b)], list(b$Group.2), mean)

a = getvals("spikeinmotif", "5utr.*txt", "")
b = getvals("spikeinmotif", "splice5utr.*txt", "")
a = rbind(a, b)
a[,5:54]=a[,5:54]-a[,4]
b=aggregate(a[,5:54], list(a$Group.1,a$Group.3), mean)
utr5=aggregate(b[,3:ncol(b)], list(b$Group.2), mean)

full=cbind(utr5, orf[,2:ncol(orf)], utr3[,2:ncol(utr3)])

full=full[full$Group.1 != 'ATTATTA',]

pdf("png/Fig5F.pdf")
matplot(t(as.matrix(full[,2:ncol(full)])), type = "l", lty = 1, lwd = 2)
abline(v=50)
abline(v=100)
legend("bottomleft", full$Group.1, col = seq_len(nrow(full)), text.col=seq_len(nrow(full)), cex=0.8, lwd=2)
dev.off()

a = getvals("spikeinmotif", "codon.*txt", "")
a[,5:54]=a[,5:54]-a[,4]
b=aggregate(a[,5:54], list(a$Group.1,a$Group.3), mean)
codon=aggregate(b[,3:ncol(b)], list(b$Group.2), mean)

cols = rep("grey",61)
cols[which(codon$Group.1=="GGT")]="brown"
cols[which(codon$Group.1=="GGG")]="black"
cols[which(codon$Group.1=="GCT")]="blue"
cols[which(codon$Group.1=="TTC")]="cyan"
cols[which(codon$Group.1=="CTG")]="magenta"
cols[which(codon$Group.1=="GAA")]="green"
cols[which(codon$Group.1=="GAC")]="orange"
cols[which(codon$Group.1=="TGT")]="darkred"
cols[which(codon$Group.1=="CGA")]="red"
cols[which(codon$Group.1=="CGG")]="lightblue"

pdf("png/Fig5G.pdf", width=10, height=5)
matplot(t(as.matrix(codon[,2:ncol(codon)])), type = "l", lty = 1, lwd = 2, col=cols)
abline(v=50)
abline(v=100)
legend("top", codon$Group.1[cols!="grey"], col = cols[cols!="grey"], text.col=cols[cols!="grey"], cex=0.8, lwd=2, ncol=10)
dev.off()

pdf("png/Fig5H.pdf")
a=read.delim("codon_stability_coefs_Forrest_et_al.txt")
"codon"
vals = data.frame(codon=codon$Group.1, meanVal = apply(codon[,2:ncol(codon)], 1, mean))
vals = merge(vals,a,by=1)
head(vals)
vals$meanCSC = apply(vals[,3:ncol(vals)], 1, mean)
plot(vals$meanVal, vals$meanCSC, xlab="Saluki mean codon influence", ylab="CSC (Forrest et al)", bty='n', las=1, pch=19)
dev.off()

cor(vals$meanVal,vals[,3:ncol(vals)])
cor(vals$meanVal,vals[,3:ncol(vals)], method='spearman')
