library(ggplot2)
library(tidyr)

getvals <- function(dir, pattern, filter){
    files = list.files(path=dir, pattern=paste0('^',pattern), full.names=T, recursive=T)
    files = files[grep(filter, files)]
    say(files)
    table <- as.data.frame(do.call("rbind", lapply(files, FUN=function(file){
      vals=readme(file,skip=95,header=T,gzip=T)
      fold = strsplit(file, "/")[[1]][2]
      foldrun = strsplit(fold, "_")[[1]]
      tmp = cbind(data.frame(fold=foldrun[1], run=foldrun[2]), vals)
      tmp
    })))
    table=aggregate(table$pred, list(table$pos, table$mut), mean)
    table
}

a = getvals("mutagenesis_predictions", ".*txt", "")
nrow(a)
nrow(a)
head(a)

colnames(a)=c("pos","Alt","Prediction")
writefile(spread(a, Alt, Prediction),"mutagenesis_MPRA_predictions.txt")

a = a[a$Prediction!=0,]

b=read.delim("mutagenesis_MPRA.txt")
colnames(b)[5] = 'T'
b <- gather(b, Alt, mpra_val, A:T, factor_key=TRUE)
b

data=merge(b,a,by=1:2)
colnames(data)[c(1,4)] = c("Position","Prediction")
cor(data$Prediction,data$mpra_val)
cor(data$Prediction,data$mpra_val,method='spearman')

pdf("cxcl_mutagenesis_scatter.pdf") #Fig6c
plot(data$Prediction,data$mpra_val, xlab="Predicted variant effect", ylab="Observed variant effect", bty='n', las=1, pch=19)
dev.off()

# pdf("cxcl_mutagenesis.pdf")
# par(mfrow=c(2,1)) #, mar=c(8,8,5,5), mgp = c(4, 3, 2), xpd=TRUE
# ticknum=300 #step by 300/3, or 100 nt
#
# mycolors <- data.frame(nuc=c("A", "C", "T", "G"),
#                            color=c("green", "blue", "orange", "red"))
#
# lab <- rep("", nrow(data))
# idxs=seq(1, nrow(data), by=ticknum)
# data$Position=-1*data$Position
# data=data[order(data$Position),]
# data
# lab[idxs] <- data$Position[idxs]-data$Position[1]
# thesecols=unlist(lapply(data$Alt,function(x) mycolors[x, "color"]))
# df.bar <- barplot(data$mpra_val,col=ifelse(data$mpra_val<0,"black","black"),border=F, names.arg=lab, las=2,cex.names=0.5,
#   ylim=c(min(data$mpra_val)-0.2,max(data$mpra_val)+0.2))
# points(x = df.bar, y = data$mpra_val, col=as.character(thesecols), pch=20, cex=0.8)
# df.bar <- barplot(data$Prediction,col=ifelse(data$Prediction<0,"black","black"),border=F, names.arg=lab, las=2,cex.names=0.5,
#   ylim=c(min(data$Prediction)-0.2,max(data$Prediction)+0.2))
# points(x = df.bar, y = data$Prediction, col=as.character(thesecols), pch=20, cex=0.8)
# dev.off()
