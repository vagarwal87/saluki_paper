library(LSD)

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
    table=aggregate(table$pred, list(table$seq), mean)
    table
}

a = getvals("fastUTR_predictions", ".*txt", "")
nrow(a)
head(a)

colnames(a)=c("seq","Prediction")
writefile(a,"fastUTR_MPRA_predictions.txt")

b=read.delim("fastUTR_mpra.txt",F)
colnames(b)[2]="mpra_val"
head(b)
data=merge(b,a,by=1)
cor(data$Prediction,data$mpra_val)
cor(data$Prediction,data$mpra_val,method='spearman')

pdf("fastUTR_scatter.pdf") #Fig6d
heatscatter(data$Prediction, data$mpra_val, xlab="Predicted 3' UTR effect", ylab="Observed 3' UTR effect", bty='n', las=1, ylim=c(1,5))
dev.off()
