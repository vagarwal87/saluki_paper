library(LSD)

getvals <- function(dir, pattern, filter, refValt=F){
    files = list.files(path=dir, pattern=paste0('^',pattern), full.names=T, recursive=T)
    files = files[grep(filter, files)]
    say(files)
    table <- as.data.frame(do.call("rbind", lapply(files, FUN=function(file){
      vals=readme(file,skip=95,header=T,gzip=T)
      fold = strsplit(file, "/")[[1]][2]
      foldrun = strsplit(fold, "_")[[1]]
      if(refValt){
        idx1=seq(1,nrow(vals),2)
        idx2=seq(2,nrow(vals),2)
        tmp = cbind(data.frame(fold=foldrun[1], run=foldrun[2]), seq=vals[idx1,1], pred=vals[idx1,2]-vals[idx2,2])
      }
      else {
        tmp = cbind(data.frame(fold=foldrun[1], run=foldrun[2]), vals)
      }
      tmp
    })))
    table=aggregate(table$pred, list(table$seq), mean)
    table
}

a = getvals("fastUTR_predictions", ".*txt", "")
nrow(a)
colnames(a)=c("seq","Prediction")
head(a)

writefile(a,"fastUTR_MPRA_predictions.txt")

a=read.delim("fastUTR_MPRA_predictions.txt")
b=read.delim("fastUTR_mpra.txt",F)
colnames(b)[2:3]=c("mpra_val_Jurkat", "mpra_val_Beas2B")
head(b)
data=merge(b,a,by=1)

head(data)
data$mean=apply(data[,2:3],1,function(x){ mean(x, na.rm=T) })

round(cor(data[,2:5], use='pairwise.complete.obs'),2)
round(cor(data[,2:5], method='spearman', use='pairwise.complete.obs'),2)

pdf("fastUTR_scatter.pdf")  #FigS7abc
heatscatter(data$Prediction[!is.na(data[,2])], data[!is.na(data[,2]),2], xlab="Predicted 3' UTR effect", ylab="Observed 3' UTR effect, Jurkat", bty='n', las=1, xlim=c(-2,0.5), ylim=c(0,1))
heatscatter(data$Prediction[!is.na(data[,3])], data[!is.na(data[,3]),3], xlab="Predicted 3' UTR effect", ylab="Observed 3' UTR effect, Beas2B", bty='n', las=1, xlim=c(-2,0.5), ylim=c(0,0.8))
heatscatter(data[complete.cases(data),2], data[complete.cases(data),3], xlab="Observed 3' UTR effect, Jurkat", ylab="Observed 3' UTR effect, Beas2B", bty='n', las=1, xlim=c(0,0.8), ylim=c(0,0.8))
dev.off()


plotscatter = function(celltype){
  a = getvals(paste0("fastUTR_",celltype,"_predictions"), ".*txt", "", refValt=T)
  colnames(a)=c("seq","Prediction")

  writefile(a,paste0("fastUTR_MPRA_predictions_",celltype,".txt"))

  b=read.delim(paste0("fastUTR_mpra_",celltype,".txt"),F)
  idx1=seq(1,nrow(b),2)
  b=b[idx1,]

  colnames(b)=c("seq","mpra_val")
  data=merge(b,a,by=1)
  head(data)

  ## THIS CODE EMULATES THE FILTERING DONE IN THE PAPER ITSELF, FOR BENCHMARKING PURPOSES
  # filterToAU = read.csv(paste0("AU-Rich-Elements/sequence_level_data_",celltype,".csv"))
  # filterToAU = filterToAU[filterToAU$parent_motif_AREs_numbered_BakheetPlus != 0 &
  #                         filterToAU$parent_motif_AREs_numbered_BakheetPlus < 7 &
  #                         filterToAU$iscontrol!=1 &
  #                         complete.cases(filterToAU[,c("effect_size_T4T0","effect_size_T4T0_GC_resid")]),
  #                         "seq"]
  # say(nrow(data))
  # data = data[data$seq %in% filterToAU,]
  # say(nrow(data))

  say(cor(data$Prediction, data$mpra_val))
  say(cor(data$Prediction, data$mpra_val, method='spearman'))

  pdf(paste0("fastUTR_scatter_",celltype,"Variants.pdf"))
  heatscatter(data$Prediction, data$mpra_val, xlab="Predicted 3' UTR effect", ylab=paste0("Observed 3' UTR variant effect, ", celltype), bty='n', las=1, xlim=c(-0.5,1.5), ylim=c(-0.5,0.5))
  dev.off()

  data[,c("seq","mpra_val")]
}

jurkat = plotscatter("Jurkat") #FigS7d
beas2b = plotscatter("Beas2B") #Fig6e & S7e (identical)

data = merge(jurkat, beas2b, by=1)

cor(data[,2], data[,3])
cor(data[,2], data[,3], method='spearman')

pdf("fastUTR_scatter_JurkatVsBeas2b_Variants.pdf") #FigS7f
heatscatter(data[,2], data[,3], xlab="Jurkat variant effect", ylab="Beas2B variant effect", bty='n', las=1, xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
dev.off()
