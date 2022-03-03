library(LSD)

getvals <- function(dir, pattern, filter, celltype){
    files = list.files(path=dir, pattern=paste0('^',pattern), full.names=T, recursive=T)
    files = files[grep(filter, files)]
    names = read.delim(paste0(celltype,".names.txt"),F)
    say(files)
    table <- as.data.frame(do.call("rbind", lapply(files, FUN=function(file){
      vals=readme(file,skip=95,header=T,gzip=T)
      fold = strsplit(file, "/")[[1]][2]
      foldrun = strsplit(fold, "_")[[1]]
      tmp = cbind(data.frame(fold=foldrun[1], run=foldrun[2]), seqids=names$V1, vals)
      tmp
    })))
    table=aggregate(table$pred, list(table$seqids), mean)
    table
}

plotscatter = function(celltype){
  a = getvals(paste0("fastUTR_",celltype,"_predictions"), ".*txt", "", celltype)
  colnames(a)=c("seqid","Prediction")

  writefile(a,paste0("fastUTR_MPRA_predictions_",celltype,".txt"))
}

plotscatter("ENCFF090JTW")
plotscatter("ENCFF770UJN")
