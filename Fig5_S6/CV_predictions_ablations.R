library(rhdf5)
library(ggplot2)

file = args[1]

getvals <- function(dir, pattern, filter){
	files = list.files(path=dir, pattern=pattern, full.names=T, recursive=T)
  files = files[grep(filter, files)]
	table <- as.data.frame(do.call("rbind", lapply(files, FUN=function(file){
    if (grepl("preds", file)) { vals = t(h5read(file,"preds")) } #(1297, 1)
    if (grepl("targets", file)) {
      vals = h5read(file,"targets") #(1297, 1, 1)
      vals = vals[1,1,]
    }
    fold = strsplit(file, "/")[[1]][2]
    foldrun = strsplit(fold, "_")[[1]]
		tmp = data.frame(vals = vals, fold=foldrun[1], run=foldrun[2])
    tmp$seqnum = 1:nrow(tmp)
		tmp
	})))
  table=aggregate(table$vals, list(table$fold, table$seqnum), mean)
	table
}

findcorrs = function(track, species){
	preds = getvals(track, "preds.h5", species)
	colnames(preds)[3]="preds"
	targets = getvals(track, "targets.h5", species)
	colnames(targets)[3]="targets"

	all=merge(preds, targets, by=c(1,2))
	corrs = do.call(rbind, lapply(unique(all$Group.1), function(f){
		idx = all$Group.1==f
		data.frame(fold=f, Pearson = cor(all$preds[idx], all$targets[idx]))
	}))
	print(corrs)
}

"Human"
"Sequence+Splice only"
findcorrs("train_woc", "test0")
"Sequence+Codon only"
findcorrs("train_wos", "test0")
"Sequence only"
findcorrs("train_wosc", "test0")

"Mouse"
"Sequence+Splice only"
findcorrs("train_woc", "test1")
"Sequence+Codon only"
findcorrs("train_wos", "test1")
"Sequence only"
findcorrs("train_wosc", "test1")
