library(rhdf5)
library(ggplot2)
library(LSD)

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

"Human"
preds = getvals("train_gru", "preds.h5", "test0")
colnames(preds)[3]="preds"
targets = getvals("train_gru", "targets.h5", "test0")
colnames(targets)[3]="targets"

all=merge(preds, targets, by=c(1,2))
head(all)
corrs = do.call(rbind, lapply(unique(all$Group.1), function(f){
	idx = all$Group.1==f
	data.frame(fold=f, Pearson = cor(all$preds[idx], all$targets[idx]))
}))

corrs

cor(all$preds, all$targets)
cor(all$preds, all$targets, method='spearman')

# Calculate R^2
avr_y_actual <- mean(all$targets)
ss_total <- sum((all$targets - avr_y_actual)^2)
ss_regression <- sum((all$preds - avr_y_actual)^2)
ss_residuals <- sum((all$targets - all$preds)^2)
"R^2:"
(r2 <- 1 - ss_residuals / ss_total)

pdf("png/DNN_vs_HL_human.pdf") #Fig5c
heatscatter(all$preds, all$targets, xlab="CV fitted model prediction", ylab="Half-life", bty='n', cex=0.3, las=1, ylim=c(-4,4))
dev.off()

"Mouse"
preds = getvals("train_gru", "preds.h5", "test1")
colnames(preds)[3]="preds"
targets = getvals("train_gru", "targets.h5", "test1")
colnames(targets)[3]="targets"

all=merge(preds, targets, by=c(1,2))
head(all)
corrs = do.call(rbind, lapply(unique(all$Group.1), function(f){
	idx = all$Group.1==f
	data.frame(fold=f, Pearson = cor(all$preds[idx], all$targets[idx]))
}))

corrs

cor(all$preds, all$targets)
cor(all$preds, all$targets, method='spearman')

pdf("png/DNN_vs_HL_mouse.pdf") #FigS6c
heatscatter(all$preds, all$targets, xlab="CV fitted model prediction", ylab="Half-life", bty='n', cex=0.3, las=1, ylim=c(-4,4))
dev.off()
