library(devtools)
library(RColorBrewer)
library(plyr)

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

c = readme("all_HLs_mouse.txt.gz", gzip=T)
c$external_gene_name=NULL

code=as.factor(unlist(lapply(strsplit(colnames(c[,2:ncol(c)]),"\\_"),function(x) x[1]))) #c(, "PC1")
type=as.factor(unlist(lapply(strsplit(colnames(c[,2:ncol(c)]),"\\_"),function(x) x[2])))
cell=as.factor(unlist(lapply(strsplit(colnames(c[,2:ncol(c)]),"\\_"),function(x) x[3])))
rep=as.factor(unlist(lapply(strsplit(colnames(c[,2:ncol(c)]),"\\_"),function(x) x[4])))

code

cell
revalue(cell, c("E14ESC" = "mESC")) -> cell
cell

colorbins=seq(0,1,0.01)

a=round(abs(cor(c[,2:ncol(c)],use='pairwise.complete.obs', method='spearman')),2)
b=round(abs(cor(c[,2:ncol(c)],use='pairwise.complete.obs', method='pearson')),2)

type.color.map = c("black","purple","blue")
cell.color.map = brewer.pal(nlevels(cell), name = "Paired")
cor.color.map = colorRampPalette(brewer.pal(n=9, name='Purples'))(100)

rlab=rbind(cell.color.map[ cell ], type.color.map[ type ])

heatmaps3=function(x) {
  heatmap.3(
    as.matrix(x),density.info="none",
    trace="none", breaks=colorbins, symkey=FALSE, cexCol=1.4, labRow=code, labCol=code,
    cexRow=1.4, scale="none", dendrogram="row", key=TRUE, RowSideColors=rlab,
    Rowv=T, Colv=T, col=matlab::jet.colors(length(colorbins)-1)
  )
  legend("topright",legend=c(levels(type),"",levels(cell)),
    fill=c(type.color.map,"white",cell.color.map), border=FALSE, bty="n", y.intersp = 0.7, cex=1)
}

nrow(c)
pdf("FigS2.pdf",height=20,width=20)
par(mar=c(10,4,4,10), oma=c(5,1,1,5))
heatmaps3(a)
heatmaps3(b)
z=apply(c,2,function(x) sum(!is.na(x)) )
barplot(z,las=2)
dev.off()
