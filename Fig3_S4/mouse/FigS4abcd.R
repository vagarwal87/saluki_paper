library(ggplot2)
library(glmnet)
library(LSD)

a=read.csv("experimentResults.csv")
round(a,3)

results = data.frame(tests = c("B vs B5","B vs BO","B vs B3","B vs BC","BC vs BCO","BC vs B5C",
"B3 vs B3M","BC vs BC3M","BC3M vs BC3MS","BC3M vs BC3MD","BC3MS vs BC3MSD"), pvals =
c(t.test(a$B, a$B5, paired=T,alternative='less')$p.value,
t.test(a$B, a$BO, paired=T,alternative='less')$p.value,
t.test(a$B, a$B3, paired=T,alternative='less')$p.value,
t.test(a$B, a$BC, paired=T,alternative='less')$p.value,
t.test(a$BC, a$BCO, paired=T,alternative='less')$p.value,
t.test(a$BC, a$B5C, paired=T,alternative='less')$p.value,
t.test(a$B3, a$B3M, paired=T,alternative='less')$p.value,
t.test(a$BC, a$BC3M, paired=T,alternative='less')$p.value,
t.test(a$BC3M, a$BC3MS, paired=T,alternative='less')$p.value,
t.test(a$BC3M, a$BC3MD, paired=T,alternative='less')$p.value,
t.test(a$BC3MS, a$BC3MSD, paired=T,alternative='less')$p.value)
)

results$p.corrected = pmin(results$pvals*nrow(results),1)
results

meltData <- reshape::melt(a)

pdf("png/violinplots_mouse.pdf",width=12,height=4) #FigS4a.pdf
ggplot(meltData, aes(x=variable, y=value)) + geom_violin()
p<-ggplot(meltData, aes(x=variable, y=value)) + geom_violin(position=position_dodge(1))
print(p)
dev.off()

names = c("BC3MSD")

topcoefs = unlist(lapply( c(names), function(sample) {
    topN=30
    load(paste("Robj/", sample,"_CV-Lasso.Robj",sep=''))
    coef.exact = coef(cvfit, s = "lambda.min", exact = TRUE)
    coefficients <- as.numeric(coef.exact)
    names(coefficients) <- coef.exact@Dimnames[[1]]
    names(coefficients)[1:9]=paste0("Basic.", names(coefficients)[1:9])
    select <- coefficients[coefficients!=0]
    select[order(abs(select),decreasing=T)]
    reorder <- select[order(abs(select),decreasing=T)]
    names(reorder)=gsub("SeqWeaver","SeqWeaver.",names(reorder))
    names(reorder)=gsub("DeepRiPe","DeepRiPe.",names(reorder))

    topcoefs=names(reorder[1:topN])
    type = as.factor(unlist(lapply(names(reorder), function(x) strsplit(x, '\\.')[[1]][1])))
    names(reorder) = unlist(lapply(names(reorder), function(x) paste0(strsplit(x, '\\.')[[1]][-1],collapse='.')))
    cols = c("darkorange", "red", "purple", "green", "cyan", "blue")[type]
    # cols = c("black")

    pdf(paste("png/", sample,"_CV-Lasso-coefficients_mouse.pdf",sep=''),width=8,height=8) #FigS4c.pdf
    par(mar = c(5, 6, 4, 2) + 0.1, oma= c(0.2, 0, 0, 0), mgp = c(3.5, 0.7, 0) )
    plot(reorder[1:topN],type='h', xaxt='n', cex.axis=1.5, ylab="", bty="n",lwd=3,col="black",xlab="",las=0, ylim=c(-0.6,0.6)) #Coefficients ordered by effect size Combined model coefficient , ylim=c(2*min(reorder),max(reorder)*2)
    text(1:topN,ifelse(reorder[1:topN] > 0,reorder[1:topN]+max(reorder)*0.05,reorder[1:topN]-max(reorder)*0.05),names(reorder[1:topN]),offset=0,adj=0,srt=90,pos=ifelse(reorder[1:topN] > 0,4,2),lwd=3,col=cols,cex=0.5)
    abline(h=0,lty=2)
    dev.off()

    topcoefs
}))

pdf("png/BC3MSD.Lasso_mouse.pdf") #FigS4b.pdf
par(oma=c(4, 4, 1, 1), mar = c(5,5,4,4))
allcorrs = do.call(rbind, lapply( names, function(sample) {
  say(sample)
	a=read.delim(paste("preds/",sample,".Lasso.yhat.tsv",sep=''),F)
  colnames(a) <- c("y","yhat")
  y=a$y
  yhat=a$yhat
  heatscatter(yhat, y, xlab="CV fitted model prediction", ylab="Half life", bty='n', cex=0.3, las=1, xlim=c(-3,2), ylim=c(-4,4))

  legend(-3,4, cex=1.7, bg="white", bty="n", legend = c(sprintf("   r = %.2f", cor(yhat, y,method="pearson")),
                                                          sprintf("rho = %.2f", cor(yhat, y,method="spearman"))))
}))
dev.off()

topcoefs

#half lives, kmer freqs, & codon freqs features
a=readme("seqFeatWithKmerFreqs.txt.gz",gzip=T)
a=a[,c(1:10,grep("Codon",colnames(a)),grep("3UTR",colnames(a)))]
colnames(a)[3:10]=paste0("Basic.", colnames(a)[3:10])
head(colnames(a),40)

# miRNA prediction
if(grepl("M",names)){
  b=readme("CWCS.txt.gz",gzip=T)
  colnames(b)[2:ncol(b)]=paste0("MIR.",colnames(b)[2:ncol(b)])
  b[2:ncol(b)] = -1*b[2:ncol(b)]
  a=merge(a,b,by=1,all.x=1)
  a[is.na(a)]=0
  say("with miRNA features")
}

seqpreds <- function(tool){
  humanorf=readme(paste0(tool,"_predictions_mm10_90/ORF_avg.txt.gz"), gzip=T)
  human3p=readme(paste0(tool,"_predictions_mm10_90/3pUTR_avg.txt.gz"), gzip=T)
  human5p=readme(paste0(tool,"_predictions_mm10_90/5pUTR_avg.txt.gz"), gzip=T)

  colnames(humanorf)[2:ncol(humanorf)]=paste0(tool,".",colnames(humanorf)[2:ncol(humanorf)],".ORF")
  colnames(human3p)[2:ncol(human3p)]=paste0(tool,".",colnames(human3p)[2:ncol(human3p)],".3UTR")
  colnames(human5p)[2:ncol(human5p)]=paste0(tool,".",colnames(human5p)[2:ncol(human5p)],".5UTR")
  a=merge(humanorf,human3p,by=1,all=T)
  a=merge(a,human5p,by=1,all=T)
  a[is.na(a)]=0
  a
}

if(grepl("S",names)){
  say("with SeqWeaver features")
  a=merge(a,seqpreds("SeqWeaver"),by=1)
}
if(grepl("D",names)){
  say("with DeepRiPe features")
  a=merge(a,seqpreds("DeepRiPe"),by=1)
}

data <- a
data = data[,3:ncol(data)]
head(colnames(data))

cormat=cor(data)
allcorrcoefs=unique(unlist(lapply(topcoefs, function(x){
  y=cormat[rownames(cormat)==x,]
  names(y[abs(y)>0.8])
})))

heatmaps=function(x) {
  typerow = as.factor(unlist(lapply(rownames(x), function(x) strsplit(x, '\\.')[[1]][1])))
  typecol = as.factor(unlist(lapply(colnames(x), function(x) strsplit(x, '\\.')[[1]][1])))
  say(levels(typerow))
  say(levels(typecol))
  colnames(x) = unlist(lapply(colnames(x), function(x) paste0(strsplit(x, '\\.')[[1]][-1],collapse='.')))
  rownames(x) = unlist(lapply(rownames(x), function(x) paste0(strsplit(x, '\\.')[[1]][-1],collapse='.')))
  colsrow = c("darkorange", "red", "purple", "green", "cyan", "blue")[typerow]
  colscol = c("darkorange", "red", "purple", "green", "cyan", "blue")[typecol]
  gplots::heatmap.2(as.matrix(x),density.info="none",
    trace="none", breaks=seq(-1,1,0.05), symkey=FALSE, cexRow=1,
    cexCol=1, key=TRUE, colRow=colsrow, colCol=colscol,
    col=colorRampPalette(c("blue","white","red"))(40)
  )
}

pdf("png/Top30Coefs_80Pearson_mouse.pdf", width=18,height=7) #FigS4d.pdf
par(mar=c(10,4,3,5), oma=c(15,4,3,5))
heatmaps(cormat[rownames(cormat) %in% topcoefs,colnames(cormat) %in% allcorrcoefs])
dev.off()
