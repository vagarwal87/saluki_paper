library(ggplot2)
library(glmnet)
library(LSD)

a=read.csv("experimentResults.csv")
round(a,3)

results = data.frame(tests = c("B vs BG","B vs Bb","Bb vs BbG"), pvals =
c(t.test(a$B, a$BG, paired=T,alternative='less')$p.value,
t.test(a$B, a$Bb, paired=T,alternative='less')$p.value,
t.test(a$Bb, a$BbG, paired=T,alternative='less')$p.value
))
results$p.corrected = pmin(results$pvals*nrow(results),1)
results

meltData <- reshape::melt(a)

pdf("png/violinplots_jointmodel.pdf",width=8,height=4) #FigS5a.pdf
ggplot(meltData, aes(x=variable, y=value)) + geom_violin()
p<-ggplot(meltData, aes(x=variable, y=value)) + geom_violin(position=position_dodge(1))
print(p)
dev.off()

names = c("BGb")

topcoefs = unlist(lapply( c(names), function(sample) {
    topN=30
    load(paste("Robj/", sample,"_CV-Lasso.Robj",sep=''))
    coef.exact = coef(cvfit, s = "lambda.min", exact = TRUE)
    coefficients <- as.numeric(coef.exact)
    names(coefficients) <- coef.exact@Dimnames[[1]]
    select <- coefficients[coefficients!=0]
    select[order(abs(select),decreasing=T)]
    reorder <- select[order(abs(select),decreasing=T)]
    names(reorder)=gsub("SeqWeaver","SeqWeaver.",names(reorder))

    topcoefs=names(reorder[1:topN])
    type = as.factor(unlist(lapply(names(reorder), function(x) strsplit(x, '\\.')[[1]][1])))
    names(reorder) = unlist(lapply(names(reorder), function(x) paste0(strsplit(x, '\\.')[[1]][-1],collapse='.')))
    cols = c("darkorange", "red", "blue", "darkorange")[type]

    pdf(paste("png/", sample,"_CV-Lasso-coefficients.pdf",sep=''),width=8,height=8) #FigS5c.pdf
    par(mar = c(5, 6, 4, 2) + 0.1, oma= c(0.2, 0, 0, 0), mgp = c(3.5, 0.7, 0) )
    cols2 <- ifelse(reorder > 0,"darkblue","darkred")
    plot(reorder[1:topN],type='h', xaxt='n', cex.axis=1.5, ylab="", bty="n",lwd=3,col="black",xlab="",las=0, ylim=c(-0.6,0.6)) #Coefficients ordered by effect size Combined model coefficient , ylim=c(2*min(reorder),max(reorder)*2)
    text(1:topN,ifelse(reorder[1:topN] > 0,reorder[1:topN]+max(reorder)*0.05,reorder[1:topN]-max(reorder)*0.05),names(reorder[1:topN]),offset=0,adj=0,srt=90,pos=ifelse(reorder[1:topN] > 0,4,2),lwd=3,col=cols,cex=0.5)
    abline(h=0,lty=2)
    dev.off()

    topcoefs
}))

pdf("png/BGb.Lasso.pdf") #FigS5b.pdf
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

load("alldataGeneticBiochemical.RData")

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
  colsrow = c("darkorange", "red", "blue", "darkorange")[typerow]
  colscol = c("darkorange", "red", "blue", "darkorange")[typecol]
  gplots::heatmap.2(as.matrix(x),density.info="none",
    trace="none", breaks=seq(-1,1,0.05), symkey=FALSE, cexRow=1,
    cexCol=1, key=TRUE, colRow=colsrow, colCol=colscol,
    col=colorRampPalette(c("blue","white","red"))(40)
  )
}

pdf("png/Top30Coefs_80Pearson_jointmodel.pdf", width=14,height=7) #FigS5d.pdf
par(mar=c(10,4,3,5), oma=c(2,4,3,5))
heatmaps(cormat[rownames(cormat) %in% topcoefs,colnames(cormat) %in% allcorrcoefs])
dev.off()
