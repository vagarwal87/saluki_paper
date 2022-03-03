library(glmnet)

database = args[1] #name of database & subset of features to evaluate

#half lives, kmer freqs, & codon freqs features
a=readme("seqFeatWithKmerFreqs.txt.gz",gzip=T)

if(database == "B") a=a[,1:10]
if(database == "BC") a=a[,c(1:10,grep("Codon",colnames(a)))]
if(database == "B5") a=a[,c(1:10,grep("5UTR",colnames(a)))]
if(database == "BO") a=a[,c(1:9,grep("ORF",colnames(a)))]
if(database == "B3" | database == "B3M") a=a[,c(1:10,grep("3UTR",colnames(a)))]
if(database == "BCO") a=a[,c(1:9,grep("Codon",colnames(a)),grep("ORF",colnames(a)))]
if(database == "B5C") a=a[,c(1:10,grep("Codon",colnames(a)),grep("5UTR",colnames(a)))]
if(grepl("BC3",database)) a=a[,c(1:10,grep("Codon",colnames(a)),grep("3UTR",colnames(a)))]

# miRNA prediction
if(grepl("M",database)){
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

  colnames(humanorf)[2:ncol(humanorf)]=paste0(tool,colnames(humanorf)[2:ncol(humanorf)],".ORF")
  colnames(human3p)[2:ncol(human3p)]=paste0(tool,colnames(human3p)[2:ncol(human3p)],".3UTR")
  colnames(human5p)[2:ncol(human5p)]=paste0(tool,colnames(human5p)[2:ncol(human5p)],".5UTR")
  a=merge(humanorf,human3p,by=1,all=T)
  a=merge(a,human5p,by=1,all=T)
  a[is.na(a)]=0
  a
}

if(grepl("S",database)){
  say("with SeqWeaver features")
  a=merge(a,seqpreds("SeqWeaver"),by=1)
}
if(grepl("D",database)){
  say("with DeepRiPe features")
  a=merge(a,seqpreds("DeepRiPe"),by=1)
}

say("With Ensembl IDs + Features: ", nrow(a))

say("Dimensions: ", dim(a))
a[1:10,1:5]
genes = a[,1]
a[,1]=NULL

x=scale(a[,2:ncol(a)])
y=scale(a[,1])

cvfit = cv.glmnet(x,y,nlambda=dim(x)[2],standardize=F,lambda=10^(seq(-10,-1,0.1)))
png(paste("png/", database, "_CV-Lasso.png",sep=''),width=1000,height=600,type="cairo")
plot(cvfit,main="meanVal")
dev.off()

cor(y,predict(cvfit,newx=x,type="response",s="lambda.min"),method="spearman")^2
cor(y,predict(cvfit,newx=x,type="response",s="lambda.min"),method="pearson")^2

save(cvfit,file=paste("Robj/", database,"_CV-Lasso.Robj",sep=''))

coef.exact = coef(cvfit, s = "lambda.min", exact = TRUE)
coefficients <- as.numeric(coef.exact)
names(coefficients) <- coef.exact@Dimnames[[1]]
select <- coefficients[coefficients!=0]
select[order(abs(select),decreasing=T)]
length(select)

srho <- cor(y,predict(cvfit,newx=x,type="response",s="lambda.min"),method="spearman")
prho <- cor(y,predict(cvfit,newx=x,type="response",s="lambda.min"),method="pearson")

# ####ANALYSIS WITH TRAINED MODEL
load(paste("Robj/", database,"_CV-Lasso.Robj",sep=''))

nrObs <- dim(x)[1]
bSize <- floor(nrObs/10)
set.seed(42)
bins <- data.frame(bin=c(sapply(1:9,function(x)rep(x,bSize)),rep(10,nrObs-9*bSize)),rows=sample(1:nrObs))
bins <- bins[order(bins$bin,bins$rows),]

binnedgenes = cbind(bins[order(bins$rows),"bin"], genes, y)
colnames(binnedgenes) = c("BIN","GENEID","HALFLIFE")
write.table(binnedgenes,file="binnedgenes.txt",sep="\t",quote=F,row.names=F,col.names=T)

yhat <- rep(NA,length(y))

for (i in 1:10) {
  selRows <- bins$rows[bins$bin != i]
  predRows <- bins$rows[bins$bin == i]
  bx <- x[selRows,]
  by <- y[selRows]
  bcvfit = cv.glmnet(bx,by,nlambda=dim(x)[2],standardize=F,lambda=10^(seq(-10,-1,0.1)))
  yhat[predRows] = predict(bcvfit,newx=x[predRows,],type="response",s=cvfit$lambda.min)
  say(cor(yhat[predRows], y[predRows]))
}

res <- data.frame(y,yhat,stringsAsFactors=F)
write.table(res,file=paste("preds/", database,".Lasso.yhat.tsv",sep=''),sep="\t",quote=F,row.names=F,col.names=F)

cor(y,yhat,method="spearman")
cor(y,yhat,method="pearson")

###############################################
## Draw unbiased scatter plots
###############################################

res <- read.table(paste("preds/", database,".Lasso.yhat.tsv",sep=''),sep="\t",as.is=T)
colnames(res) <- c("y","yhat")

y=res$y
yhat=res$yhat

colors <- c("#377eb8","#4daf4a","#984ea3","#ff7f00") # "gray18","#e41a1c",
png(paste("png/", database,"_scatter.Lasso.png",sep=''),width=600,height=600,type="cairo")
plot(y,yhat, main=sprintf("Half life vs combined model\n(p-rho %.2f, s-rho %.2f)", cor(y,yhat,method="pearson"), cor(y,yhat,method="spearman")), xlab="Half-life", ylab="CV fitted model prediction", pch=19, cex=0.5, col = colors[1], las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
abline(0,1,lty=2)
dev.off()
