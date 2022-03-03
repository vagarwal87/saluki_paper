library(glmnet)

database = args[1] #name of database {choose a code below}

load("alldata.RData")

if(database == "B") a=a[,1:10]
if(database == "BE") a=a[,c(1:10,grep("^ENCORI",colnames(a)))]
if(database == "BP") a=a[,c(1:10,grep("^PARCLIP",colnames(a)))]
if(database == "BR") a=a[,c(1:10,grep("^RIP|^TE",colnames(a)))]
if(database == "BM") a=a[,c(1:10,grep("^M6ACLIP",colnames(a)))]
if(database == "Be") a=a[,c(1:10,grep("^eCLIP",colnames(a)))]
if(database == "BER") a=a[,c(1:10,grep("^ENCORI|^RIP|^TE",colnames(a)))]
if(database == "BEM") a=a[,c(1:10,grep("^ENCORI|^M6ACLIP",colnames(a)))]
if(database == "BEe") a=a[,c(1:10,grep("^ENCORI|^eCLIP",colnames(a)))]
if(database == "BEeM") a=a[,c(1:10,grep("^ENCORI|^M6ACLIP|^eCLIP",colnames(a)))]

colnames(a)

say("With Ensembl IDs + Features: ", nrow(a))

say("Dimensions: ", dim(a))
a[1:10,1:5]
genes = a[,1]
a[,1]=NULL

x=scale(a[,2:ncol(a)])
y=scale(a[,1])

# summary(lm(Z[,1]~SEQ.UTR5LEN+SEQ.CDSLEN+SEQ.INTRONLEN+SEQ.UTR3LEN+SEQ.UTR5GC+SEQ.CDSGC+SEQ.UTR3GC+SEQ.ORFEXONDENSITY, data=Z))

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

png(paste("png/", database,"_CV-Lasso-coefficients.png",sep=''),width=2000,height=600,type="cairo")
par(mar = c(5, 6, 4, 2) + 0.1, oma= c(0.2, 0, 0, 0), mgp = c(3.5, 0.7, 0) )
reorder <- select[order(abs(select),decreasing=T)]
cols <- ifelse(reorder > 0,"darkblue","darkred")
cols["(Intercept)"] <- "darkgray"
plot(reorder,type='h',ylab="Combined model coefficient",ylim=c(1.7*min(reorder),max(reorder)*1.4),lwd=3,col=cols,xlab="Coefficients ordered by effect size",main="meanVal Lasso model starting with all annotations and individual LS-GKM scores",las=1) #sprintf("(p-rho %.4f / s-rho = %.4f / s-rho² = %.2f)",prho,srho,srho^2)
text(1:length(select),ifelse(reorder > 0,reorder+max(reorder)*0.05,reorder-max(reorder)*0.05),names(reorder),offset=0,adj=0,cex=1.1,srt=90,pos=ifelse(reorder > 0,4,2),lwd=3,col=cols)
abline(h=0,lty=2)
dev.off()

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

cor(y,yhat,method="spearman")^2
cor(y,yhat,method="pearson")^2

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
