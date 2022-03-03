library(rhdf5)
library(ggplot2)
library(abind)

files = args

getdata = function(track){
  do.call(abind, lapply(files, function(x){
    h5read(x,track)
  }))
}

values = getdata("ism") #(1, 4, 12288, nseqs)
coding = getdata("coding")  #(12288, nseqs)
seqs = getdata("seqs")  #(12288, nseqs)
preds = getdata("ref")  #(12288, nseqs)

NBINS = 100

binnedval = function(x){
  a = data.frame(vals=x, bin = 1:length(x))
  a$bin = as.numeric(cut_number(a$bin, NBINS))
  a = aggregate(a$vals,list(a$bin),function(x) mean(x, na.rm=T))
  b=rep(NA,NBINS)
  b[a[,1]]=a[,2]
  b
}

allutr5=c()
allorf=c()
allutr3=c()
allpreds=c()

for (i in 1:dim(values)[4]){ # iterate through all sequences
  if(i%%1000 == 0) say(i)

  X=abs(values[1,,,i])*(seqs[,,i]=="TRUE")
  myvals = apply(X,2,sum)

  txstart = which(myvals!=0)[1]
  txend = 12288
  orf = which(coding[,i] == "TRUE")
  orfstart = orf[1]
  orfend = min(orf[length(orf)]+3, txend)

  utr5 = myvals[txstart:(orfstart-1)]
  orf = myvals[orfstart:orfend]
  utr3 = myvals[orfend:txend]

  if(length(utr3) < 100 || length(orf) < 100 || length(utr5) < 100) next

  if(i==1){
    allutr5=binnedval(utr5)
    allorf=binnedval(orf)
    allutr3=binnedval(utr3)
  }
  else{
    allutr5=rbind(allutr5,binnedval(utr5))
    allorf=rbind(allorf,binnedval(orf))
    allutr3=rbind(allutr3,binnedval(utr3))
  }
  allpreds = c(allpreds, preds[1,i])
}

master_table = as.data.frame(cbind(allutr5, allorf, allutr3))
master_table$bin <- cut(allpreds, breaks=quantile(allpreds), include.lowest = TRUE) #4 prediction bins
save(master_table, file="master_table.Rdata")

load("master_table.Rdata")
table(master_table$bin)
dim(master_table)
xin = 1:(NBINS*3)

meanvals = lapply(levels(master_table$bin), function(x) apply(master_table[master_table$bin==x,xin],2,function(x) mean(x, na.rm=T)))

pdf("png/Fig5E.pdf")
plot(xin, predict(loess(meanvals[[4]]~xin, span=0.05)), col='cyan', type='l') #, ylim=c(-0.004,0.004), axes = FALSE
lines(xin, predict(loess(meanvals[[3]]~xin, span=0.05)), col='blue', type='l')
lines(xin, predict(loess(meanvals[[2]]~xin, span=0.05)), col='red', type='l')
lines(xin, predict(loess(meanvals[[1]]~xin, span=0.05)), col='black', type='l')
abline(v=NBINS)
abline(v=NBINS*2)
plot(xin, predict(loess(meanvals[[4]]~xin, span=0.02)), col='cyan', type='l') #, ylim=c(-0.004,0.004), axes = FALSE
lines(xin, predict(loess(meanvals[[3]]~xin, span=0.02)), col='blue', type='l')
lines(xin, predict(loess(meanvals[[2]]~xin, span=0.02)), col='red', type='l')
lines(xin, predict(loess(meanvals[[1]]~xin, span=0.02)), col='black', type='l')
abline(v=NBINS)
abline(v=NBINS*2)
 plot(xin, predict(loess(meanvals[[4]]~xin, span=0.01)), col='cyan', type='l') #, ylim=c(-0.004,0.004), axes = FALSE
lines(xin, predict(loess(meanvals[[3]]~xin, span=0.01)), col='blue', type='l')
lines(xin, predict(loess(meanvals[[2]]~xin, span=0.01)), col='red', type='l')
lines(xin, predict(loess(meanvals[[1]]~xin, span=0.01)), col='black', type='l')
abline(v=NBINS)
abline(v=NBINS*2)
 plot(xin, meanvals[[4]], col='cyan', type='l')
lines(xin, meanvals[[3]], col='blue', type='l')
lines(xin, meanvals[[2]], col='red', type='l')
lines(xin, meanvals[[1]], col='black', type='l')
abline(v=NBINS)
abline(v=NBINS*2)
dev.off()
