library(ggplot2)
library(glmnet)
library(LSD)

a=read.csv("experimentResults_ablations.csv")
round(a,3)

results = data.frame(tests = c("S vs SCs human","SC vs SCs human","Ss vs SCs human",
"S vs SCs mouse","SC vs SCs mouse","Ss vs SCs mouse"), pvals =
c(t.test(a[,1], a[,4], paired=T,alternative='less')$p.value,
t.test(a[,2], a[,4], paired=T,alternative='less')$p.value,
t.test(a[,3], a[,4], paired=T,alternative='less')$p.value,
t.test(a[,5], a[,8], paired=T,alternative='less')$p.value,
t.test(a[,6], a[,8], paired=T,alternative='less')$p.value,
t.test(a[,7], a[,8], paired=T,alternative='less')$p.value))

results$p.corrected = pmin(results$pvals*nrow(results),1)
results

meltData <- reshape::melt(a)

pdf("png/violinplots_ablations.pdf",width=4,height=4) #FigS6b
p<-ggplot(meltData, aes(x=variable, y=value)) + geom_violin(position=position_dodge(1))
print(p)
dev.off()
