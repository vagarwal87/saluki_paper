library(ggplot2)
library(glmnet)
library(LSD)

a=read.csv("experimentResults.csv")
round(a,3)

results = data.frame(tests = c("BC3MS vs Saluki human","BeEM vs Saluki human","BC3MSD vs Saluki mouse"), pvals =
c(t.test(a[,1], a[,3], paired=T)$p.value,
t.test(a[,2], a[,3], paired=T)$p.value,
t.test(a[,4], a[,5], paired=T)$p.value))
results

meltData <- reshape::melt(a)

pdf("png/violinplots.pdf",width=4,height=4) #Fig5b
p<-ggplot(meltData, aes(x=variable, y=value)) + geom_violin(position=position_dodge(1))
print(p)
dev.off()
