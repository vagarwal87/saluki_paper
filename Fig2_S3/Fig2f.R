library(LSD)

a=read.delim("all_HLs_mouse_PC1.txt")
b=read.delim("all_HLs_human_PC1.txt")
c=read.delim("1to1_orthologs_expression.txt",F)
c=c[,1:2]
colnames(c)=c("humangene","mousegene")

c=merge(c,b,by=1)
c=merge(c,a,by.x=2,by.y=1)

"Sample size:"
nrow(c)
"Pearson"
cor(c$human_PC1, c$mouse_PC1)
"Spearman"
cor(c$human_PC1, c$mouse_PC1, method='spearman')

pdf("Fig2f.pdf")
heatscatter(c$human_PC1, c$mouse_PC1, xlab="Human half-life, PC1", ylab="Mouse half-life, PC1", bty='n', cex=0.3, las=1)
dev.off()
