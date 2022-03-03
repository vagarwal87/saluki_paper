library("readxl")

a = read_excel("41587_2014_BFnbt2851_MOESM12_ESM.xls")
# xlsx files
b = read_excel("41587_2014_BFnbt2851_MOESM13_ESM.xls",skip=3)

head(a)
head(b)

c=merge(a,b,by=1)
writefile(c[c[,8]!="nd",c(7,8)], "fastUTR_mpra.txt", col.names=F)
