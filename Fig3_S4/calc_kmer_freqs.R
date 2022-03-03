library(data.table)
library(Biostrings)
library(reshape2)

kmerlen = 7

file1 = "all_HLs_human_featTable.txt.gz" #replace with "all_HLs_mouse_featTable.txt.gz" for mouse

inp.tbl <- fread(paste("zcat", file1),header=T,data.table=F,sep="\t")
rownames(inp.tbl) = inp.tbl[,1]
inp.tbl[,1] = NULL

#get codon frequencies
inp.tbl=cbind(inp.tbl, do.call(rbind, lapply(inp.tbl$ORF, function(x){
    y = oligonucleotideFrequency(DNAStringSet(x), 3, step=3)
    colnames(y)=paste0("Codon.",colnames(y))
    y/sum(y)
    # y
})))

#get kmer frequencies in ORF
for(k in 1:kmerlen){
  inp.tbl=cbind(inp.tbl, do.call(rbind, lapply(inp.tbl$ORF, function(x){
      y = oligonucleotideFrequency(DNAStringSet(x), k, step=1)
      colnames(y)=paste0("ORF.",colnames(y))
      y/sum(y)
      # y
  })))
}

#get kmer frequencies in 5UTR
for(k in 1:kmerlen){
  inp.tbl=cbind(inp.tbl, do.call(rbind, lapply(inp.tbl$"5UTR", function(x){
      y = oligonucleotideFrequency(DNAStringSet(x), k, step=1)
      colnames(y)=paste0("5UTR.",colnames(y))
      y/sum(y)
      # y
  })))
}

#get kmer frequencies in 3UTR
for(k in 1:kmerlen){
  inp.tbl=cbind(inp.tbl, do.call(rbind, lapply(inp.tbl$"3UTR", function(x){
      y = oligonucleotideFrequency(DNAStringSet(x), k, step=1)
      colnames(y)=paste0("3UTR.",colnames(y))
      y/sum(y)
      # y
  })))
}

inp.tbl$ORF = NULL
inp.tbl$"5UTR" = NULL
inp.tbl$"3UTR" = NULL

inp.tbl=inp.tbl[,apply(inp.tbl,2,var)!=0]
inp.tbl[,c(2:5, 9)] = log10(inp.tbl[,c(2:5, 9)]+0.1)

writefile(inp.tbl,"seqFeatWithKmerFreqs.txt", row.names=T)
