Rscript prepare_features.R
#output stored in datasets/

# Benchmark the following models:
# B: basic general regional features - 8 features
# BE: basic+codon
# BP: basic+kmers(in 5'UTR regions)
# BR: basic+kmers(in ORF regions)
# BM: basic+kmers(in 3'UTR regions)
# Be: basic+codon+kmers(in ORF regions)
# BER: basic+kmers(in 3'UTR regions)+miRNA
# BEM: basic+codon+kmers(in 3'UTR regions)+miRNA
# BeEM: basic+codon+kmers(in 3'UTR regions)+miRNA

#request 20Gb per job
for x in B BE BP BR BM Be BEM BEe BER BEeM; do { sbatch --mem 20000 -J $x -o $x.txt --wrap="Rscript lasso_model.R $x"; } done

Rscript Fig4.R
