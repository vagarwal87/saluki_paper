#result stored in datasets/
./generate_training_input.pl all_HLs_mouse_PC1.txt | gzip -c >all_HLs_mouse_featTable.txt.gz

Rscript ../calc_kmer_freqs.R #replace input file with that of mouse
#manually added 1st column name to seqFeatWithKmerFreqs.txt
gzip seqFeatWithKmerFreqs.txt
#result stored in datasets/

# Benchmark the following models:
# B: basic general regional features - 8 features
# BC: basic+codon - 70 features
# B5: basic+kmers(in 5'UTR regions) - 21852 features
# B5C: basic+kmers(in 5'UTR regions)+codon
# BO: basic+kmers(in ORF regions) - 21852 features
# B3: basic+kmers(in 3'UTR regions) - 21852 features
# BCO: basic+codon+kmers(in ORF regions) - 21914 features
# B3M: basic+kmers(in 3'UTR regions)+miRNA - 22171 features
# BC3M: basic+codon+kmers(in 3'UTR regions)+miRNA - 22233 features
# BC3MS: basic+codon+kmers(in 3'UTR regions)+seqweaver - 23013 features
# BC3MD: basic+codon+kmers(in 3'UTR regions)+deepripe - 22410 features
# BC3MSD: basic+codon+kmers(in 3'UTR regions)+miRNA+seqweaver+deepripe - 23190 features

#request 90Gb per job
for x in B BC B5 B5C BO B3 BCO B3M BC3M BC3MS BC3MD BC3MSD; do { sbatch --mem 90000 -J $x -o $x.txt --wrap="Rscript lasso_model.R $x"; } done

#after running the same models for human
Rscript cross_species_performance.R

Rscript FigS4abcd.R
