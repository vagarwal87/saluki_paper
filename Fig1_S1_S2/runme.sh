#prepare human and mouse tables (Table S1)
Rscript master_comparison_human.R
Rscript master_comparison_mouse.R

gzip all_HLs_human.txt
gzip all_HLs_mouse.txt

#visualize Spearman correlations & bar plots
Rscript Fig1_S1.R
Rscript FigS2.R
