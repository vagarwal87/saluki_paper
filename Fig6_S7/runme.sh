#plot benchmarking results of neural net vs lasso regression for human/mouse

#Zhao et al 2014 dataset
cd cxcl2_Erle_testSet

#prepare table of variant effects
Rscript joint_table.R
#prepare fastUTR data
Rscript process_seqs_and_mpraVals.R

#FastUTR predictions from all models
for x in {0..9}; do { for y in {0..4}; do { echo $x, $y; sbatch --mem 20000 -J B$x$y -e fastUTR_ENCFF090JTW_predictions/$x\_c$y.err \
  -o fastUTR_ENCFF090JTW_predictions/f$x\_c$y.txt --wrap=". /home/drk/anaconda3/etc/profile.d/conda.sh; conda activate tf2.6-rna; python fastUTR_predict.py train_gru/params.json train_gru/f$x\_c$y/train/model0_best.h5 ENCFF090JTW.sequences.txt"; } done } done

#Mutagenesis predictions from all models
for x in {0..9}; do { for y in {0..4}; do { echo $x, $y; sbatch --mem 20000 -J $x$y -e f$x\_c$y.err -o f$x\_c$y.txt \
  --wrap=". /home/drk/anaconda3/etc/profile.d/conda.sh; conda activate tf2.6-rna; \
  python mutagenesis_predict.py train_gru/params.json train_gru/f$x\_c$y/train/model0_best.h5"; } done } done

#plot of variant effects for CXCL2 saturation mutagenesis experiment
Rscript Fig6c.R

#plot of 3' UTR effects for fastUTR experiment
Rscript Fig6d.R

#For data shown in Fig6a
python sequences_viz_by_ism_fullcxcl2.py

#For plot shown in Fig6b
python sequences_viz_by_ism.py
