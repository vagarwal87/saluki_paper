#plot benchmarking results of neural net vs lasso regression for human/mouse
Rscript Fig5bc.R
Rscript Fig5c_S6c.R

#metagene plot of ISM analysis
Rscript Fig5d.R

#visualize a few examples of ISM results
python Fig5e.py

#ablation analysis
Rscript CV_predictions_ablations.R
Rscript FigS6b.R

#plot neural network model structure
python FigS6a.py

#insertional analysis of motifs detected by TF-MODISCO
for x in {0..9}; do { for y in {0..4}; do { REG="5utr"; echo $x, $y; sbatch --mem 20000 -J $REG$x$y -e spikeinmotif/$REG\_f$x\_c$y.err -o spikeinmotif/$REG\_f$x\_c$y.txt --wrap=". /home/drk/anaconda3/etc/profile.d/conda.sh; conda activate tf2.6-rna; python spikeinmotif.py --region $REG --seqs ISM/f$x\_c$y\_scores.h5 train_gru/params.json train_gru/f$x\_c$y/train/model0_best.h5"; } done } done
for x in {0..9}; do { for y in {0..4}; do { REG="orf"; echo $x, $y; sbatch --mem 20000 -J $REG$x$y -e spikeinmotif/$REG\_f$x\_c$y.err -o spikeinmotif/$REG\_f$x\_c$y.txt --wrap=". /home/drk/anaconda3/etc/profile.d/conda.sh; conda activate tf2.6-rna; python spikeinmotif.py --region $REG --seqs ISM/f$x\_c$y\_scores.h5 train_gru/params.json train_gru/f$x\_c$y/train/model0_best.h5"; } done } done
for x in {0..9}; do { for y in {0..4}; do { REG="3utr"; echo $x, $y; sbatch --mem 20000 -J $REG$x$y -e spikeinmotif/$REG\_f$x\_c$y.err -o spikeinmotif/$REG\_f$x\_c$y.txt --wrap=". /home/drk/anaconda3/etc/profile.d/conda.sh; conda activate tf2.6-rna; python spikeinmotif.py --region $REG --seqs ISM/f$x\_c$y\_scores.h5 train_gru/params.json train_gru/f$x\_c$y/train/model0_best.h5"; } done } done

#codon insertional analysis
for x in {0..9}; do { for y in {0..4}; do { REG="orf"; echo $x, $y; sbatch --mem 20000 -J codon0$REG$x$y -e spikeinmotif/codon_frame0_$REG\_f$x\_c$y.err -o spikeinmotif/codon_frame0_$REG\_f$x\_c$y.txt --wrap=". /home/drk/anaconda3/etc/profile.d/conda.sh; conda activate tf2.6-rna; python spikeinmotif.py --codons --frame 0 --region $REG --seqs ISM/f$x\_c$y\_scores.h5 train_gru/params.json train_gru/f$x\_c$y/train/model0_best.h5"; } done } done

#splice site insertional analysis
for x in {0..9}; do { for y in {0..4}; do { REG="5utr"; echo $x, $y; sbatch --mem 20000 -J S$REG$x$y -e spikeinmotif/splice$REG\_f$x\_c$y.err -o spikeinmotif/splice$REG\_f$x\_c$y.txt --wrap=". /home/drk/anaconda3/etc/profile.d/conda.sh; conda activate tf2.6-rna; python spikeinmotif.py --addsplicesites --region $REG --seqs ISM/f$x\_c$y\_scores.h5 train_gru/params.json train_gru/f$x\_c$y/train/model0_best.h5"; } done } done
for x in {0..9}; do { for y in {0..4}; do { REG="orf"; echo $x, $y; sbatch --mem 20000 -J S$REG$x$y -e spikeinmotif/splice$REG\_f$x\_c$y.err -o spikeinmotif/splice$REG\_f$x\_c$y.txt --wrap=". /home/drk/anaconda3/etc/profile.d/conda.sh; conda activate tf2.6-rna; python spikeinmotif.py --addsplicesites --region $REG --seqs ISM/f$x\_c$y\_scores.h5 train_gru/params.json train_gru/f$x\_c$y/train/model0_best.h5"; } done } done
for x in {0..9}; do { for y in {0..4}; do { REG="3utr"; echo $x, $y; sbatch --mem 20000 -J S$REG$x$y -e spikeinmotif/splice$REG\_f$x\_c$y.err -o spikeinmotif/splice$REG\_f$x\_c$y.txt --wrap=". /home/drk/anaconda3/etc/profile.d/conda.sh; conda activate tf2.6-rna; python spikeinmotif.py --addsplicesites --region $REG --seqs ISM/f$x\_c$y\_scores.h5 train_gru/params.json train_gru/f$x\_c$y/train/model0_best.h5"; } done } done

#plot results of insertional analysis
Rscript Fig5fgh.R

#For FigS6d
#changed code to run for orf, 5utr, and 3utr independently
bash run_modisco.sh
#modisco results, mapped to known motifs via tomtom, saved in *_scores/out/tomtom.html
