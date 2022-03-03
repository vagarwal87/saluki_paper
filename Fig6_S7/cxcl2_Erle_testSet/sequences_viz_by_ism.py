#!/usr/bin/env python

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from basenji import dna_io
import logomaker

def plot_logo(seq_1hot, scores, viz_len=None, ymin=-np.inf, ymax=np.inf, outfile=None):
    if viz_len is not None and viz_len < seq_1hot.shape[0]:
        mid_pos = seq_1hot.shape[0]//2
        viz_start = mid_pos - viz_len//2
        viz_end = viz_start + viz_len
        scores = scores[viz_start:viz_end]
        seq_1hot = seq_1hot[viz_start:viz_end]
    scores = np.clip(scores, ymin, ymax)

    fig, axs = plt.subplots(3, 1, figsize=(25,5))

    tmpscores = np.mean(scores,1)
    tmpscores = np.repeat(tmpscores[:,np.newaxis], 4, axis=1)
    scores2 = tmpscores - scores

    seq_df = pd.DataFrame(seq_1hot*scores2, columns=['A','C','G','T'])
    print(seq_df)
    seq_logo = logomaker.Logo(seq_df, ax=axs[0])
    if not np.isinf(ymin):
        seq_logo.ax.set_ylim(bottom=ymin)
    if not np.isinf(ymax):
        seq_logo.ax.set_ylim(top=ymax)

    sns.heatmap(scores.T, center=0, cbar=False, ax=axs[1], cmap = "RdBu_r")
    axs[1].set_yticklabels('ACGT')
    axs[1].tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are of
    sns.heatmap(scores.T, center=0, cbar=True, ax=axs[2], cmap = "RdBu_r")

    if outfile is not None:
        fig.savefig(outfile,format="pdf")
        plt.close()

mpra = pd.read_table("mutagenesis_MPRA.txt", header=0)
mpra = np.array(mpra.iloc[:,1:5])
mpra = np.flip(mpra,0)

preds = pd.read_table("mutagenesis_MPRA_predictions.txt", header=0)
preds = np.array(preds.iloc[:,1:5])
seq = preds==0
seq = seq.astype("float32")
seq = np.flip(seq,0)
preds = np.flip(preds,0)

print(dna_io.hot1_dna(seq))

plot_logo(seq, mpra, outfile='mutagenesis_mpra.pdf', ymin=-0.25, ymax=0.5)
plot_logo(seq, preds, outfile='mutagenesis_mpra_predictions.pdf', ymin=-0.1, ymax=0.25)
