#!/usr/bin/env python

import h5py
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from basenji.plots import *
from basenji.bin.basenji_sat_plot import *
from basenji import dna_io

# take the scores
h5 = h5py.File('ISM/f9_scores.h5','r')
# h5 = h5py.File('ISM/f0_scores.h5','r')
seqs = h5['seqs'][:]
scores = h5['ism'][:] #scores
h5.close()

cov_index = 0
target_scores = scores[:,:,:,cov_index]

def plot_heat(ax, sat_delta_ti, min_limit, cbar=False):
  """ Plot satmut deltas.

    Args:
        ax (Axis): matplotlib axis to plot to.
        sat_delta_ti (4 x L_sm array): Single target delta matrix for saturated mutagenesis region,
        min_limit (float): Minimum heatmap limit.
    """

  vlim = max(min_limit, np.nanmax(np.abs(sat_delta_ti)))
  sns.heatmap(
      sat_delta_ti,
      linewidths=0,
      cmap='RdBu_r',
      vmin=-vlim,
      vmax=vlim,
      xticklabels=False, cbar=cbar,
      ax=ax)
  ax.yaxis.set_ticklabels('ACGT', rotation='horizontal')

def plot_seq(si, distance, outfile=None):
    seqlen = 12288
    plot_start = 12288-distance
    plot_end = seqlen

    seqs_tmp =            seqs[si, plot_start:plot_end]
    scores_tmp = target_scores[si, plot_start:plot_end]
    ref_scores = scores_tmp[seqs_tmp]
    ref_scores = np.repeat(ref_scores[:,np.newaxis], 4, axis=1)

    # difference from reference
    delta_ti = scores_tmp - ref_scores

    # compute loss and gain
    delta_mean = delta_ti.mean(axis=1)
    delta_loss = delta_ti.min(axis=1)
    delta_gain = delta_ti.max(axis=1)

    print(si)
    print(dna_io.hot1_dna(seqs_tmp))

    f, axs = plt.subplots(figsize=(40, 10),nrows=6, ncols=1)
    plot_seqlogo(axs[0], seqs_tmp, delta_mean, pseudo_pct=0) # plot seqlogo, mean
    plot_seqlogo(axs[1], seqs_tmp, -delta_loss, pseudo_pct=0) # plot seqlogo, loss
    plot_seqlogo(axs[2], seqs_tmp, delta_gain, pseudo_pct=0) # plot seqlogo, gain
    plot_sad(axs[3], delta_loss, delta_gain) # plot loss and gain
    plot_heat(axs[4], delta_ti.T, 0.01) #, cbar=False plot heatmap
    plot_heat(axs[5], delta_ti.T, 0.01, cbar=True) #, cbar=False plot heatmap
    plt.tight_layout()

    if outfile is not None:
        f.savefig(outfile,format="pdf")
        plt.close()

# 10 random examples
np.random.seed(5)
ids = np.random.choice(500, 20, replace=False)
# ids = [136]
print(ids)
for i in ids:
    plot_seq(i, 200, outfile='png/example_%d.pdf'%i)
    # plot_seq(i, 725, outfile='png/example_%d.pdf'%i)
