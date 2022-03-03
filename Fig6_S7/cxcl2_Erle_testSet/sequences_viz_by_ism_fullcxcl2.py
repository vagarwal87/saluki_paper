#!/usr/bin/env python

import h5py
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from basenji.plots import *
from basenji.bin.basenji_sat_plot import *
from basenji import dna_io
import pandas as pd

# take the scores
h5 = h5py.File('ISM/f0_scores.h5','r')
seqs = h5['seqs'][:]
scores = h5['ism'][:]
h5.close()

cov_index = 0
target_scores = scores[:,:,:,cov_index]

def get_vals(si, distance, outfile=None):
    buffer=32
    seqlen = 12288
    plot_start = 12288-distance-buffer
    plot_end = seqlen-buffer

    seqs_tmp =            seqs[si, plot_start:plot_end]
    scores_tmp = target_scores[si, plot_start:plot_end]
    ref_scores = scores_tmp[seqs_tmp]
    ref_scores = np.repeat(ref_scores[:,np.newaxis], 4, axis=1)

    # difference from reference
    delta_ti = scores_tmp - ref_scores
    print(dna_io.hot1_dna(seqs_tmp))

    # compute mean effect of mutation
    delta_mean = delta_ti.mean(axis=1)
    df = pd.DataFrame({'Position':pd.Series(np.flip(np.arange(74962788,74963462)), dtype='int'),
                       'delta_HL_pred':pd.Series(delta_mean, dtype='float')})
    print(df)
    df = df.sort_values('Position')
    df['delta_HL_pred'] = df['delta_HL_pred'].rolling(window=8, center=True).mean()
    df = df.dropna()
    print(df)
    df.to_csv("cxcl2_predictions.wig", sep = "\t", header=False, index=False)

# 10 random examples
i = 136
get_vals(i, 674)
