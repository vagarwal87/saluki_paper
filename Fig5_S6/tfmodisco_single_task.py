# this run in tf114-gpu envir

import matplotlib
matplotlib.use('pdf')
import os, sys
import argparse
from collections import OrderedDict
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow
import modisco
from modisco.visualization import viz_sequence
import glob

######################################################################
##### Script modified from one written by Han Yuan, Calico Labs
######################################################################

##########
# inputs #
##########
parser = argparse.ArgumentParser(description='tfmodisco wrapper for custom scores.')

parser.add_argument('--out', dest='out', default='out',
                    help='output directory')
parser.add_argument('--mspmc', dest='mspmc', type=int,
                    default=20000,
                    help='maximum seqlets per metacluster')
parser.add_argument('--target_id', dest='id', type=int,
                    default=0,
                    help='target task id')
parser.add_argument('--region', dest='region',
                    default='full', help='region = full, 5utr, orf, 3utr')

args = parser.parse_args()

out_dir = args.out
mspmc = args.mspmc
target_id = args.id
region = args.region

files = glob.glob("ISM/*scores.h5")
seqs = []
coding = []
scores = []

for f in files:
    h5 = h5py.File(f,'r')
    seqs.append(h5['seqs'][:])
    coding.append(h5['coding'][:])
    try:
        scores.append(h5['scores'][:,:,:,target_id])
    except:
        scores.append(h5['ism'][:,:,:,target_id])
    h5.close()

seqs = np.concatenate(seqs)
coding = np.concatenate(coding)
scores = np.concatenate(scores)

os.makedirs(out_dir, exist_ok=True)
os.chdir(out_dir)

# run tfmodisco
target_importance = OrderedDict()
target_hypothetical = OrderedDict()

emptyidx=[]
locations=[]

if region != 'full':
    print(region)
    for i in range(scores.shape[0]): # iterate through all sequences
        myscore = np.sum(scores[i,], axis=1)
        txstart = np.where(myscore!=0)[0][0]
        txend = 12288
        orf = np.where(coding[i,]==True)[0]
        orfstart = orf[0]
        orfend = min(orf[-1]+3, txend)
        if region == '5utr':
            idxs=np.arange(txstart, orfstart)
        elif region == 'orf':
            idxs=np.arange(orfstart, orfend)
        else:
            idxs=np.arange(orfend, txend)
        if len(idxs) < 20:
            emptyidx.append(i)
        else:
            locations.append(idxs)

print("Threw away following indices:")
print(emptyidx)

scores = np.delete(scores, emptyidx, axis=0)
seqs = np.delete(seqs, emptyidx, axis=0)

##########################
# run tfmodisco workflow #
##########################
tfm_results = modisco.tfmodisco_workflow.workflow.TfModiscoWorkflow(
                  sliding_window_size=8,
                  flank_size=8,
                  min_metacluster_size=20,
                  target_seqlet_fdr=0.1,
                  max_seqlets_per_metacluster=mspmc, # don't put constrain on this
                  seqlets_to_patterns_factory=
                      modisco.tfmodisco_workflow.seqlets_to_patterns.TfModiscoSeqletsToPatternsFactory(
                          n_cores=16, # use 16 cores
                          trim_to_window_size=30,
                          # min_num_to_trim_to=15, #added later
                          initial_flank_to_add=10,
                          kmer_len=8, num_gaps=3,
                          num_mismatches=2,
                          final_min_cluster_size=30)
                  )(
                task_names=["task"], #list(target_importance.keys()),
                contrib_scores={"task": [x[idxs] for (x,idxs) in zip(seqs * scores,locations)]}, #target_importance,
                hypothetical_contribs={"task": [x[idxs] for (x,idxs) in zip(scores,locations)]}, #target_hypothetical,
                revcomp=False,
                one_hot=[x[idxs] for (x,idxs) in zip(seqs,locations)]) #seqs)

h5_out = h5py.File('tfmodisco_out_' + region + '.h5', 'w')
tfm_results.save_hdf5(h5_out)
h5_out.close()
