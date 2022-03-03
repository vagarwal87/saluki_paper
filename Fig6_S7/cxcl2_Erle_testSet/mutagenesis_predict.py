import os, sys
import argparse, json, h5py, time
import numpy as np
import tensorflow as tf
from basenji import dataset
from basenji import dna_io
import pandas as pd
try:
  import rnann
except:
  from basenji import rnann

if tf.__version__[0] == '1':
  tf.compat.v1.enable_eager_execution()

MAXLEN = 12288

##########
# inputs #
##########
parser = argparse.ArgumentParser(description='In silico MPRA experiment')
parser.add_argument(dest='pfile', help='params file')
parser.add_argument(dest='mfile', help='model file')
args = parser.parse_args()

params_file = args.pfile #train_gru/params.json
model_file = args.mfile #train_gru/f0_c0/train/model0_best.h5

# read model parameters
with open(params_file) as params_open:
    params = json.load(params_open)
params_model = params['model']
params_train = params['train']

# initialize model
seqnn_model = rnann.RnaNN(params_model)
seqnn_model.restore(model_file)

construct = pd.read_table("BTV_construct.txt", index_col=0, header=None).values
aa_len = int(len(construct[1][0])/3)
coding = np.append(np.zeros(len(construct[0][0])), np.tile([1,0,0], aa_len))

reporter = construct[0]+construct[1]+construct[2]
seq = pd.read_table("cxcl2_hg19.fa", header=0).values[0][0]

batch = np.zeros((1,MAXLEN,6))
wtseq = (reporter+seq+construct[3])[0]
batch[0,0:len(wtseq),0:4] = dna_io.dna_1hot(wtseq)
batch[0,0:len(coding),4] = coding
wtpred = seqnn_model.predict(batch)
wtpred = wtpred[0][0]

print('%s\t%s\t%s' % ("pos", "mut", "pred"))
pos = 74963332
for i in range(len(seq)): # iterate through all positions
    for mut in ['A','C','G','T']:
        if seq[i] == mut: print('%s\t%s\t0' % (pos-i, mut))
        else:
            batch = np.zeros((1,MAXLEN,6))
            batch[0,0:len(wtseq),0:4] = dna_io.dna_1hot((reporter+seq[:i]+mut+seq[i+1:]+construct[3])[0])
            batch[0,0:len(coding),4] = coding
            pred = seqnn_model.predict(batch)
            print('%s\t%s\t%s' % (pos-i, mut, pred[0][0]-wtpred))
