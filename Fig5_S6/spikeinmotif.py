import os, sys
import argparse, json, h5py, time
import numpy as np
import tensorflow as tf
import deeplift.dinuc_shuffle as ds
from basenji import dataset
from basenji import plots
from basenji import dna_io
try:
  import rnann
except:
  from basenji import rnann
# from basenji import trainer

if tf.__version__[0] == '1':
  tf.compat.v1.enable_eager_execution()


#python spikeinmotif.py --region 5utr --seqs ISM/f3_c2_scores.h5 train_gru/params.json \
#   train_gru/f8_c4/train/model0_best.h5 >spikeinmotif_predictions_5utr_f8c4.txt

#python spikeinmotif.py --codons --frame 0 --region orf --seqs ISM/f3_c2_scores.h5 train_gru/params.json \
#   train_gru/f8_c4/train/model0_best.h5 >spikeinmotif_predictions_orf_codons_frame0.txt

# for x in {0..9}; do { for y in {0..4}; do { REG="5utr"; echo $x, $y; sbatch --mem 20000 -J S$REG$x$y -e spikeinmotif/splice$REG\_f$x\_c$y.err -o spikeinmotif/splice$REG\_f$x\_c$y.txt --wrap=". /home/drk/anaconda3/etc/profile.d/conda.sh; conda activate tf2.6-rna; python spikeinmotif.py --addsplicesites --region $REG --seqs ISM/f$x\_c$y\_scores.h5 train_gru/params.json train_gru/f$x\_c$y/train/model0_best.h5"; } done } done


##########
# inputs #
##########
parser = argparse.ArgumentParser(description='spike in motifs at different positions and evaluate effect on prediction.')

parser.add_argument('--seqs', dest='seqs',
                    help='hdf5 for seqs.')
parser.add_argument('--region', dest='region',
                    default='3utr', help='region = 5utr, orf, 3utr')
parser.add_argument('--frame', dest='frame',
                    default=0, help='frame = 0, 1, 2', type=int)
parser.add_argument('--shuf', dest='shuf',
                    default=False, help='shuffle or not', action="store_true")
parser.add_argument('--codons', dest='codons',
                    default=False, help='use codons, exluding stop codons [Modisco motifs by default]', action="store_true")
parser.add_argument('--stopcodons', dest='stopcodons',
                    default=False, help='use stop codons [Modisco motifs by default]', action="store_true")
parser.add_argument('--addsplicesites', dest='splicesites',
                    default=False, help='insert splice site only', action="store_true")
parser.add_argument(dest='pfile', help='params file')
parser.add_argument(dest='mfile', help='model file')


args = parser.parse_args()
frame = args.frame
codons = args.codons
stopcodons = args.stopcodons
shuf = args.shuf
input_seqs = args.seqs
region = args.region
splicesites = args.splicesites
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

#######################################################
# evaluation

if codons:
    motifs = []
    for a in ['A','T','G','C']:
        for b in ['A','T','G','C']:
            for c in ['A','T','G','C']:
                motifs.append(a+b+c)
    motifs.remove('TAG')
    motifs.remove('TAA')
    motifs.remove('TGA')
elif stopcodons:
    motifs = ['TAG', 'TAA', 'TGA']
elif splicesites:
    motifs = ["SS"]
else:
    motifs = ["GGACT","CCCCC","TGTACAT","TATTTAT","AATAAA", "ATTATTA"]


# os.makedirs(out_dir, exist_ok=True)
# os.chdir(out_dir)

# print(input_seqs)
h5 = h5py.File(input_seqs)
seqs = h5['seqs']
splice = h5['splice']
coding = h5['coding']

# print('%s\t%s\t%s\t%s\t%s' % ("seq", "shuf_num", "motif", "preds"))
# print(seqs.shape[0])

def left_justify(data, txstart, txend):
    tmp = np.copy(data[txstart:])
    data[:,:] = 0
    data[0:(txend-txstart),:] = tmp
    return data

bins = 50
numshufs = 1

count = 0
for i in range(seqs.shape[0]): # iterate through all sequences
    if shuf:
        numshufs = 10
        if count==10:
           break
    myscore = np.sum(seqs[i,], axis=1)
    txstart = np.where(myscore!=0)[0][0]
    txend = 12288
    orf = np.where(coding[i,]==True)[0]
    orfstart = orf[0]
    orfend = min(orf[-1]+3, txend)
    idxs5utr=np.arange(txstart, orfstart)
    idxsorf=np.arange(orfstart, orfend)
    idxs3utr=np.arange(orfend, txend)
    if region == '5utr':
        idxs = idxs5utr
    elif region == 'orf':
        idxs = idxsorf
    else:
        idxs = idxs3utr
    if len(idxs5utr) >= 100 and len(idxsorf) >= 500 and len(idxs3utr) >= 500:
        count += 1
        if shuf:
            myshufs = ds.dinuc_shuffle(seqs[i,idxs], numshufs) #10 dinuc shuffles
        else:
            myshufs = seqs[i,idxs]
        for j in range(numshufs):
            wt = seqs[i,]
            if shuf:
                wt[idxs] = myshufs[j]
            for m in motifs:
                batch = np.zeros((bins+1,txend,6))
                data = np.concatenate((wt, np.expand_dims(coding[i,], -1), np.expand_dims(splice[i,], -1)), axis=-1)
                data = data.astype('float32')
                data = left_justify(data, txstart, txend)
                batch[0] = np.expand_dims(data, 0)
                for p in range(bins):
                    shuffseq = np.copy(wt)
                    splicetrack = np.copy(splice[i,])
                    codontrack = np.copy(coding[i,])
                    regionpos = int(p/bins * (idxs[-1]-idxs[0]) )
                    if codons or stopcodons:
                        thisframe = regionpos % 3
                        while thisframe != frame:
                            regionpos += 1
                            thisframe = regionpos % 3
                    pos = idxs[0] + regionpos
                    if splicesites:
                        splicetrack[pos:(pos+2)] = True #insert 5'/3' splice sites
                    else:
                        shuffseq[pos:(pos+len(m))] = dna_io.dna_1hot(m) #replace with motif
                        if stopcodons:
                            codontrack[(pos+1):] = False
                    # print(dna_io.hot1_dna(shuffseq[idxs]))
                    data = np.concatenate((shuffseq, np.expand_dims(codontrack, -1), np.expand_dims(splicetrack, -1)), axis=-1)
                    data = data.astype('float32')
                    data = left_justify(data, txstart, txend) #left justify
                    batch[p+1] = np.expand_dims(data, 0)
                pred = seqnn_model.predict(batch)
                print('seq%s\t%s\t%s' % (count, j, m), end = '')
                for p in pred:
                    print('\t%s' % (p[0]), end = '')
                print('')
h5.close()
