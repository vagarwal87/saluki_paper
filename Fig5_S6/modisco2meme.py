import os, sys
import argparse
import h5py
import numpy as np
import pandas as pd
import pickle
import subprocess

######################################################################
##### Script modified from one written by Han Yuan, Calico Labs
######################################################################

##########
# inputs #
##########
parser = argparse.ArgumentParser(description='convert tfmodisco output to meme.')
parser.add_argument('--h5', dest='h5',
                    help='hdf5 for tfmodisco output.')
parser.add_argument('--thres', dest='threshold',
                    default=0.5,
                    help='information content threshold for trimming the long pwms.')
parser.add_argument('--out', dest='out',
                    default='out',
                    help='output directory')
parser.add_argument('--meme_db', dest='db',
                   default='/home/vagar/databases/meme_databases/motif_databases/CISBP-RNA/Homo_sapiens.dna_encoded.meme',
                   help='meme database')


args = parser.parse_args()
h5 = args.h5
threshold = args.threshold
out = args.out
meme_db = args.db
os.makedirs(out, exist_ok=True) # create directory
background = np.array([0.27, 0.23, 0.23, 0.27])
# background = np.array([0.276, 0.217, 0.219, 0.288])

h5_res = h5py.File(h5,'r')
metacluster_names = [x.decode("utf-8") for x in list(h5_res["metaclustering_results"]["all_metacluster_names"][:])]

def ic_clip(pwm, background, threshold):
    odds_ratio = ((pwm+0.001)/(1.004))/(background[None,:])
    ic = (np.log((pwm+0.001)/(1.004))/np.log(2))*pwm - (np.log(background)*background/np.log(2))[None,:]
    ic_total = np.sum(ic,axis=1)[:,None]

    # no bp pass threshold
    if ~np.any(ic_total.flatten()>threshold): return None

    left = np.where(ic_total>threshold)[0][0]
    right = np.where(ic_total>threshold)[0][-1]
    return pwm[left:(right+1)]

all_pwms = dict()
for metacluster_name in metacluster_names:
    print(metacluster_name)
    metacluster_grp = (h5_res["metacluster_idx_to_submetacluster_results"][metacluster_name])
    print("activity pattern:", metacluster_grp["activity_pattern"][:])
    all_patterns = metacluster_grp["seqlets_to_patterns_result"]["patterns"]["all_pattern_names"][:]
    all_pattern_names = [x.decode("utf-8") for x in list(all_patterns)]
    if (len(all_pattern_names)==0):
        print("No motifs found in this metacluster")
    for pattern_name in all_pattern_names:
        pattern_id = (metacluster_name+'_'+pattern_name)
        pattern = metacluster_grp["seqlets_to_patterns_result"]["patterns"][pattern_name]
        fwd = np.array(pattern["sequence"]["fwd"])
        clip_pwm = ic_clip(fwd, background, threshold)
        if clip_pwm is None or clip_pwm.shape[0] < 3:
            print('pattern_id: %s is skipped because motif too small.' %pattern_id)
        else:
            all_pwms[pattern_id] = clip_pwm
            print('pattern_id: %s is converted to pwm.' %pattern_id)

# save pwm to pickle
f = open(out+'/pwms.pkl', 'wb')
pickle.dump(all_pwms, f)
f.close()

# save to meme
meme_file = out+'/pwms.meme'
with open(meme_file,'w+') as f:
    f.write('MEME version 4\n')
    f.write('\n')
    f.write('ALPHABET= ACGT\n')
    f.write('\n')
    f.write('strands: + -\n')
    f.write('\n')
    f.write('Background letter frequencies\n')
    f.write('A 0.27 C 0.23 G 0.23 T 0.27\n')
    f.write('\n')

for key in all_pwms.keys():
    with open(meme_file,'a') as f:
        f.write('MOTIF '+key+'\n')
        f.write('letter-probability matrix:\n')
    with open(meme_file,'ab') as f:
        np.savetxt(f, all_pwms[key])
    with open(meme_file,'a') as f:
        f.write('\n')

# run tomtom
subprocess.call('tomtom -norc -dist pearson -thresh 0.1 -oc %s %s %s' %(out, meme_file, meme_db), shell=True)
