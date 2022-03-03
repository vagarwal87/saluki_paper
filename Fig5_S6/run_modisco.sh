#!/bin/bash

#SBATCH -p gpu
#SBATCH -n 16
#SBATCH --gres=gpu:gtx1080ti:1
#SBATCH --mem 200G
#SBATCH -J orf
#SBATCH -o orf.%j.out

######################################################################
##### Script modified from one written by Han Yuan, Calico Labs
######################################################################

source /home/yuanh/.bashrc
source activate tf114-gpu

python tfmodisco_single_task.py --region orf --out orf_scores

cd orf_scores
python ../modisco2meme.py --h5 tfmodisco_out_orf.h5
