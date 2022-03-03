import tensorflow as tf
import argparse, json

from basenji.basenji import rnann

if tf.__version__[0] == '1':
  tf.compat.v1.enable_eager_execution()

##########
# inputs #
##########
parser = argparse.ArgumentParser(description='spike in motifs at different positions and evaluate effect on prediction.')
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
model = seqnn_model.get_model()
# model._layers = model._layers[:-1].

tf.keras.utils.plot_model(model, 'png/model.pdf', show_shapes=True, show_layer_names=False)
