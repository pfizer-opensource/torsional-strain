import os, sys
import math
import argparse
import shutil

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold

from keras.models import Sequential
from keras.layers import Dense, Dropout, LocallyConnected1D, Activation, \
    GaussianNoise, GaussianDropout
from keras.layers.normalization import BatchNormalization
from keras.wrappers.scikit_learn import KerasRegressor
from keras.utils import multi_gpu_model
from keras.callbacks import EarlyStopping
from keras.callbacks import ModelCheckpoint
from keras.optimizers import Adam
from keras.models import load_model
from keras.callbacks import Callback

import timeit
import pickle

from openeye import oechem

from torsion.model import get_sf_elements
from torsion.analysis import get_dihedral_inchi_key

import matplotlib.pyplot as plt

# fix random seed for reproducibility
seed = 7
np.random.seed(seed)


def get_model(num_feat=294, lr=1e-3, drop_out=0.1, layer_dims=''):
    model = Sequential()
    act_fn = 'relu'

    if len(layer_dims) == 0:
        layer_dims = [10, 5, 0.2]
    else:
        layer_dims = [float(d) for d in layer_dims.split('-')]

    model.add(
        Dense(
            int(num_feat * layer_dims[0]), input_dim=num_feat,
            kernel_initializer='normal'))
    model.add(Activation(act_fn))
    model.add(BatchNormalization())
    model.add(Dropout(drop_out))

    for layer_dim in layer_dims[1:-1]:
        model.add(Dense(int(num_feat * layer_dim)))
        model.add(Activation(act_fn))
        model.add(BatchNormalization())
        model.add(Dropout(drop_out))

    model.add(Dense(int(num_feat * layer_dims[-1])))
    model.add(Activation(act_fn))
    model.add(Dropout(drop_out))

    model.add(Dense(1))

    adam = Adam(lr=lr)
    model.compile(loss='logcosh', optimizer=adam)

    return model


ENERGY_KEY = 'ENERGY'
INCHI_KEY = 'Inchi'

def generate_training_input(mol_file):
    '''


    :param mol_file: str
    :return: pd.DataFrame
    '''
    ifs = oechem.oemolistream(mol_file)
    training_data = []
    for mol in ifs.GetOEGraphMols():
        energy = float(oechem.OEGetSDData(mol, ENERGY_KEY))
        sf_elements = get_sf_elements(mol)
        dihe_inchi = get_dihedral_inchi_key(mol)

        data = [dihe_inchi, energy]
        data.extend(sf_elements)
        training_data.append(data)

    ifs.close()

    columns = [INCHI_KEY, ENERGY_KEY]
    num_sf_elements = len(training_data[0]) - 2
    sf_columns = ['sf_%d'%(i+1) for i in range(num_sf_elements)]
    columns.extend(sf_columns)

    df = pd.DataFrame(training_data, columns=columns)

    # calculate relative energy for each profile
    grouped = df.loc[:,[INCHI_KEY, ENERGY_KEY]].groupby(INCHI_KEY)
    df2 = grouped.transform(lambda x: x - x.min())
    df[ENERGY_KEY] = df2[ENERGY_KEY]

    return df



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Train neural network model to predict torsional relative energy')
    parser.add_argument('--input', type=str, help='sd file containing MM structures alongwith '
                                                  'sd properties with torsion atom indices and QM energy')
    parser.add_argument('--num_epoch', default=5000, type=int, help='number of epoch (default = 2000)')
    parser.add_argument('--batch_size', default=256, type=int, help='batch size (default: 256)')
    parser.add_argument('--layer_dims', default='10-5-1-0.2', type=str, help='layer dimensions')
    parser.add_argument('--lr', default=0.0001, type=float, help='learning rate (default: 1e-r)')
    parser.add_argument('--dropout', default=0.2, type=float, help='dropout (default: 0.2)')
    parser.add_argument('--val_split', default=0.1, type=float, help='validation split (default: 0.1)')

    parser.add_argument('--scalar', default='scaler.pkl', type=str, help='output file with standard scaler')
    parser.add_argument('--model', default='model.h5', type=str, help='output file with trained model')

    parser.add_argument('-v', '--verbose', action='count', default=0)
    args = parser.parse_args()

    input_file = args.input

    num_epoch = args.num_epoch
    batch_size = args.batch_size
    lr = args.lr
    dropout = args.dropout
    layer_dims = args.layer_dims

    # generate training data using the molecules in the input file
    # for each molecule in the input file, extract the QM energy from SD property "ENERGY"
    # and generate symmetry function elements around the specified torsion (SD property "TORSION_ATOMS_FRAGMENT")
    df = generate_training_input(input_file)

    # cap the relative energy
    tmp_idx = df.ENERGY > 30
    df.ENERGY[tmp_idx] = 30.0 + np.exp(30 - df.ENERGY[tmp_idx])

    dihe_inchis = df[INCHI_KEY].unique()
    print('Number of profiles: %d' % len(dihe_inchis))

    desc_bgn_idx = df.columns.get_loc('sf_1')

    Xtrain = df.as_matrix(columns=df.columns[desc_bgn_idx:])
    ytrain = df.ENERGY

    # feature transformation
    scaler = StandardScaler().fit(Xtrain)
    Xtrain = scaler.transform(Xtrain)

    print('Xtrain.shape ', Xtrain.shape)

    # save feature transformation
    with open(args.scalar, 'wb') as fptr:
        pickle.dump(scaler, fptr)

    _, num_feat = Xtrain.shape

    # early stopping criteria
    earlystop = EarlyStopping(monitor='val_loss', min_delta=0.001, patience=100, \
                              verbose=1, mode='auto')

    model_file = args.model
    # create DNN model
    model = get_model(num_feat, lr, dropout, layer_dims)

    print(model.summary())

    checkpointer = ModelCheckpoint(
        filepath=model_file, verbose=1, save_best_only=True)
    callbacks_list = [checkpointer]

    # train DNN model
    model.fit(
        Xtrain,
        ytrain,
        epochs=num_epoch,
        batch_size=batch_size,
        validation_split=args.val_split,
        callbacks=callbacks_list,
        verbose=1)

    print('Training complete')
    print('Standard scalar is saved in %s' % args.scalar)
    print('Model is saved in %s' % args.model)




