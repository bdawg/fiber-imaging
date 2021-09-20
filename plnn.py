"""
A class to do neural network things, mostly with simulated data from lanterfiber.py.
This should probably be merged with the plwfs class.
"""

import matplotlib.pyplot as plt
import numpy as np
from tensorflow import keras
from sklearn.model_selection import train_test_split
plt.ion()

class plnn:
    def __init__(self, datadir=None):
        self.datadir = datadir
        self.datafilename = None
        self.pdict = {}
        self.savefilename = None
        self.model = None
        self.X_train = None
        self.X_test = None
        self.y_train = None
        self.y_test = None
        self.history_val_loss = None
        self.history_loss = None
        self.hparam_scan_results = {}

        self.init_pdict()


    def init_pdict(self, show_pdict=True):
        self.pdict['actFunc'] = 'relu'
        self.pdict['batchSize'] = 32
        self.pdict['learningRate'] = 0.0001
        self.pdict['lossFunc'] = 'mean_squared_error'
        self.pdict['n_units'] = 2000
        self.pdict['epochs'] = 100
        self.pdict['dropout_rate'] = 0.2
        self.pdict['bottleneck'] = 25
        if show_pdict:
            print(self.pdict)


    def modify_pdict(self, key, value, show_pdict=True, ndata_limit = None):
        self.pdict[key] = value
        if show_pdict:
            print(self.pdict)


    def load_data(self, datafilename, test_split=0.2, indata_key='simdata_input_arr',
                  outdata_key='simdata_outputI_arr', ndata_limit = None, shuffle=False):
        self.datafilename = datafilename
        dfile = np.load(self.datadir+datafilename)
        data_input = dfile[indata_key]
        data_output = dfile[outdata_key]
        if ndata_limit is not None:
            data_input = data_input[:ndata_limit, :]
            data_output = data_output[:ndata_limit, :]
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(data_input, data_output,
                test_size=test_split, shuffle=shuffle)


    def build_network(self):
        input_dim = self.X_train.shape[1]
        output_dim = self.y_train.shape[1]

        self.model = keras.Sequential([
            keras.layers.Dense(self.pdict['n_units'], activation=self.pdict['actFunc'], input_shape = (input_dim,)),
            # keras.layers.Dropout(self.pdict['dropout_rate']),
            keras.layers.Dense(self.pdict['bottleneck'], activation=self.pdict['actFunc']),
            # keras.layers.Dropout(pdict['dropout_rate']),
            keras.layers.Dense(self.pdict['n_units'], activation=self.pdict['actFunc']),
            # keras.layers.Dropout(self.pdict['dropout_rate']),
            keras.layers.Dense(output_dim)
        ])

        self.model.summary()
        self.model.compile(optimizer=keras.optimizers.Adam(learning_rate=self.pdict['learningRate']),
                      loss=self.pdict['lossFunc'])


    def fit_network(self, save_output=False, savefilename=None, show_plots=True, wg_to_plot=0,
                    n_predictplot = 100):
        history = self.model.fit(self.X_train, self.y_train, validation_data=(self.X_test, self.y_test),
                            epochs=self.pdict['epochs'], batch_size=self.pdict['batchSize'])

        predictions = self.model.predict(self.X_test[:n_predictplot, :])
        self.history_loss = history.history['loss']
        self.history_val_loss = history.history['val_loss']

        if save_output:
            self.savefilename = savefilename
            np.savez(self.datadir+self.savefilename+'.npz', predictions=predictions,
                     y_test=self.y_test[:n_predictplot,0], datafilename=self.datafilename,
                     history_loss=self.history_loss, history_val_loss=self.history_val_loss, pdict=self.pdict)
            self.model.save(self.datadir+self.savefilename+'.h5')

        if show_plots:
            # Plot training & validation loss values
            plt.figure(1)
            plt.plot(history.history['loss'])
            plt.plot(history.history['val_loss'])
            plt.title('Model loss')
            plt.ylabel('Loss')
            plt.xlabel('Epoch')
            plt.legend(['Train', 'Test'], loc='upper left')
            plt.show()

            # Plot real vs predicted data
            plt.figure(2)
            plt.clf()
            plt.plot(self.y_test[:n_predictplot,wg_to_plot], alpha=0.7)
            plt.plot(predictions[:,wg_to_plot], alpha=0.7)


    def scan_hyperparam(self, key, network_savefilename_prefix=None, scan_values=None, value_min=None, value_max=None,
                        n_values=None, show_plots=True, wg_to_plot=0, n_predictplot=100, save_network_output=False,
                        scan_outfilename='scantemp'):
        if scan_values is not None:
            all_hparam_vals = scan_values
            n_values = len(all_hparam_vals)
        else:
            all_hparam_vals = np.linspace(value_min, value_max, n_values)
        all_training_loss_history = []
        all_val_loss_history = []
        all_best_training_loss = []
        all_best_val_loss = []
        for k in range(n_values):
            curval = all_hparam_vals[k]
            print('Starting iteration %d, with value %f' % (k, curval))
            if save_network_output:
                filename = network_savefilename_prefix + '_%.3d'%k
                print('Saving to file ' + filename)
            else:
                filename = None
            self.modify_pdict(key, curval)
            self.build_network()
            self.fit_network(save_output=save_network_output, savefilename=filename, show_plots=show_plots,
                             wg_to_plot=wg_to_plot, n_predictplot=n_predictplot)

            all_training_loss_history.append(self.history_loss)
            all_val_loss_history.append(self.history_val_loss)
            all_best_training_loss.append(np.min(self.history_loss))
            all_best_val_loss.append(np.min(self.history_val_loss))

            if scan_outfilename is not None:
                np.savez(self.datadir+scan_outfilename+'_TEMPSAVE_'+'.npz',
                         all_training_loss_history=all_training_loss_history, all_val_loss_history=all_val_loss_history,
                         all_best_training_loss=all_best_training_loss, all_best_val_loss=all_best_val_loss,
                         all_hparam_vals=all_hparam_vals)

        if scan_outfilename is not None:
            np.savez(self.datadir+scan_outfilename+'.npz',
                     all_training_loss_history=all_training_loss_history, all_val_loss_history=all_val_loss_history,
                     all_best_training_loss=all_best_training_loss, all_best_val_loss=all_best_val_loss,
                     all_hparam_vals=all_hparam_vals)

        if show_plots:
            plt.figure(3)
            plt.plot(all_hparam_vals, all_best_val_loss)
            plt.xlabel('Hyperparamter value')
            plt.ylabel('Best validation loss')

        self.hparam_scan_results = {}
        self.hparam_scan_results['all_training_loss_history'] = all_training_loss_history
        self.hparam_scan_results['all_val_loss_history'] = all_val_loss_history
        self.hparam_scan_results['all_best_training_loss'] = all_training_loss_history
        self.hparam_scan_results['all_best_val_loss'] = all_best_val_loss
        self.hparam_scan_results['all_hparam_vals'] = all_hparam_vals


