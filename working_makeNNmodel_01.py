import matplotlib.pyplot as plt
import numpy as np
from lanternfiber import lanternfiber
from tensorflow import keras
from sklearn.model_selection import train_test_split
plt.ion()

# datadir = '/Users/bnorris/Dropbox/data/fiber-imaging/'
# datadir = '/Users/bnorris/Dropbox/data/fiber-imaging/'
datadir = '../../data/fiber-imaging/sim_data/'
datafilename = 'simdata_amp0-1_ph0-2pi_01.npz'


save_output = False
savefilename = 'nnsimdata_a0-1p0-2p01_2x1000_n18'

datafilename = 'simdata_amp0-1_ph0-2pi_01_19FakeOuts.npz'
savefilename = 'nnsimdata_a0-1p0-2p01_2x1000_n18'


ndata_limit = None
# ndata_limit = 1000
n_predictplot = 100

pdict = {}
pdict['actFunc'] = 'relu'
# pdict['actFunc'] = tf.keras.layers.LeakyReLU(alpha=0.3) #'relu'
pdict['batchSize'] = 32
pdict['learningRate'] = 0.0001
pdict['lossFunc'] = 'mean_squared_error'
pdict['n_units'] = 2000
pdict['epochs'] = 200
# pdict['dropout_rate'] = 0.2
pdict['bottleneck'] = 25



dfile = np.load(datadir+datafilename)
data_input = dfile['simdata_input_arr']
data_output = dfile['simdata_outputI_arr']

if ndata_limit is not None:
    data_input = data_input[:ndata_limit, :]
    data_output = data_output[:ndata_limit, :]

X_train, X_test, y_train, y_test = train_test_split(data_input, data_output, test_size=0.2, shuffle=False)
input_dim = X_train.shape[1]
output_dim = y_train.shape[1]

model = keras.Sequential([
    keras.layers.Dense(pdict['n_units'], activation=pdict['actFunc'], input_shape = (input_dim,)),
    # keras.layers.Dropout(pdict['dropout_rate']),
    keras.layers.Dense(pdict['bottleneck'], activation=pdict['actFunc']),
    # keras.layers.Dropout(pdict['dropout_rate']),
    keras.layers.Dense(pdict['n_units'], activation=pdict['actFunc']),
    # keras.layers.Dropout(pdict['dropout_rate']),
    keras.layers.Dense(output_dim)
])
model.summary()

model.compile(optimizer=keras.optimizers.Adam(learning_rate=pdict['learningRate']),
              loss=pdict['lossFunc'])


history = model.fit(X_train, y_train, validation_data=(X_test, y_test),
                    epochs=pdict['epochs'], batch_size=pdict['batchSize'])


predictions = model.predict(X_test[:n_predictplot, :])
history_loss = history.history['loss']
history_val_loss = history.history['val_loss']

if save_output:
    np.savez(datadir+savefilename+'.npz', predictions=predictions, y_test=y_test[:n_predictplot,0],
             datafilename=datafilename, history_loss=history_loss, history_val_loss=history_val_loss,
             pdict=pdict)
    model.save(datadir+savefilename+'.h5')

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
plt.plot(y_test[:n_predictplot,25], alpha=0.7)
plt.plot(predictions[:,25], alpha=0.7)
