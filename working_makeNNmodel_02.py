import matplotlib.pyplot as plt
import numpy as np
from lanternfiber import lanternfiber
from tensorflow import keras
from sklearn.model_selection import train_test_split
from plnn import plnn
plt.ion()


datadir = '../../data/fiber-imaging/sim_data/'

datafilename = 'simdata_amp0-1_ph0-2pi_01.npz'
savefilename = 'TEST_nnsimdata_a0-1p0-2p01_2x1000_n25'

# datafilename = 'simdata_amp0-1_ph0-2pi_01_19FakeOuts.npz'
# savefilename = 'TEST_nnsimdata_19fk_a0-1p0-2p01_2x1000_n100'

save_output = True

epochs = 200

ndata_limit = None
# ndata_limit = 1000
n_predictplot = 100

pn = plnn(datadir)
pn.modify_pdict('epochs', epochs)
pn.load_data(datafilename, ndata_limit=ndata_limit)
pn.build_network()
pn.fit_network(save_output=save_output, savefilename=savefilename, show_plots=True, n_predictplot=n_predictplot)

