import matplotlib.pyplot as plt
import numpy as np
from lanternfiber import lanternfiber
from tensorflow import keras
from sklearn.model_selection import train_test_split
from plnn import plnn
plt.ion()


datadir = '../../data/fiber-imaging/sim_data/'

# datafilename = 'simdata_amp0-1_ph0-2pi_01.npz'
# # savefilename = 'nnsimdata_a0-1p0-2p01_2x1000_n25'
# savefilename = 'HPSCAN1_nnsimdata_a0-1p0-2p01_2x1000_n25'

datafilename = 'simdata_amp0-1_ph0-2pi_01_19FakeOuts.npz'
# savefilename = 'TEST_nnsimdata_19fk_a0-1p0-2p01_2x1000_n100'
savefilename = 'HPSCAN1_nnsimdata_a0-1p0-2p01_2x2000_19FakeOuts'

save_output = True

epochs = 200

ndata_limit = None
n_predictplot = 100

pn = plnn(datadir)
pn.modify_pdict('epochs', epochs)
# pn.modify_pdict('n_units', 100)

pn.load_data(datafilename, ndata_limit=ndata_limit)
# pn.build_network()
# pn.fit_network(save_output=save_output, savefilename=savefilename, show_plots=True, n_predictplot=n_predictplot)

# Allow about 3s/epoch, or 10 mins per 200 epoch.
scan_key = 'bottleneck'
# scan_values = [100, 50, 30, 20, 15, 10]
# scan_values = [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100, 90, 80, 70, 60, 50, 45, 40, 35, 30, 25, 20,
#                15, 10, 5, 1]
# scan_values = np.arange(50,0,-1)
scan_values = np.array([1000, 900, 800, 700, 600, 500, 400, 300, 200, 100, 90, 80, 70, 60])
scan_values = np.hstack((scan_values, np.arange(50,0,-1)))

network_savefilename_prefix = savefilename + '_HPSCAN'
scan_outfilename = savefilename + '_HPScanResults'
pn.scan_hyperparam(scan_key, scan_values=scan_values, save_network_output=True,
                   network_savefilename_prefix=network_savefilename_prefix,
                   scan_outfilename=scan_outfilename)
