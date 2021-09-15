"""
Generate a set of simulated PL test data, using a transfer matrix built from saved rsoft probe measurements.
The input data will be a vector of input complex amplitudes (saved as a vector of alternating real and imaginary
components) and the output data is the wavegude ouput *intensities*.
"""



import matplotlib.pyplot as plt
import numpy as np
from lanternfiber import lanternfiber
plt.ion()

datadir = '/Users/bnorris/DataAnalysis/fiber-imaging/'
data_filename = 'extractedvals_probeset_19LP__Good202107.npz'

f = lanternfiber(datadir=datadir)
f.load_savedvalues(data_filename)
f.make_transfer_matrix_mm2sm(show_plots=False)


# Make simulated data
ndata = 100000
outfilename = 'sim_data/temp_simdata_amp0-1_ph0-2pi_01.npz'
amp_range = [0,1]
phase_range = [0, 2*np.pi]

f.generate_sim_I_data(ndata=ndata, amp_range=amp_range, phase_range=phase_range, outfilename=outfilename)




