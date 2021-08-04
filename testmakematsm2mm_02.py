
from lanternfiber import lanternfiber
import matplotlib.pyplot as plt
import numpy as np
plt.ion()


rsoft_datadir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix_19cPLJ21_allmodesrun20210802/'
rsoft_fileprefix = '19cPLJ21_sm2mm_srun20210802_wg'
processed_datadir = '/Users/bnorris/DataAnalysis/fiber-imaging/'
data_filename = 'extractedvals_19cPLJ21_sm2mm_srun20210802_wg.npz'
show_plots = False

f = lanternfiber(datadir=processed_datadir)
# f.load_rsoft_data_sm2mm(rsoft_datadir, rsoft_fileprefix, show_plots=show_plots)
f.load_savedvalues(data_filename)

# Make a matrix
# This requires that for each measurement, only one waveguide is excited
nwgs = len(f.all_smpowers)
Cmat = np.zeros((nwgs,nwgs), dtype=np.cfloat)
sm_phase = 0 # Set the phase of all SM waveguides (at injection) to zero
for col_num in range(nwgs):
    sm_val = f.all_smpowers[col_num][col_num] * np.exp(1j*sm_phase)
    col = f.all_mmpowers[col_num] * np.exp(1j*f.all_mmphases[col_num]/180*np.pi) / sm_val
    Cmat[:,col_num] = col


# Validate matrix:
plt.figure(1)
for n in range(nwgs):
    plt.clf()
    plt.subplot(211)
    plt.plot(f.all_mmpowers[n],'-x')
    plt.plot(np.abs(np.matmul(Cmat, f.all_smpowers[n])),'-+')
    plt.subplot(212)
    plt.plot(f.all_mmphases[n],'-x')
    newangs = np.angle(np.matmul(Cmat, f.all_smpowers[n]))/np.pi*180
    newangs[newangs < 0] = newangs[newangs < 0]+360
    plt.plot(newangs,'-+')
    plt.pause(0.001)
    a=input('')