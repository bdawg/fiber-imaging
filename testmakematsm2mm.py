
from viewRSoftData import *
import matplotlib.pyplot as plt
import numpy as np
plt.ion()


rsoft_datadir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix_19cPLJ21_allmodesrun20210802/'
rsoft_fileprefix = '19cPLJ21_sm2mm_srun20210802_wg'
processed_datadir = '/Users/bnorris/DataAnalysis/fiber-imaging/'
nwgs = 19
show_plots = True
av_fluxes = 100 #100 # If >0, average n flux measurements from monitor
offset_sm_meas = 100 # 100 # Measure SM fluxes starting at this index.
                     # Useful to skip first part of monitor since still coupling.


# Specify relevant indices of MONdata
sm_power_monrange = (0, 19)
mm_power_monrange = (19, 38)
mm_phase_monrange = (38, 57)




# Extract required data from rsoft files
all_smpowers = []
all_mmpowers = []
all_mmphases = []

for wgnum in range(nwgs):
    rsoft_filename = rsoft_fileprefix + '%.2d' % (wgnum+1)
    print('Reading rsoft files ' + rsoft_filename)
    r = Rsoftdata(rsoft_datadir)
    r.readall(filename=rsoft_filename)
    if show_plots:
        r.plotall()
        plt.pause(0.001)

    smpower_mons = r.MONdata[:, sm_power_monrange[0]:sm_power_monrange[1]]
    mmpower_mons = r.MONdata[:, mm_power_monrange[0]:mm_power_monrange[1]]
    mmphase_mons = r.MONdata[:, mm_phase_monrange[0]:mm_phase_monrange[1]]
    if av_fluxes > 0:
        f = smpower_mons[offset_sm_meas:av_fluxes+offset_sm_meas, :]
        smpower = f.mean(axis=0)
        f = mmpower_mons[-av_fluxes:, :]
        mmpower = f.mean(axis=0)
    else:
        smpower = smpower_mons[offset_sm_meas, :]
        mmpower = mmpower_mons[-1, :]
    mmphase = mmphase_mons[-1, :]

    all_smpowers.append(smpower)
    all_mmpowers.append(mmpower)
    all_mmphases.append(mmphase)

outfilename = processed_datadir + 'extractedvals_' + rsoft_fileprefix + '.npz'
np.savez(outfilename, all_smpowers=all_smpowers, all_mmpowers=all_mmpowers, all_mmphases=all_mmphases)

if show_plots:
    plt.figure(2)
    plt.clf()
    plt.subplot(211)
    plt.imshow(np.asarray(all_mmpowers))
    plt.colorbar()
    plt.title('Output mode power')
    plt.ylabel('Excited waveguide no.')
    plt.xlabel('Mode no.')
    plt.subplot(212)
    plt.imshow(np.asarray(all_mmphases), cmap='twilight_shifted')
    plt.colorbar()
    plt.title('Output mode phase')
    plt.ylabel('Excited waveguide no.')
    plt.xlabel('Mode no.')
    plt.tight_layout()


