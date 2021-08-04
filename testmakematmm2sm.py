
from viewRSoftData import *
import matplotlib.pyplot as plt
import numpy as np
plt.ion()


rsoft_datadir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix_19cPLJ21_MM2SM_run20210803/'
rsoft_fileprefix = '19cPLJ21_mm2sm_run20210803_'
processed_datadir = '/Users/bnorris/DataAnalysis/fiber-imaging/'
nwgs = 19 #19
show_plots = True
av_fluxes = 100 #100 # If >0, average n flux measurements from monitor
# offset_sm_meas = 100 # 100 # Measure SM fluxes starting at this index.
#                      # Useful to skip first part of monitor since still coupling.



# Specify mode indices
LP_modes = np.array([[0,1],
                     [0,2],
                     [0,3],
                     [1,1],
                     [-1,1],
                     [1,2],
                     [-1,2],
                     [2,1],
                     [-2,1],
                     [2,2],
                     [-2,2],
                     [3,1],
                     [-3,1],
                     [3,2],
                     [-3,2],
                     [4,1],
                     [-4,1],
                     [5,1],
                     [-5,1]
                     ])


# Specify relevant indices of MONdata
sm_power_monrange = (0, 19)
mm_power_monrange = (19, 38)
mm_phase_monrange = (38, 57)
sm_phase_monrange = (57, 76)


# Extract required data from rsoft files
all_smpowers = []
all_mmpowers = []
all_mmphases = []
all_smphases = []

for wgnum in range(nwgs):
    rsoft_suffix = 'LP%d%d' % (LP_modes[wgnum,0], LP_modes[wgnum,1])
    rsoft_filename = rsoft_fileprefix + rsoft_suffix
    print('Reading rsoft files ' + rsoft_filename)
    r = Rsoftdata(rsoft_datadir)
    r.readall(filename=rsoft_filename)
    if show_plots:
        r.plotall()
        plt.pause(0.001)

    smpower_mons = r.MONdata[:, sm_power_monrange[0]:sm_power_monrange[1]]
    mmpower_mons = r.MONdata[:, mm_power_monrange[0]:mm_power_monrange[1]]
    mmphase_mons = r.MONdata[:, mm_phase_monrange[0]:mm_phase_monrange[1]]
    smphase_mons = r.MONdata[:, sm_phase_monrange[0]:sm_phase_monrange[1]]

    if av_fluxes > 0:
        f = smpower_mons[-av_fluxes:, :]
        smpower = f.mean(axis=0)
        f = mmpower_mons[0:av_fluxes, :]
        mmpower = f.mean(axis=0)
    else:
        smpower = smpower_mons[-1, :]
        mmpower = mmpower_mons[0, :]
    mmphase = mmphase_mons[0, :]
    smphase = smphase_mons[-1, :]

    all_smpowers.append(smpower)
    all_mmpowers.append(mmpower)
    all_mmphases.append(mmphase)
    all_smphases.append(smphase)

    plt.figure(2)
    plt.clf()
    plt.plot(smpower_mons)
    plt.pause(0.001)

# outfilename = processed_datadir + 'extractedvals_' + rsoft_fileprefix + '.npz'
# np.savez(outfilename, all_smpowers=all_smpowers, all_mmpowers=all_mmpowers, all_mmphases=all_mmphases)

if show_plots:
    plt.figure(3)
    plt.clf()
    plt.subplot(211)
    plt.imshow(np.asarray(all_smpowers))
    plt.colorbar()
    plt.title('Output mode power')
    plt.ylabel('Excited waveguide no.')
    plt.xlabel('Mode no.')
    plt.subplot(212)
    plt.imshow(np.asarray(all_smphases), cmap='twilight_shifted')
    plt.colorbar()
    plt.title('Output mode phase')
    plt.ylabel('Excited waveguide no.')
    plt.xlabel('Mode no.')
    plt.tight_layout()


