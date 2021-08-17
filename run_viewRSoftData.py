#%%
from viewRSoftData import *

# plt.ion()

# datapath = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/19CorePL/'
# datapath = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/19CorePL/run19_ZernSet01/'
# datapath = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/19CorePL/run19_ZernFoc-0.2/'
# datapath = './rsoftOutputs/'
# scanset = ""

# datapath = "/Users/bnorris/Dropbox/Win-Mac Share/rsoft/19CorePL/scan_focusScan01/"
# scanset = "focusScan01"
# datapath = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/19CorePL/scan_comaScan01/'
# scanset = 'comaScan01'
# datapath = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/19CorePL/scan_tipScan01/'
# scanset = 'tipScan01'
# datapath = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/19CorePL/scan_astigScan01/'
# scanset = 'astigScan01'
# datapath = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/19CorePL/scan_sphericalScan01/'
# scanset = 'sphericalScan01'

# datapath = "/Volumes/silo4/snert/FMF_PL_rsoft/best_gridsize/"

# datapath = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/19CorePL/scan_astigScan02_hex/'
# scanset = 'astigScan01'

# datapath = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/19CorePL/scan_foctest02/'
datapath = "/Volumes/silo4/snert/FMF_PL_rsoft/zernike/"
filepref = 'zernike'

datapath = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix/'
filepref = '19cPLJ21_sm2mm_noOS_wg1_grid0-5'

datapath = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix_19cPLJ21_MM2SM_run20210803/'
filepref = '19cPLJ21_mm2sm_run20210803_randlaunches_nophase_01'

# filepref = 'bptmp4'

# datapath = '/Users/bnorris/Dropbox/Win-Mac Share/fiber_sims/'
# filepref = 'bptmp2'

r = Rsoftdata(datapath)
r.readall(filename=filepref)
# r.readall('focusScan02zernikePSFs_0.0000_0.0000_0.0000_-0.0500_realhex')

# r.readMulti(filepref, keepAllResults=True)

# r.plotall(clim=[-200, 200])
r.plotall(clim=[0,0.3])

# r.readMulti(scanset, keepAllResults=False, readOne=10)

# print('Final fluxes:')
# print(r.finalFluxes())

plt.figure()
ff = r.finalFluxes()
plt.plot(np.log10(ff[19:]), 'x')

# for k in range(r.MONdata.shape[1]):
#     plt.plot(r.MONdata[:, k])
#     plt.ylim([0, 1])
#     print(k)
#     plt.pause(0.5)

# # Test effect of averaging final points in monitors
# plt.figure()
# fluxes = []
# for k in range(1, 4000, 10):
#     fluxes.append(r.finalFluxes(useNPts=k))
# p=plt.plot(fluxes)

# np.savez('testsave.npz', allResults=r.allResults)


# %%
