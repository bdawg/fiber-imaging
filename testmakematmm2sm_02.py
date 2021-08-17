
from viewRSoftData import *
import matplotlib.pyplot as plt
import numpy as np
from lanternfiber import lanternfiber
plt.ion()


rsoft_datadir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix_19cPLJ21_MM2SM_run20210803/'
rsoft_fileprefix = '19cPLJ21_mm2sm_run20210803_'
processed_datadir = '/Users/bnorris/DataAnalysis/fiber-imaging/'
nwgs = 19 #19
show_plots = True
av_fluxes = 100 #100 # If >0, average n flux measurements from monitor
data_filename = 'extractedvals_19cPLJ21_mm2sm_run20210803_.npz'

f = lanternfiber(datadir=processed_datadir)
f.load_rsoft_data_mm2sm(rsoft_datadir, rsoft_fileprefix, show_plots=show_plots, save_output=True)
# f.load_savedvalues(data_filename)

f.set_mmvals_nominal(square=True)
# f.make_transfer_matrix_mm2sm(mm_phase=None)
f.make_transfer_matrix_mm2sm(mm_phase=0)
f.test_matrix(f.Dmat, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases, pausetime=0.1)


rsoft_fileprefix = '19cPLJ21_mm2sm_run20210803_randlaunches_nophase_powers1_01'
# rsoft_fileprefix = '19cPLJ21_mm2sm_run20210803_randlaunch03'
# rsoft_fileprefix = 'rough_19cPLJ21_mm2sm_run20210803_randlaunches_nophase_01_onlyLP03LP21'
# rsoft_fileprefix = 'midres_19cPLJ21_mm2sm_run20210803_randlaunches_nophase_01_severalPosLP_nomons'
# rsoft_fileprefix = 'midres_19cPLJ21_mm2sm_run20210803_randlaunches_nophase_01_LP41only_nomons'

f.load_rsoft_data_mm2sm_single(rsoft_datadir, rsoft_fileprefix, show_plots=False)
# f.all_mmphases[0]=np.zeros(19)
f.test_matrix(f.Dmat, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases,
              pausetime=0.1, fignum=2)


f.load_rsoft_data_fldonly(rsoft_datadir, rsoft_fileprefix)
ind_filename = '19cPLJ21_mm2sm_run20210803_randlaunches_nophase_powers1_01.ind'
# ind_filename = 'midres_19cPLJ21_mm2sm_run20210803_randlaunches_nophase_01_LP41only_nomons.ind'
# ind_filename = 'rough_19cPLJ21_mm2sm_run20210803_randlaunches_nophase_01_onlyLP01.ind'
all_mmpowers, all_mmphases, all_wgposns = f.read_ind_file(rsoft_datadir, ind_filename, skipfirst=True, getWGposns=True)

# # all_mmpowers = np.ones(19)
# all_mmpowers = all_mmpowers / np.sqrt(np.sum(np.array(all_mmpowers)**2))
# f.all_mmpowers[0] = all_mmpowers
# f.all_mmphases[0] = np.array(all_mmphases)
# f.all_mmphases[0]=np.zeros(19)

# f.show_outfield()

f.measure_wg_fields(fignum=1, show_plots=False)
f.normalise_smpowers()
f.test_matrix(f.Dmat, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases,
              pausetime=0.1, fignum=2)




# plt.clf()
# plt.imshow(f.out_field_ampl[0], cmap='jet')
# plt.plot(all_wgposns[0,:], all_wgposns[1,:], 'x')

# ind_filename = '19cPLJ21_mm2sm_run20210803_LP02.ind'
# f.all_mmpowers, f.all_mmphases = f.read_ind_file(rsoft_datadir, ind_filename)
# all_mmpowers, all_mmphases = f.read_ind_file(rsoft_datadir, ind_filename, skipfirst=True)

# o_amps = f.all_smpowers
# o_phases = f.all_smphases
# f.measure_wg_fields()
# plt.plot(f.all_smpowers[0] / np.sum(f.all_smpowers[0]))
# plt.plot(o_amps[0] / np.sum(o_amps[0]))
# plt.plot(f.all_smphases[0])
# plt.plot(o_phases[0])

# ap_rad = 5
# ampl_im = f.out_field_ampl[0]
# phase_im = f.out_field_phase[0]
# wg_posns = f.wg_posns
# imsz = ampl_im.shape[0]
# Y, X = np.ogrid[:imsz, :imsz]
# for k in range(1):
#     dist = np.sqrt((X-f.wg_posns[0,k])**2 + (Y-f.wg_posns[1,k])**2)
#     mask = dist <= ap_rad
#     im = np.copy(phase_im)
#     im[~mask] = 0
#     plt.imshow(im)














