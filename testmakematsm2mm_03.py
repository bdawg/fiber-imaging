
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
# f.load_rsoft_data_sm2mm(rsoft_datadir, rsoft_fileprefix, show_plots=show_plots, zero_phases=True, save_output=True)
f.load_savedvalues(data_filename)

f.make_transfer_matrix_sm2mm()
f.test_matrix(f.Cmat, f.all_smpowers, f.all_smphases, f.all_mmpowers, f.all_mmphases, pausetime=0.001)

Cmat = f.Cmat
rsoft_fileprefix = '19cPLJ21_sm2mm_randlaunches02'
# rsoft_fileprefix = '19cPLJ21_sm2mm_srun20210802_wg02'
f.load_rsoft_data_sm2mm_single(rsoft_datadir, rsoft_fileprefix, show_plots=True)

ind_filename = '19cPLJ21_sm2mm_randlaunches02.ind'
old_smpowers = f.all_smpowers
old_smphases = f.all_smphases
# f.all_smpowers, f.all_smphases = f.read_ind_file(rsoft_datadir, ind_filename)
dummy, f.all_smphases = f.read_ind_file(rsoft_datadir, ind_filename)

f.test_matrix(f.Cmat, f.all_smpowers, f.all_smphases, f.all_mmpowers, f.all_mmphases, pausetime=0.1)


# rsoft_datadir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix_19cPLJ21_MM2SM_run20210803/'
# ind_filename = '19CorePL_July2021_noOuter_MMtoSM_extraMM_randlaunches01.ind'
# f.all_mmpowers, f.all_mmphases = f.read_ind_file(rsoft_datadir, ind_filename)

