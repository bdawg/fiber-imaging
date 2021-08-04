
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
# f.load_rsoft_data_mm2sm(rsoft_datadir, rsoft_fileprefix, show_plots=show_plots, save_output=True)
f.load_savedvalues(data_filename)

f.make_transfer_matrix_mm2sm(mm_phase=0)
# f.test_matrix(f.Dmat, f.all_mmpowers, f.all_smpowers, f.all_smphases, pausetime=1)


# rsoft_fileprefix = '19cPLJ21_mm2sm_run20210803_randlaunch01'
# f.load_rsoft_data_single(rsoft_datadir, rsoft_fileprefix, show_plots=True)


ind_filename = '19cPLJ21_mm2sm_run20210803_LP02.ind'
# f.all_mmpowers, f.all_mmphases = f.read_ind_file(rsoft_datadir, ind_filename)
all_mmpowers, all_mmphases = f.read_ind_file(rsoft_datadir, ind_filename, skipfirst=False)