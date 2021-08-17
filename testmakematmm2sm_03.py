
from viewRSoftData import *
import matplotlib.pyplot as plt
import numpy as np
from lanternfiber import lanternfiber
plt.ion()

rsoft_datadir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix_19cPLJ21_MM2SM_makeInputFLDs/'
rsoft_fileprefix = 'probeset_19LP_01probeset_19LP_01_'
indfile = '19cPLJ21_mm2sm_run20210803_launchfile_fewmons.ind'
processed_datadir = '/Users/bnorris/DataAnalysis/fiber-imaging/'
nwgs = 19 #19
show_plots = True
av_fluxes = 100 #100 # If >0, average n flux measurements from monitor
data_filename = 'extractedvals_probeset_19LP_01probeset_19LP_01_.npz'

indfile = '19cPLJ21_mm2sm_run20210803_launchfile_SMmons.ind'
rsoft_fileprefix = 'probeset_19LP_02probeset_19LP_02_'
data_filename = 'extractedvals_probeset_19LP_02_fromPathwayMons.npz'

rsoft_datadir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix_19cPLJ21_MM2SM_makeInputFLDs_monobj/'
indfile = '19cPLJ21_mm2sm_run20210803_launchfile_MonObj_19mons.ind'
rsoft_fileprefix = 'probeset_19LP_02probeset_19LP_02_'
data_filename = ''

f = lanternfiber(datadir=processed_datadir)
f.load_rsoft_data_customfld(rsoft_datadir, rsoft_fileprefix, indfile, np_fileprefix='probeset_19LP_02_',
                            show_plots=show_plots, save_output=True, show_indivmasks=False,
                            use_pathway_mons=False, use_monitor_objects=True)
# f.load_savedvalues(data_filename)

# f.set_mmvals_nominal(square=True)
f.make_transfer_matrix_mm2sm(mm_phase=None)

f.test_matrix(f.Dmat, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases, pausetime=0.1)



# rsoft_fileprefix = 'probeset_19LP_01randset3probeset_19LP_01randset3_'
# indfile = '19cPLJ21_mm2sm_run20210803_launchfile_fewmons.ind'

# indfile = '19cPLJ21_mm2sm_run20210803_launchfile_SMmons.ind'
rsoft_fileprefix = 'probeset_19LP_02randset5-anyampphaseprobeset_19LP_02randset5-anyampphase_'


f.load_rsoft_data_customfld(rsoft_datadir, rsoft_fileprefix, indfile, np_fileprefix='probeset_19LP_02randset5-anyampphase_',
                            show_plots=show_plots, save_output=False, nwgs=4, show_indivmasks=False,
                            use_pathway_mons=False, use_monitor_objects=True)
f.test_matrix(f.Dmat, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases, pausetime=0.1,
              num_to_show=4, unnormalise_smpower=True)


