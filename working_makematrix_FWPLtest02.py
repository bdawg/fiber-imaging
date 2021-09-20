"""
Make a transmission matrix going from the MM to SM end, from beamprop runs.
Assumes fields (and .npz) files have been generated )ie by lanternfiber..make_rsoft_launch_fields()
and beamprop has been run on them.

This example uses the preferred method of injecting modes into beamprop and measuring teh outputs:
- Input a single .fld file containing the total mode field
- Read SM output amps & phases using Monitor Objects (not pathway monitors, or encircled fluxes)
"""


import matplotlib.pyplot as plt
import numpy as np
from lanternfiber import lanternfiber
plt.ion()

processed_datadir = '/Users/bnorris/DataAnalysis/fiber-imaging/'
rsoft_datadir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/FWPL_test02/'
indfile = 'FWtest02_19cPLJ21.ind'
show_plots = True


### Calculate a matrix using the probe fields
rsoft_fileprefix = 'bp_probeset_FWtest02_run01_'
data_filename = ''

f = lanternfiber(datadir=processed_datadir)
f.load_rsoft_data_customfld(rsoft_datadir, rsoft_fileprefix, indfile, np_fileprefix='probeset_FWtest02_run01_',
                            show_plots=show_plots, save_output=True, show_indivmasks=False,
                            use_pathway_mons=False, use_monitor_objects=True, reorder_monobjs=True)
# f.load_savedvalues(data_filename) # Use instead of load_rsoft_data_* in subequent runs, if save_output was True
# f.make_transfer_matrix_mm2sm()
# f.test_matrix(f.Cmat, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases, pausetime=0.1)




# ### Test the matrix using some random test field propagations
# rsoft_fileprefix = 'probeset_19LP_02randset5-anyampphaseprobeset_19LP_02randset5-anyampphase_'
# num_test_files = 4
# f.load_rsoft_data_customfld(rsoft_datadir, rsoft_fileprefix, indfile, np_fileprefix='probeset_19LP_02randset5-anyampphase_',
#                             show_plots=False, save_output=False, nwgs=num_test_files, show_indivmasks=False,
#                             use_pathway_mons=False, use_monitor_objects=True, reorder_monobjs=True)
# f.test_matrix(f.Cmat, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases, pausetime=0.1,
#               num_to_show=num_test_files, unnormalise_smpower=True)

