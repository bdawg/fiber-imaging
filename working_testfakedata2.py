import matplotlib.pyplot as plt
import numpy as np
from lanternfiber import lanternfiber
plt.ion()

processed_datadir = '/Users/bnorris/DataAnalysis/fiber-imaging/'
rsoft_fileprefix = 'probeset_19LP_02probeset_19LP_02_'
data_filename = 'extractedvals_probeset_19LP__Good202107.npz'

f = lanternfiber(datadir=processed_datadir)
f.load_savedvalues(data_filename)
f.make_transfer_matrix_mm2sm()
# f.test_matrix(f.Cmat, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases, pausetime=0.1)

# ######### Test the matrix using some random test field propagations
# rsoft_datadir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix_19cPLJ21_MM2SM_makeInputFLDs_monobj/'
# rsoft_fileprefix = 'probeset_19LP_02randset5-anyampphaseprobeset_19LP_02randset5-anyampphase_'
# indfile = '19cPLJ21_mm2sm_run20210803_launchfile_MonObj_19mons.ind'
# num_test_files = 4
# f.load_rsoft_data_customfld(rsoft_datadir, rsoft_fileprefix, indfile,
#                             np_fileprefix='probeset_19LP_02randset5-anyampphase_',
#                             show_plots=False, save_output=False, nwgs=num_test_files, show_indivmasks=False,
#                             use_pathway_mons=False, use_monitor_objects=True, reorder_monobjs=True)
# f.test_matrix(f.Cmat, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases, pausetime=0.1,
#               num_to_show=num_test_files, unnormalise_smpower=True)

def unpack_cvec(input_vec):
    n_out = len(input_vec)*2
    out_vec = np.zeros(n_out)
    for k in range(n_out):
        inmode_num = np.int(k/2)
        if k%2 == 0:
            out_vec[k] = np.real(input_vec[inmode_num])
        else:
            out_vec[k] = np.imag(input_vec[inmode_num])
    return out_vec


# Make simulated data
ndata = 100000
amp_range = [0,1]
phase_range = [0, np.pi]
simdata_inputs, simdata_outputs = f.make_sim_data(ndata, amp_range=amp_range, phase_range=phase_range)

# Unpack input data into real,imag pairs
simdata_input_upk = []
for invec in simdata_inputs:
    invec_upk = unpack_cvec(invec)
    simdata_input_upk.append(invec_upk)

# Make output data into intensities
simdata_outputI = []
for outvec in simdata_outputs:
    # outvec_upk = unpack_cvec(outvec)
    simdata_outputI.append(np.abs(outvec)**2)

# Save to fiel as numpy arrays
outfilename = 'simdata_amp0-1_ph0-pi_01.npz'
simdata_input_arr = np.array(simdata_input_upk)
simdata_outputI_arr = np.array(simdata_outputI)
np.savez(outfilename, simdata_input_arr=simdata_input_arr, simdata_outputI_arr=simdata_outputI_arr)



