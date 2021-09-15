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


ndata = 10000
amp_range = [0,0.5]
phase_range = [0, 0]
# simdata_inputs, simdata_outputs = f.make_sim_data(ndata)
simdata_inputs, simdata_outputs = f.make_sim_data(ndata, amp_range=amp_range, phase_range=phase_range)

all_mmpowers, all_mmphases = f.cx2ap(simdata_inputs)
all_smpowers, all_smphases = f.cx2ap(simdata_outputs)
# f.test_matrix(f.Cmat, all_mmpowers, all_mmphases, all_smpowers, all_smphases, pausetime=0.1,
#               num_to_show=1, unnormalise_smpower=False, unwrap_phase=False)

plt.clf()
# plt.plot(np.sum(all_mmpowers,axis=1), np.sum(all_smpowers,axis=1), '.', markersize=2)
plt.plot(all_mmpowers[:,0], all_smpowers[:,1], '.', markersize=2)
# plt.plot(all_mmpowers.ravel(), all_smpowers.ravel()**2, '.', markersize=2)

nmodes = f.nmodes

npts = 100
invals = np.linspace(0,1,npts)
out_coeffs = []
in_coeffs = []
plt.clf()
for k in range(npts):
    cur_amps = np.ones(nmodes)
    cur_phases = np.zeros(nmodes)
    cur_phases[1] = invals[k]
    input_modecoeffs, output_smvals = f.make_sim_data(in_amp_phase=[cur_amps, cur_phases])
    in_coeffs.append(input_modecoeffs)
    out_coeffs.append(output_smvals)
out_coeffs = np.array(out_coeffs)
in_coeffs = np.array(in_coeffs)

plt.clf()
outvals = np.abs(out_coeffs)**2
plt.plot(outvals)
# outvals = outvals - np.mean(outvals,axis=0)
# plt.imshow(outvals, aspect='auto', cmap='gray')

# plt.clf()
# for k in range(f.nmodes):
#     plt.clf()
#     plt.plot(all_mmpowers[:,0], all_smpowers[:,k], '.', markersize=2)
#     plt.pause(1)

# Rmat = f.matrix_complex2real(f.Cmat)