import matplotlib.pyplot as plt
import numpy as np
from lanternfiber import lanternfiber
plt.ion()

processed_datadir = '/Users/bnorris/DataAnalysis/fiber-imaging/'
rsoft_datadir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix_19cPLJ21_MM2SM_makeInputFLDs_monobj/'
indfile = '19cPLJ21_mm2sm_run20210803_launchfile_MonObj_19mons.ind'



######## Calculate a matrix using the probe fields
rsoft_fileprefix = 'probeset_19LP_02probeset_19LP_02_'
data_filename = ''

f = lanternfiber(datadir=processed_datadir)
f.load_rsoft_data_customfld(rsoft_datadir, rsoft_fileprefix, indfile, np_fileprefix='probeset_19LP_02_',
                            show_plots=False, save_output=False, show_indivmasks=False,
                            use_pathway_mons=False, use_monitor_objects=True, reorder_monobjs=True)
f.make_transfer_matrix_mm2sm()
# f.test_matrix(f.Cmat, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases, pausetime=0.1)




# ######### Test the matrix using some random test field propagations
# rsoft_fileprefix = 'probeset_19LP_02randset5-anyampphaseprobeset_19LP_02randset5-anyampphase_'
# num_test_files = 1
# f.load_rsoft_data_customfld(rsoft_datadir, rsoft_fileprefix, indfile,
#                             np_fileprefix='probeset_19LP_02randset5-anyampphase_',
#                             show_plots=False, save_output=False, nwgs=num_test_files, show_indivmasks=False,
#                             use_pathway_mons=False, use_monitor_objects=True, reorder_monobjs=True)
# # f.test_matrix(f.Cmat, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases, pausetime=0.1,
# #               num_to_show=num_test_files, unnormalise_smpower=True)


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



nmodes = f.nmodes
Cmat = f.Cmat

# Linear approximation of intensity will be centred around a zero point, formed when all modes, all with coeffs 1+0j,
# are injected.
incoeffs_0 = np.ones(nmodes, dtype=complex)
outcoeffs_0 = Cmat @ incoeffs_0
outintens_0 = np.abs(outcoeffs_0)**2

# Make I matrix mapping delta input complex coeffs to delta output intensities
Imat = np.zeros((nmodes, nmodes*2))
for k in range(nmodes*2):
    # For even columns, input mode is 1+0j, and for odd columns it's 0+1j
    inmode_num = np.int(k/2)
    if k%2 == 0:
        inmode_val = 1+0j
    else:
        inmode_val = 0+1j
    incoeffs = np.zeros(f.nmodes, dtype=complex)
    incoeffs[inmode_num] = inmode_val

    d_incoeffs = incoeffs - incoeffs_0
    # d_outintens = np.abs(Cmat @ d_incoeffs)**2 / np.sum(np.abs(d_incoeffs)**2)
    # d_outintens = np.abs(Cmat @ d_incoeffs / np.sum(d_incoeffs))**2
    d_outintens = np.abs(Cmat @ d_incoeffs)**2 / np.sum(np.abs(d_incoeffs))**2
    Imat[:, k] = d_outintens


# Find an output intensity for given input coeffs
incoeffs = np.ones(nmodes, dtype=complex)
incoeffs[2] = 3 + 0j

# loval = 0.9
# hival = 1.1
# loval_phase = 0
# hival_phase = 0.2
# mode_amps = np.random.rand(nmodes) * (hival-loval) + loval
# mode_phases = np.random.rand(nmodes) * (hival_phase-loval_phase) + loval_phase
# incoeffs = mode_amps * np.exp(1j*mode_phases)

d_incoeffs = incoeffs - incoeffs_0
d_incoeffs_unpacked = unpack_cvec(d_incoeffs)
pred_d_outintens = Imat @ d_incoeffs_unpacked
pred_outintens = pred_d_outintens + outintens_0

pred_outintens_cx = np.abs(Cmat @ incoeffs)**2
pred_d_outintens_cx = pred_outintens_cx - outintens_0

plt.clf()
# plt.plot(pred_outintens_cx, '-x', label='From Cmat')
# plt.plot(pred_outintens, '-+', label='From Imat')
# plt.plot(pred_d_outintens_cx, '-x', label='From Cmat')
plt.plot(pred_d_outintens, '-+', label='From Imat')
plt.legend()


# Run in a loop to test
ntests = 1000
all_pred_outintens_cx = []
all_pred_outintens = []
for k in range(ntests):
    loval = 0.9
    hival = 1.1
    loval_phase = 0
    hival_phase = 0.2
    mode_amps = np.random.rand(nmodes) * (hival-loval) + loval
    mode_phases = np.random.rand(nmodes) * (hival_phase-loval_phase) + loval_phase
    incoeffs = mode_amps * np.exp(1j*mode_phases)

    d_incoeffs = incoeffs - incoeffs_0
    d_incoeffs_unpacked = unpack_cvec(d_incoeffs)
    pred_d_outintens = Imat @ d_incoeffs_unpacked
    pred_outintens = pred_d_outintens + outintens_0
    pred_outintens_cx = np.abs(Cmat @ incoeffs)**2
    pred_d_outintens_cx = pred_outintens_cx - outintens_0

    all_pred_outintens_cx.append(pred_outintens_cx)
    all_pred_outintens.append(pred_outintens)

    # plt.clf()
    # plt.plot(pred_outintens_cx, '-x', label='From Cmat')
    # plt.plot(pred_outintens, '-+', label='From Imat')
    # plt.pause(0.1)

all_pred_outintens_cx = np.array(all_pred_outintens_cx)
all_pred_outintens = np.array(all_pred_outintens)



# plt.plot(all_pred_outintens_cx.ravel(), all_pred_outintens.ravel(), '.', markersize=0.5)


















