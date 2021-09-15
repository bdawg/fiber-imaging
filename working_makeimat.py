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





######### Test the matrix using some random test field propagations
rsoft_fileprefix = 'probeset_19LP_02randset5-anyampphaseprobeset_19LP_02randset5-anyampphase_'
num_test_files = 1
f.load_rsoft_data_customfld(rsoft_datadir, rsoft_fileprefix, indfile,
                            np_fileprefix='probeset_19LP_02randset5-anyampphase_',
                            show_plots=False, save_output=False, nwgs=num_test_files, show_indivmasks=False,
                            use_pathway_mons=False, use_monitor_objects=True, reorder_monobjs=True)
# f.test_matrix(f.Cmat, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases, pausetime=0.1,
#               num_to_show=num_test_files, unnormalise_smpower=True)


ind=0
in_power = f.all_mmpowers[ind]
in_phase = f.all_mmphases[ind]
output_power = f.all_smpowers[ind]
output_phases = f.all_smphases[ind]

### Test with forwards or backwards directions
M = f.Cmat
f.test_matrix(M, f.all_mmpowers, f.all_mmphases, f.all_smpowers, f.all_smphases, pausetime=0.1,
              num_to_show=num_test_files, unnormalise_smpower=True)
# M = np.asarray(np.matrix(f.Cmat).H)
# f.test_matrix(M, f.all_smpowers, f.all_smphases, f.all_mmpowers, f.all_mmphases, pausetime=0.1,
#               num_to_show=num_test_files, unnormalise_smpower=True)





# ######## Test using real-valued matrix
#
# Rmat = f.matrix_complex2real(f.Cmat)
#
# # Cmat = np.zeros((Imat.shape[0]//2, Imat.shape[1]//2), dtype='complex')
# # for m in range(Cmat.shape[0]):
# #     for n in range(Cmat.shape[1]):
# #         block = Imat[m*2:m*2+2, n*2:n*2+2]
# #         a = block[0,0]
# #         b = block[1,0]
# #         Cmat[m,n] = a + b*1j
#
# input_val_complex = in_power * np.exp(1j*in_phase/180*np.pi)
# input_val_complex = input_val_complex.reshape((-1,1))
# input_val_real = f.matrix_complex2real(input_val_complex)
#
# output_val_real = np.matmul(Rmat, input_val_real)
#
# output_val_complex = f.matrix_real2complex(output_val_real)
#
# pred_pow = np.abs(output_val_complex)
# newangs = np.angle(output_val_complex) /np.pi*180
# newangs[newangs < 0] = newangs[newangs < 0]+360
#
# output_power = output_power * np.sqrt(np.sum(in_power**2))
#
# plt.figure(2)
# plt.clf()
# plt.subplot(211)
# plt.plot(output_power,'-x')
# plt.plot(pred_pow,'-+')
# plt.subplot(212)
# plt.plot(output_phases,'-x')
# plt.plot(newangs,'-+')




########## Make output-intensity version
# Make a new set of columns of output intensities, each with unit excitation of each complex mode
nmodes = f.nmodes
Imat = np.zeros((f.nmodes, f.nmodes*2))
Cmat = f.Cmat
testmat = np.zeros((f.nmodes, f.nmodes*2), dtype=complex)
test_compImat = np.zeros((f.nmodes, f.nmodes*2), dtype=complex)
for k in range(f.nmodes*2):
    # For even columns, mode is 1+0j, and for odd columns it's 0+1j
    inmode_num = np.int(k/2)
    if k%2 == 0:
        inmode_val = 1+0j
    else:
        inmode_val = 0+1j
    inmode = np.zeros(f.nmodes, dtype=complex)
    inmode[inmode_num] = inmode_val
    testmat[:,k] = inmode

    output_ampl = Cmat @ inmode
    test_compImat[:, k] = output_ampl
    Imat[:, k] = np.abs(output_ampl)**2


# fake_Imat = np.zeros((nmodes*2,nmodes*2))
# fake_Imat[:nmodes,:] = Imat
# ncom = 3
# for k in range(nmodes):
#     newrow = np.zeros(nmodes*2)
#     for l in range(ncom):
#         ind = int(np.random.rand()*19)
#         newrow = Imat[ind,:] * np.random.rand()
#     # cfs = np.random.rand(nmodes)
#     # newrow = cfs @ Imat
#     fake_Imat[nmodes+k,:] = newrow
#
# plt.imshow(fake_Imat)
#
# Rmat = f.matrix_complex2real(Cmat)
# u, s, vh = np.linalg.svd(Rmat, full_matrices=True)
# plt.plot(s,'o-')



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


## Test propagating through the Imat

# Get zero point
modecoeff = np.ones(nmodes).astype(complex)
modecoeff = np.zeros(nmodes).astype(complex)
modecoeff[4] = 1
modecoeff[3] = 1
modecoeffI = unpack_cvec(modecoeff)
zp_fromcplx = np.abs(Cmat @ modecoeff)**2
zp_fromImat = Imat @ modecoeffI

plt.clf()
plt.plot(zp_fromcplx, '-x')
plt.plot(zp_fromImat, 'o')


# Make a random set of input coeffs
nmodes = f.nmodes
loval = 0.9
hival = 1.1
loval_phase = 0
hival_phase = 0#2*np.pi
mode_amps = np.random.rand(nmodes) * (hival-loval) + loval
mode_phases = np.random.rand(nmodes) * (hival_phase-loval_phase) + loval_phase
modecoeff = mode_amps * np.exp(1j*mode_phases)

# modecoeff = f.all_mmpowers[0] * np.exp(1j*f.all_mmphases[0]/180*np.pi)

pred_from_complex = Cmat @ modecoeff
pred_from_complex_intens = np.abs(pred_from_complex)**2 - zp_fromcplx

modecoeffI = unpack_cvec(modecoeff)
pred_I = Imat @ modecoeffI - zp_fromImat

plt.clf()
plt.plot(pred_from_complex_intens, 'r')
plt.plot(pred_I, 'b')
plt.tight_layout()



# modecoeffs = []
# for k in range(nmodes):
#     mode_amps = np.random.rand(nmodes) * (hival-loval) + loval
#     mode_phases = np.random.rand(nmodes) * (hival_phase-loval_phase) + loval_phase
#     modecoeff = mode_amps * np.exp(1j*mode_phases)
#     modecoeffs.append(modecoeff)
# modecoeffs = np.array(modecoeffs)






