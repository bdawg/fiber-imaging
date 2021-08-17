from lanternfiber import lanternfiber
import numpy as np
import matplotlib.pyplot as plt

output_dir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix_19cPLJ21_MM2SM_makeInputFLDs/'
output_filename = 'probeset_19LP_02'

### Fibre parameters
n_core = 1.44
n_cladding = 1.4345
wavelength = 1.5 # microns
core_radius = 32.8/2 # microns

### Scale paramters
max_r = 2 # Maximum radius to calculate mode field, where r=1 is the core diameter
npix = 200 # Half-width of mode field calculation in pixels

array_size_microns = max_r * core_radius * 2
microns_per_pixel = array_size_microns / (npix*2)

show_plots = False

# indfile = '19cPLJ21_mm2sm_run20210803_launchfile_fewmons.ind'
indfile = '19cPLJ21_mm2sm_run20210803_launchfile_SMmons.ind'

f = lanternfiber(n_core, n_cladding, core_radius, wavelength)
f.find_fiber_modes()
f.make_fiber_modes(show_plots=show_plots, npix=npix)

# Make one of each mode
modecoeffs = []
for k in range(f.nmodes):
    coeff = np.zeros(f.nmodes)
    coeff[k] = 1
    modecoeffs.append(coeff)

f.save_multiple_rsoft(modecoeffs, outpath=output_dir, outPrefix=output_filename, size_data=array_size_microns/2,
                      indFile=indfile, makeBatFile=True, numBatfiles=4)



# # Make some random modes
# nrandsets = 4
# loval = 0.5
# hival = 1
# modecoeffs = []
# for k in range(nrandsets):
#     modecoeffs.append(np.random.rand(f.nmodes) * (hival-loval) + loval)
# output_filename = output_filename+'randset2'
# # f.save_multiple_rsoft(modecoeffs, outpath=output_dir, outPrefix=output_filename, size_data=array_size_microns/2,
# #                       indFile=indfile)
# f.save_multiple_rsoft(modecoeffs, outpath=output_dir, outPrefix=output_filename, size_data=array_size_microns/2,
#                       indFile=indfile, makeBatFile=True, numBatfiles=4)



# Make some random modes with random phases
nrandsets = 10
loval = 0.5
hival = 1
loval_phase = 0
hival_phase = 2*np.pi
modecoeffs = []
for k in range(nrandsets):
    mode_amps = np.random.rand(f.nmodes) * (hival-loval) + loval
    mode_phases = np.random.rand(f.nmodes) * (hival_phase-loval_phase) + loval_phase
    modecoeff = mode_amps * np.exp(1j*mode_phases)
    modecoeffs.append(modecoeff)
output_filename = output_filename+'randset5-limitedampanyphase'
f.save_multiple_rsoft(modecoeffs, outpath=output_dir, outPrefix=output_filename, size_data=array_size_microns/2,
                      indFile=indfile, makeBatFile=True, numBatfiles=3)






