from lanternfiber import lanternfiber
# from make_rsoft_fld_batch import RSOFTPSF
import numpy as np
import matplotlib.pyplot as plt

output_dir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_getmatrix_19cPLJ21_MM2SM_run20210803/'
output_filename = 'testlaunch_LP21_02'


### Fibre parameters
n_core = 1.44
n_cladding = 1.4345
wavelength = 1.5 # microns
core_radius = 32.8/2 # microns

### Scale paramters
max_r = 2 # Maximum radius to calculate mode field, where r=1 is the core diameter
npix = 100 # Half-width of mode field calculation in pixels

mode_of_interest = 5
show_plots = False

array_size_microns = max_r * core_radius * 2
microns_per_pixel = array_size_microns / (npix*2)

f = lanternfiber(n_core, n_cladding, core_radius, wavelength)
f.find_fiber_modes()
f.make_fiber_modes(show_plots=show_plots, npix=npix)
f.plot_fiber_modes(mode_of_interest)


raw_ampl = f.allmodefields_cos_cart[mode_of_interest]
fld_ampl = np.abs(raw_ampl)
fld_phase = np.zeros_like(raw_ampl)
fld_phase[raw_ampl < 0] = np.pi
fld_complex = fld_ampl * np.exp(1j * fld_phase)
fld_complex = np.transpose(fld_complex)

# r = RSOFTPSF()
# r.complex_psf = fld_complex
# r.saveToRSoft(outfile=output_dir+output_filename, size_data=array_size_microns/2)
# f.saveToRSoft(fld_complex, outfile=output_dir+output_filename, size_data=array_size_microns/2)

# r.makeMultiplePSFs()
