# Make a set of launch fields and batch files to run through beamprop


from lanternfiber import lanternfiber
import numpy as np

output_dir = '/Users/bnorris/temp/testing/'
# output_dir = '/Users/bnorris/Dropbox/Win-Mac Share/rsoft/PL_roundtriptest/'

### Fibre parameters
n_core = 1.44
n_cladding = 1.4345
wavelength = 1.5 # microns
core_radius = 32.8/2 # microns

### Scale parameters
max_r = 2 # Maximum radius to calculate mode field, where r=1 is the core diameter
npix = 200 # Half-width of mode field calculation in pixels

show_plots = False
indfile = '19cPLJ21_mm2sm_run20210803_launchfile_MonObj_19mons.ind'


f = lanternfiber(n_core, n_cladding, core_radius, wavelength)
# Make a set of probe fields:
output_filename = 'probeset_19LP_run02'
f.make_rsoft_launch_fields('probe', npix=npix, max_r=max_r, indfile=indfile, show_plots=show_plots,
                           make_bat_file=True, num_bat_files=3, outpath=output_dir,
                           outprefix=output_filename)

# # Make a set of 4 random test fields:
# output_filename = 'randomset_19LP_run02'
# f.make_rsoft_launch_fields('randampphase', npix=npix, indfile=indfile, show_plots=show_plots,
#                            make_bat_file=True, num_bat_files=3, outpath=output_dir,
#                            outprefix=output_filename, num_outs=4)