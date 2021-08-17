from lanternfiber import lanternfiber
import numpy as np
import matplotlib.pyplot as plt



n_core = 1.43
n_cladding = 1.425
wavelength = 0.5 # microns
core_radius = 5 # microns

n_core = 1.44
n_cladding = 1.4345
wavelength = 1.55 # microns
core_radius = 32.8/2 # microns


f = lanternfiber(n_core, n_cladding, core_radius, wavelength)
f.find_fiber_modes()
f.make_fiber_modes(show_plots=True)
# f.plot_fiber_modes(0)




