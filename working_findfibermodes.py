from lanternfiber import lanternfiber
import numpy as np
import matplotlib.pyplot as plt



# n_core = 1.43
# n_cladding = 1.425
# wavelength = 0.5 # microns
# core_radius = 5 # microns

n_core = 1.44
n_cladding = 1.4345
wavelength = 1.5 # microns
core_radius = 32.8/2 # microns
core_radius = 6.5/2 # microns

core_radius = 20.2/2 # microns # 8 modes
# wavelength = 2.45

f = lanternfiber(n_core, n_cladding, core_radius, wavelength)
# print(f.V)
f.find_fiber_modes()
f.make_fiber_modes(show_plots=True)
# f.plot_fiber_modes(0)

# core_radii = np.linspace(15, 40, 100) / 2
# n_unique = []
# for core_radius in core_radii:
#     f = lanternfiber(n_core, n_cladding, core_radius, wavelength)
#     nu = f.find_fiber_modes(return_n_unique=True)
#     n_unique.append(nu)
#
# plt.clf()
# plt.plot(core_radii, n_unique)


# wl = np.linspace(1., 3, 100)
# n_unique = []
# for wavelength in wl:
#     f = lanternfiber(n_core, n_cladding, core_radius, wavelength)
#     nu = f.find_fiber_modes(return_n_unique=True)
#     n_unique.append(nu)
#
# plt.plot(wl, n_unique)



