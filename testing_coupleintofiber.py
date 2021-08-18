from lanternfiber import lanternfiber
import numpy as np
import matplotlib.pyplot as plt




n_core = 1.44
n_cladding = 1.4345
wavelength = 1.5 # microns
core_radius = 32.8/2 # microns

### Scale parameters
max_r = 2 # Maximum radius to calculate mode field, where r=1 is the core diameter
npix = 200 # Half-width of mode field calculation in pixels

show_plots = True


f = lanternfiber(n_core, n_cladding, core_radius, wavelength)
f.find_fiber_modes()
f.make_fiber_modes(npix=npix, show_plots=False, max_r=max_r)
f.plot_fiber_modes(0)


### Define input field

# Small Gaussian point source
power = 0.5
posn_microns = [0, 0]
sigma = 10

f.make_arb_input_field('gaussian', power=power, location=posn_microns, sigma=sigma,
                       add_to_existing=False, show_plots=show_plots)

# power = 0.5
# posn_microns = [13, 0]
# f.make_arb_input_field('gaussian', power=power, location=posn_microns, sigma=sigma,
#                        add_to_existing=True, show_plots=show_plots, phase=np.pi/2)


mode_num = 1
f.plot_injection_field(f.input_field)
coupling = f.calc_injection(mode_field_number=1, verbose=True)



# input_field = f.input_field
# mode_field = f.make_complex_fld(f.allmodefields_cos_cart[1])
#
# f.plot_injection_field(input_field)
# plt.pause(1)
# f.plot_injection_field(mode_field)
#
# # plt.imshow(np.abs(input_field*mode_field))
#
# overlap_int = np.abs(np.sum(input_field*mode_field))**2 / \
#               ( np.sum(np.abs(mode_field)**2) * np.sum(np.abs(input_field)**2) )
# # overlap_int = np.abs(np.sum(input_field*mode_field)**2 / np.sum(mode_field**2))**2
# print(overlap_int)

# o = input_field*mode_field
# f.plot_injection_field(o)
# print(np.sum(np.abs(o)**2))
#
# plt.imshow(np.abs(o)**2)


