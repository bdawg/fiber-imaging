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
# f.plot_fiber_modes(0)


### Define input field

# Small Gaussian point source
power = 1
posn_microns = [0, 0]
sigma = 3

f.make_arb_input_field('gaussian', power=power, location=posn_microns, sigma=sigma,
                       add_to_existing=False, show_plots=show_plots)

power = 2
posn_microns = [4, 10]
sigma = 6
f.make_arb_input_field('gaussian', power=power, location=posn_microns, sigma=sigma,
                       add_to_existing=True, show_plots=show_plots)

f.plot_injection_field(f.input_field)