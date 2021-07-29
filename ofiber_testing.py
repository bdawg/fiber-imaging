import numpy as np
import ofiber
import matplotlib.pyplot as plt
import polarTransform
plt.ion()

n_core = 1.43
n_cladding = 1.425
wavelength = 0.5 # microns
core_radius = 5 # microns

NA = ofiber.numerical_aperture(n_core, n_cladding)
V = ofiber.V_parameter(core_radius, NA, wavelength)

# Find all modes for this fiber
max_l = 100 # Maximum number of l modes to find (can be arbitrarily large)
all_modes_b = []
all_l = []
all_m = []
for l in range(max_l):
    cur_b = ofiber.LP_mode_values(V, l)
    if len(cur_b) == 0:
        break
    else:
        all_modes_b.extend(cur_b)
        ls = (np.ones_like(cur_b)) * l
        all_l.extend(ls.astype(int))
        ms = np.arange(len(cur_b))+1
        all_m.extend(ms)

all_modes_b = np.asarray(all_modes_b)
nmodes = len(all_modes_b)
print('Total number of modes found: %d' % nmodes)


# Plot a mode field
# mode_to_plot = 3
max_r = 2
numpix = 100
zlim = 0.04

r = np.linspace(0, max_r, numpix) # Radial positions, normalised so core_radius = 1

for mode_to_plot in range(nmodes):
    field_1d = ofiber.LP_radial_field(V, all_modes_b[mode_to_plot], all_l[mode_to_plot], r)

    phivals = np.linspace(0, 2*np.pi, numpix)
    phi_cos = np.cos(all_l[mode_to_plot] * phivals)
    phi_sin = np.sin(all_l[mode_to_plot] * phivals)

    rgrid, phigrid = np.meshgrid(r, phivals)
    field_r_cos, field_phi = np.meshgrid(phi_cos, field_1d)
    field_r_sin, field_phi = np.meshgrid(phi_sin, field_1d)
    field_cos = field_r_cos * field_phi
    field_sin = field_r_sin * field_phi

    # Normalise each field so its total intensity is 1
    field_cos = field_cos / np.sqrt(np.sum(field_cos**2))
    field_sin = field_sin / np.sqrt(np.sum(field_sin**2))
    field_cos = np.nan_to_num(field_cos)
    field_sin = np.nan_to_num(field_sin)

    field_cos_cart, d = polarTransform.convertToCartesianImage(field_cos.T)
    field_sin_cart, d = polarTransform.convertToCartesianImage(field_sin.T)

    plt.clf()
    plt.subplot(121)
    sz = max_r * core_radius
    plt.imshow(field_cos_cart, extent=(-sz, sz, -sz, sz), cmap='bwr', vmin=-zlim, vmax=zlim)
    plt.xlabel('Position ($\mu$m)')
    plt.ylabel('Position ($\mu$m)')
    plt.title('Mode l=%d, m=%d' % (all_l[mode_to_plot], all_m[mode_to_plot]))
    core_circle = plt.Circle((0,0), core_radius, color='k', fill=False, linestyle='--', alpha=0.2)
    plt.gca().add_patch(core_circle)
    plt.subplot(122)
    sz = max_r * core_radius
    plt.imshow(field_sin_cart, extent=(-sz, sz, -sz, sz), cmap='bwr', vmin=-zlim, vmax=zlim)
    plt.xlabel('Position ($\mu$m)')
    plt.title('Mode l=%d, m=%d' % (all_l[mode_to_plot], all_m[mode_to_plot]))
    core_circle = plt.Circle((0,0), core_radius, color='k', fill=False, linestyle='--', alpha=0.2)
    plt.gca().add_patch(core_circle)

    plt.pause(1)

