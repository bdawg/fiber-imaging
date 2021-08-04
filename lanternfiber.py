"""
A class to do useful fiber and photonic lantern things, such as mode finding, coupling, etc.
"""
import numpy as np
import ofiber
import matplotlib.pyplot as plt
import polarTransform
from viewRSoftData import *

class lanternfiber:
    def __init__(self, n_core=None, n_cladding=None, core_radius=None, wavelength=None, nmodes=19, datadir='./'):
        self.n_core = n_core
        self.n_cladding = n_cladding
        self.core_radius = core_radius
        self.wavelength = wavelength
        self.allmodes_b = None
        self.allmodes_l = None
        self.allmodes_m = None
        self.nmodes = nmodes
        self.max_r = None
        self.datadir = datadir

        self.all_smpowers = []
        self.all_mmpowers = []
        self.all_mmphases = []
        self.all_smphases = []
        self.Cmat = None # Transfer matrix from SM to MM
        self.Dmat = None # Transfer matrix from MM to SM

        if n_core is not None:
            self.NA = ofiber.numerical_aperture(n_core, n_cladding)
            self.V = ofiber.V_parameter(core_radius, self.NA, wavelength)


    def find_fiber_modes(self, max_l=100):
        """
        Finds LP modes for the specified fiber.

        Parameters
        ----------
        max_l
            Maximum number of l modes to find (can be arbitrarily large)
        """
        allmodes_b = []
        allmodes_l = []
        allmodes_m = []
        for l in range(max_l):
            cur_b = ofiber.LP_mode_values(self.V, l)
            if len(cur_b) == 0:
                break
            else:
                allmodes_b.extend(cur_b)
                ls = (np.ones_like(cur_b)) * l
                allmodes_l.extend(ls.astype(int))
                ms = np.arange(len(cur_b))+1
                allmodes_m.extend(ms)

        allmodes_b = np.asarray(allmodes_b)
        nLPmodes = len(allmodes_b)
        # print('Total number of LP modes found: %d' % nLPmodes)
        l = np.asarray(allmodes_l)
        total_unique_modes = len(np.where(l == 0)[0]) + len(np.where(l > 0)[0])*2
        print('Total number of unique modes found: %d' % total_unique_modes)
        self.allmodes_b = allmodes_b
        self.allmodes_l = allmodes_l
        self.allmodes_m = allmodes_m
        self.nLPmodes = nLPmodes


    def make_fiber_modes(self, max_r=2, npix=100, zlim=0.04, show_plots=False):
        """
        Calculate the LP mode fields, and store as polar and cartesian amplitude maps

        Parameters
        ----------
        max_r
            Maximum radius to calculate mode field, where r=1 is the core diameter
        npix
            Width of mode field calculation in pixels
        zlim
            Maximum value to plot
        show_plots : bool
            Whether to produce a plot for each mode
        """

        r = np.linspace(0, max_r, npix) # Radial positions, normalised so core_radius = 1
        self.max_r = max_r
        self.allmodefields_cos_polar = []
        self.allmodefields_cos_cart = []
        self.allmodefields_sin_polar = []
        self.allmodefields_sin_cart = []

        for mode_to_calc in range(self.nLPmodes):
            field_1d = ofiber.LP_radial_field(self.V, self.allmodes_b[mode_to_calc],
                                              self.allmodes_l[mode_to_calc], r)

            phivals = np.linspace(0, 2*np.pi, npix)
            phi_cos = np.cos(self.allmodes_l[mode_to_calc] * phivals)
            phi_sin = np.sin(self.allmodes_l[mode_to_calc] * phivals)

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

            self.allmodefields_cos_polar.append(field_cos)
            self.allmodefields_cos_cart.append(field_cos_cart)
            self.allmodefields_sin_polar.append(field_sin)
            self.allmodefields_sin_cart.append(field_sin_cart)

            if show_plots:
                self.plot_fiber_modes(mode_to_calc, zlim)
                plt.pause(0.1)


    def plot_fiber_modes(self, mode_to_plot, zlim=0.04, fignum=1):
        """
        Make a plot of the cos and sin amplitudes of a given mode

        Parameters
        ----------
        mode_to_plot
            Number of mode to plot
        zlim
            Maximum value to plot
        """
        plt.figure(1)
        plt.clf()
        plt.subplot(121)
        sz = self.max_r * self.core_radius
        plt.imshow(self.allmodefields_cos_cart[mode_to_plot], extent=(-sz, sz, -sz, sz), cmap='bwr',
                   vmin=-zlim, vmax=zlim)
        plt.xlabel('Position ($\mu$m)')
        plt.ylabel('Position ($\mu$m)')
        plt.title('Mode l=%d, m=%d (cos)' % (self.allmodes_l[mode_to_plot], self.allmodes_m[mode_to_plot]))
        core_circle = plt.Circle((0,0), self.core_radius, color='k', fill=False, linestyle='--', alpha=0.2)
        plt.gca().add_patch(core_circle)
        plt.subplot(122)
        sz = self.max_r * self.core_radius
        plt.imshow(self.allmodefields_sin_cart[mode_to_plot], extent=(-sz, sz, -sz, sz), cmap='bwr',
                   vmin=-zlim, vmax=zlim)
        plt.xlabel('Position ($\mu$m)')
        plt.title('Mode l=%d, m=%d (sin)' % (self.allmodes_l[mode_to_plot], self.allmodes_m[mode_to_plot]))
        core_circle = plt.Circle((0,0), self.core_radius, color='k', fill=False, linestyle='--', alpha=0.2)
        plt.gca().add_patch(core_circle)
        plt.pause(0.001)
        print('LP mode %d, %d' % (self.allmodes_l[mode_to_plot], self.allmodes_m[mode_to_plot]))


    def load_rsoft_data_sm2mm(self, rsoft_datadir, rsoft_fileprefix, show_plots=False,
                        av_fluxes=100, offset_sm_meas=100, save_output=False, zero_phases=True,
                        fignum=2):
        """
        Load rsoft outputs for a single-mode to multi-mode propagation. Assumes one SM waveguide
        excited at a time.

        Parameters
        ----------
        rsoft_datadir
        rsoft_fileprefix
        show_plots : bool
        av_fluxes
            If >0, average n flux measurements from monitor
        offset_sm_meas
            Measure SM fluxes starting at this index. Useful to skip first part of monitor
            since still coupling.
        save_output : bool
            Save the measured powers and phases to a npz file
        zero_phases : bool
            Set the SM input phases to zero
        fignum
            Figure number in which to display plots
        """

        # Specify relevant indices of MONdata:
        sm_power_monrange = (0, 19)
        mm_power_monrange = (19, 38)
        mm_phase_monrange = (38, 57)

        nwgs = self.nmodes
        for wgnum in range(nwgs):
            rsoft_filename = rsoft_fileprefix + '%.2d' % (wgnum+1)
            print('Reading rsoft files ' + rsoft_filename)
            r = Rsoftdata(rsoft_datadir)
            r.readall(filename=rsoft_filename)
            if show_plots:
                r.plotall()
                plt.pause(0.001)

            smpower_mons = r.MONdata[:, sm_power_monrange[0]:sm_power_monrange[1]]
            mmpower_mons = r.MONdata[:, mm_power_monrange[0]:mm_power_monrange[1]]
            mmphase_mons = r.MONdata[:, mm_phase_monrange[0]:mm_phase_monrange[1]]
            if av_fluxes > 0:
                f = smpower_mons[offset_sm_meas:av_fluxes+offset_sm_meas, :]
                smpower = f.mean(axis=0)
                f = mmpower_mons[-av_fluxes:, :]
                mmpower = f.mean(axis=0)
            else:
                smpower = smpower_mons[offset_sm_meas, :]
                mmpower = mmpower_mons[-1, :]
            mmphase = mmphase_mons[-1, :]

            if zero_phases:
                smphases = np.zeros(nwgs)

            self.all_smpowers.append(smpower)
            self.all_smphases.append(smphases)
            self.all_mmpowers.append(mmpower)
            self.all_mmphases.append(mmphase)

        if save_output:
            outfilename = self.datadir + 'extractedvals_' + rsoft_fileprefix + '.npz'
            np.savez(outfilename, all_smpowers=self.all_smpowers, all_mmpowers=self.all_mmpowers,
                     all_mmphases=self.all_mmphases, all_smphases=self.all_smphases)

        if show_plots:
            plt.figure(fignum)
            plt.clf()
            plt.subplot(211)
            plt.imshow(np.asarray(self.all_mmpowers))
            plt.colorbar()
            plt.title('Output mode power')
            plt.ylabel('Excited waveguide no.')
            plt.xlabel('Mode no.')
            plt.subplot(212)
            plt.imshow(np.asarray(self.all_mmphases), cmap='twilight_shifted')
            plt.colorbar()
            plt.title('Output mode phase')
            plt.ylabel('Excited waveguide no.')
            plt.xlabel('Mode no.')
            plt.tight_layout()


    def load_rsoft_data_mm2sm(self, rsoft_datadir, rsoft_fileprefix, show_plots=False,
                              av_fluxes=100, save_output=False, fignum=2):
        """
        Load rsoft outputs for a multi-mode to single-mode propagation. Assumes one MM mode
        excited at a time.

        Parameters
        ----------
        rsoft_datadir
        rsoft_fileprefix
        show_plots : bool
        av_fluxes
            If >0, average n flux measurements from monitor
        save_output : bool
            Save the measured powers and phases to a npz file
        fignum
            Figure number in which to display plots
        """

        # Specify mode indices
        LP_modes = np.array([[0,1],
                             [0,2],
                             [0,3],
                             [1,1],
                             [-1,1],
                             [1,2],
                             [-1,2],
                             [2,1],
                             [-2,1],
                             [2,2],
                             [-2,2],
                             [3,1],
                             [-3,1],
                             [3,2],
                             [-3,2],
                             [4,1],
                             [-4,1],
                             [5,1],
                             [-5,1]
                             ])


        # Specify relevant indices of MONdata
        sm_power_monrange = (0, 19)
        mm_power_monrange = (19, 38)
        mm_phase_monrange = (38, 57)
        sm_phase_monrange = (57, 76)

        nwgs = self.nmodes
        for wgnum in range(nwgs):
            rsoft_suffix = 'LP%d%d' % (LP_modes[wgnum,0], LP_modes[wgnum,1])
            rsoft_filename = rsoft_fileprefix + rsoft_suffix
            print('Reading rsoft files ' + rsoft_filename)
            r = Rsoftdata(rsoft_datadir)
            r.readall(filename=rsoft_filename)
            if show_plots:
                r.plotall()
                plt.pause(0.001)

            smpower_mons = r.MONdata[:, sm_power_monrange[0]:sm_power_monrange[1]]
            mmpower_mons = r.MONdata[:, mm_power_monrange[0]:mm_power_monrange[1]]
            mmphase_mons = r.MONdata[:, mm_phase_monrange[0]:mm_phase_monrange[1]]
            smphase_mons = r.MONdata[:, sm_phase_monrange[0]:sm_phase_monrange[1]]

            if av_fluxes > 0:
                f = smpower_mons[-av_fluxes:, :]
                smpower = f.mean(axis=0)
                f = mmpower_mons[0:av_fluxes, :]
                mmpower = f.mean(axis=0)
            else:
                smpower = smpower_mons[-1, :]
                mmpower = mmpower_mons[0, :]
            mmphase = mmphase_mons[0, :]
            smphase = smphase_mons[-1, :]

            self.all_smpowers.append(smpower)
            self.all_mmpowers.append(mmpower)
            self.all_mmphases.append(mmphase)
            self.all_smphases.append(smphase)

        if save_output:
            outfilename = self.datadir + 'extractedvals_' + rsoft_fileprefix + '.npz'
            np.savez(outfilename, all_smpowers=self.all_smpowers, all_mmpowers=self.all_mmpowers,
                     all_mmphases=self.all_mmphases, all_smphases=self.all_smphases)

        if show_plots:
            plt.figure(fignum)
            plt.clf()
            plt.subplot(211)
            plt.imshow(np.asarray(self.all_smpowers))
            plt.colorbar()
            plt.title('Output waveguide power')
            plt.ylabel('Excited mode no.')
            plt.xlabel('Waveguide no.')
            plt.subplot(212)
            plt.imshow(np.asarray(self.all_smphases), cmap='twilight_shifted')
            plt.colorbar()
            plt.title('Output waveguide phase')
            plt.ylabel('Excited mode no.')
            plt.xlabel('Waveguide no.')
            plt.tight_layout()



    def load_savedvalues(self, filename):
        """
        Load previously saved powers and phases from npz file

        Parameters
        ----------
        filename
        """
        f = np.load(self.datadir+filename)
        self.all_smpowers = f['all_smpowers']
        self.all_mmpowers = f['all_mmpowers']
        self.all_mmphases = f['all_mmphases']
        try:
            self.all_smphases = f['all_smphases']
        except:
            self.all_smphases = None


    def set_mmvals_nominal(self, square=False):
        """
        Set all appropriate mm powers and phases to 1 and 0 respectively.
        """
        if square:
            self.all_mmpowers = np.diag(np.ones(self.nmodes))
            self.all_mmphases = np.zeros((self.nmodes, self.nmodes))
        else:
            self.all_mmpowers = np.zeros(self.nmodes)
            self.all_mmpowers[0] = 1
            self.all_mmphases = np.zeros(self.nmodes)


    def set_smvals_nominal(self, square=False):
        """
        Set all appropriate sm powers and phases to 1 and 0 respectively.
        """
        if square:
            self.all_smpowers = np.diag(np.ones(self.nmodes))
            self.all_smphases = np.zeros((self.nmodes, self.nmodes))
        else:
            self.all_smpowers = np.zeros(self.nmodes)
            self.all_smpowers[0] = 1
            self.all_smphases = np.zeros(self.nmodes)


    def make_transfer_matrix_sm2mm(self, sm_phase=0):
        """
        Calculate the transfer matrix Cmat from SM-to-MM simulations.
        This requires that for each measurement, only one waveguide is excited

        Parameters
        ----------
        sm_phase
            Phase to set the SM waveguide excitation to
        """
        nwgs = len(self.all_smpowers)
        Cmat = np.zeros((nwgs,nwgs), dtype=np.cfloat)
        for col_num in range(nwgs):
            sm_val = self.all_smpowers[col_num][col_num] * np.exp(1j*sm_phase)
            col = self.all_mmpowers[col_num] * np.exp(1j*self.all_mmphases[col_num]/180*np.pi) / sm_val
            Cmat[:,col_num] = col
        self.Cmat = Cmat


    def make_transfer_matrix_mm2sm(self, mm_phase=None):
        """
        Calculate the transfer matrix Dmat from MM-to-SM simulations.
        This requires that for each measurement, only one mode is excited

        Parameters
        ----------
        mm_phase
            Phase to set MM mode phase to. If None, use phase from data.
        """
        nwgs = len(self.all_mmpowers)
        Dmat = np.zeros((nwgs,nwgs), dtype=np.cfloat)
        for col_num in range(nwgs):
            if mm_phase is not None:
                mm_val = self.all_mmpowers[col_num][col_num] * np.exp(1j*mm_phase)
            else:
                mm_val = self.all_mmpowers[col_num][col_num] * np.exp(1j*self.all_mmphases[col_num][col_num])
            col = self.all_smpowers[col_num] * np.exp(1j*self.all_smphases[col_num]/180*np.pi) / mm_val
            Dmat[:,col_num] = col
        self.Dmat = Dmat


    def load_rsoft_data_sm2mm_single(self, rsoft_datadir, rsoft_fileprefix, show_plots=False,
                              av_fluxes=100, offset_sm_meas=100, zero_phases=True):
        """
        Load rsoft outputs for a single propagation.

        Parameters
        ----------
        rsoft_datadir
        rsoft_fileprefix
        show_plots : bool
        av_fluxes
            If >0, average n flux measurements from monitor
        offset_sm_meas
            Measure SM fluxes starting at this index. Useful to skip first part of monitor
            since still coupling.
        """

        self.all_smpowers = []
        self.all_mmpowers = []
        self.all_mmphases = []
        self.all_smphases = []

        # Specify relevant indices of MONdata:
        sm_power_monrange = (0, 19)
        mm_power_monrange = (19, 38)
        mm_phase_monrange = (38, 57)

        rsoft_filename = rsoft_fileprefix
        print('Reading rsoft files ' + rsoft_filename)
        r = Rsoftdata(rsoft_datadir)
        r.readall(filename=rsoft_filename)
        if show_plots:
            r.plotall()
            plt.pause(0.001)

        self.MONdata = r.MONdata ## TESTING

        smpower_mons = r.MONdata[:, sm_power_monrange[0]:sm_power_monrange[1]]
        mmpower_mons = r.MONdata[:, mm_power_monrange[0]:mm_power_monrange[1]]
        mmphase_mons = r.MONdata[:, mm_phase_monrange[0]:mm_phase_monrange[1]]
        if av_fluxes > 0:
            f = smpower_mons[offset_sm_meas:av_fluxes+offset_sm_meas, :]
            smpower = f.mean(axis=0)
            f = mmpower_mons[-av_fluxes:, :]
            mmpower = f.mean(axis=0)
        else:
            smpower = smpower_mons[offset_sm_meas, :]
            mmpower = mmpower_mons[-1, :]
        mmphase = mmphase_mons[-1, :]

        if zero_phases:
            smphases = np.zeros(self.nmodes)

        self.all_smpowers.append(smpower)
        self.all_smphases.append(smphases)
        self.all_mmpowers.append(mmpower)
        self.all_mmphases.append(mmphase)


    def test_matrix(self, matrix, input_powers, input_phases, output_powers, output_phases, fignum=1, pausetime=1):
        """
        Plot the outputs from a set of inputs using the transfer matrix, and overplot
        the 'true' values measured.

        Parameters
        ----------
        matrix
            Matrix to use for testing
        input_powers
            List or array of input powers (one row / item per measurement)
        input_phases
            List or array of input phases (one row / item per measurement)
        output_powers
            List or array of output powers (one row / item per measurement)
        output_phases
            List or array of output phases (one row / item per measurement)
        fignum
            Figure number in which to display plots
        pausetime
            Time (s) to wait between plots
        """
        plt.figure(fignum)
        nmeas = len(input_powers)
        for n in range(nmeas):
            plt.clf()
            plt.subplot(211)
            plt.plot(output_powers[n],'-x')
            input_val = input_powers[n] * np.exp(1j*input_phases[n]/np.pi*180)
            # plt.plot(np.abs(np.matmul(matrix, input_powers[n])),'-+')
            plt.plot(np.abs(np.matmul(matrix, input_val)),'-+')
            plt.subplot(212)
            plt.plot(output_phases[n],'-x')
            # newangs = np.angle(np.matmul(matrix, input_powers[n]))/np.pi*180
            newangs = np.angle(np.matmul(matrix, input_val))/np.pi*180
            newangs[newangs < 0] = newangs[newangs < 0]+360
            plt.plot(newangs,'-+')
            plt.pause(pausetime)


    def read_ind_file(self, rsoft_datadir, ind_filename, skipfirst=True):
        """
        Read useful info from an rsoft .ind file, such as the power and phase of launch
        fields.

        Parameters
        ----------
        rsoft_datadir
        rsoft_fileprefix

        Returns
        -------
        all_powers
        all_phases
        """
        with open(rsoft_datadir+ind_filename, 'r') as indfile:
            indfile_data = indfile.readlines()

        all_powers = []
        all_phases = []
        if skipfirst:
            skip_this_power = True # Skip the first occurrences, since they are repeated.
            skip_this_phase = True
        else:
            skip_this_power = False
            skip_this_phase = False
        for line in indfile_data:
            if 'launch_power' in line:
                if skip_this_power:
                    skip_this_power = False
                else:
                    ll = line.split('=')
                    power = np.float(ll[1].strip())
                    all_powers.append(power)
            if 'launch_phase' in line:
                if skip_this_phase:
                    skip_this_phase = False
                else:
                    ll = line.split('=')
                    phase = np.float(ll[1].strip())
                    all_phases.append(phase)

        return all_powers, all_phases












