from ctypes import byref, c_double, c_int, c_uint64, CDLL, POINTER, Structure
from glob import glob
from logging import getLogger
from pathlib import Path

from numpy import asarray, ones, zeros
from numpy.ctypeslib import ndpointer

from pyrad.lbl.hitran import Hitran, Voigt
from pyrad.lbl.tips import TotalPartitionFunction


info = getLogger(__name__).info
library = glob(str(Path(__file__).parent / "libgrtoptics*.so"))[0]
gas_optics = CDLL(library)
HOST = -1


class Gas(object):
    def __init__(self, formula, hitran_database=None, isotopologues=None,
                 line_profile=Voigt(), tips_database=None, device="host"):
        self.device = -1 if device.lower() == "host" else device
        self.formula = formula
        database = Hitran(formula, line_profile, isotopologues, hitran_database)
        partition_function = TotalPartitionFunction(formula, tips_database)
        self.spectral_lines = database.spectral_lines(partition_function)
        self.num_iso = len(database.isotopologues)
        self.mol_id = database.molecule_id

    def absorption_coefficient(self, pressure, temperature, volume_mixing_ratio, grid):
        """Calculates absorption coefficients for the gas using GRTCODE.
        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            spectral_grid: Wavenumber grid [cm-1].
            line_cut_off: Cut-off from spectral line center [cm-1].

        Returns:
            Absorption coefficients [m2].
        """
        info("Calculating GRTCODE line-by-line spectra for {}.".format(self.formula))
        #Convert pressure from Pa to atm.
        p = asarray([pressure,])*9.86923e-6
        t = asarray([temperature,])
        vmr = asarray([volume_mixing_ratio,])

        #Calcuate line center.
        center = pressure_shift_line_centers(self.spectral_lines.v,
                                             self.spectral_lines.d_air, p)

        #Calculate total partition functions.
        partition_function = total_partition_functions(self.mol_id, self.num_iso, t)

        #Calculate temperature-corrected line strengths.
        strength = correct_line_strengths(self.num_iso, self.spectral_lines.iso,
                                          self.spectral_lines.s, self.spectral_lines.v,
                                          self.spectral_lines.en, t, partition_function)

        #Calcluate lorentz half-widths.
        gamma = lorentz_halfwidths(self.spectral_lines.n_air, self.spectral_lines.gamma_air,
                                   self.spectral_lines.gamma_self, t, p, vmr)

        #Calculate doppler half-widths.
        alpha = doppler_halfwidths(self.spectral_lines.mass[0], self.spectral_lines.v, t)

        #Calculate the molecule's absorption coefficients.
        bins = create_spectral_bins(p.size, grid[0], grid.size, grid[1] - grid[0], 1.5)
        k = absorption_coefficients(center, strength, gamma, alpha, bins)
        destroy_spectral_bins(bins)
        return k


class SpectralBins(Structure):
    _fields_ = [("num_layers", c_int),
                ("w0", c_double),
                ("wres", c_double),
                ("num_wpoints", c_uint64),
                ("n", c_uint64),
                ("width", c_double),
                ("isize", c_uint64),
                ("ppb", c_int),
                ("do_interp", c_int),
                ("last_ppb", c_int),
                ("do_last_interp", c_int),
                ("w", POINTER(c_double)),
                ("tau", POINTER(c_double)),
                ("l", POINTER(c_uint64)),
                ("r", POINTER(c_uint64)),
                ("device", c_int)]


def create_spectral_bins(num_layers, w0, n, wres, bin_width, device=HOST):
    gas_optics.create_spectral_bins.argtypes = [POINTER(SpectralBins), c_int, c_double,
                                                c_uint64, c_double, c_double, c_int]
    bins = SpectralBins()
    rc = gas_optics.create_spectral_bins(byref(bins), c_int(num_layers), c_double(w0),
                                         c_uint64(n), c_double(wres), c_double(bin_width),
                                         c_int(device))
    return bins


def destroy_spectral_bins(bins):
    gas_optics.destroy_spectral_bins.argtypes = [POINTER(SpectralBins),]
    gas_optics.destroy_spectral_bins(byref(bins))


def pressure_shift_line_centers(center, shift, pressure):
    num_layers, num_lines = pressure.size, center.size
    shifted_center = zeros((num_layers, num_lines))
    gas_optics.calc_line_centers.argtypes = [c_uint64, c_int] + \
                                            4*[ndpointer(c_double, flags="C_CONTIGUOUS")]
    gas_optics.calc_line_centers(c_uint64(num_lines), c_int(num_layers),
                                 center, shift, pressure, shifted_center)
    return shifted_center


def total_partition_functions(mol_id, num_iso, temperature):
    num_layers = temperature.size
    partition_function = zeros((num_layers, num_iso))
    gas_optics.calc_partition_functions.argtypes = 3*[c_int] + \
                                                   2*[ndpointer(c_double, flags="C_CONTIGUOUS")]
    gas_optics.calc_partition_functions(c_int(num_layers), c_int(mol_id), c_int(num_iso),
                                        temperature, partition_function)
    return partition_function


def correct_line_strengths(num_iso, iso, strength_t0, center, energy, temperature,
                           partition_function):
    num_layers, num_lines = temperature.size, center.size
    strength = zeros((num_layers, num_lines))
    iso_id = asarray(iso, dtype=c_int)
    gas_optics.calc_line_strengths.argtypes = [c_uint64, c_int, c_int,
                                               ndpointer(c_int, flags="C_CONTIGUOUS")] + \
                                              6*[ndpointer(c_double, flags="C_CONTIGUOUS")]
    gas_optics.calc_line_strengths(c_uint64(num_lines), c_int(num_layers), c_int(num_iso),
                                   iso_id, strength_t0, center, energy, temperature,
                                   partition_function, strength)
    return strength


def lorentz_halfwidths(n, yair, yself, temperature, pressure, volume_mixing_ratio):
    num_layers, num_lines = temperature.size, yair.size
    gamma = zeros((num_layers, num_lines))
    ps = volume_mixing_ratio[:]*pressure[:]
    gas_optics.calc_lorentz_hw.argtypes = [c_uint64, c_int] + \
                                          7*[ndpointer(c_double, flags="C_CONTIGUOUS")]
    gas_optics.calc_lorentz_hw(c_uint64(num_lines), c_int(num_layers), n,
                               yair, yself, temperature, pressure, ps, gamma)
    return gamma


def doppler_halfwidths(mass, center, temperature):
    num_layers, num_lines = temperature.size, center.size
    alpha = zeros((num_layers, num_lines))
    gas_optics.calc_doppler_hw.argtypes = [c_uint64, c_int, c_double] + \
                                          3*[ndpointer(c_double, flags="C_CONTIGUOUS")]
    gas_optics.calc_doppler_hw(c_uint64(num_lines), c_int(num_layers), c_double(mass),
                               center, temperature, alpha)
    return alpha


def sort_lines(center, strength, gamma, alpha):
    num_layers, num_lines = center.shape
    gas_optics.sort_lines.argtypes = [c_uint64, c_int] + \
                                     4*[ndpointer(c_double, flags="C_CONTIGUOUS")]
    gas_optics.sort_lines(c_uint64(num_lines), c_int(num_layers), center, strength,
                          gamma, alpha)


def absorption_coefficients(center, strength, gamma, alpha, bins):
    num_layers, num_lines = center.shape
    absorption_coefficient = zeros((num_layers, bins.num_wpoints))
    n = ones(num_layers)
    sort_lines(center, strength, gamma, alpha)
    gas_optics.calc_optical_depth_line_sample.argtypes = [c_uint64, c_int] + \
                                                        5*[ndpointer(c_double, flags="C_CONTIGUOUS")] + \
                                                        [SpectralBins,
                                                         ndpointer(c_double, flags="C_CONTIGUOUS")]
    gas_optics.calc_optical_depth_line_sample(c_uint64(num_lines), c_int(num_layers),
                                              center, strength, gamma, alpha, n,
                                              bins, absorption_coefficient)
    return absorption_coefficient