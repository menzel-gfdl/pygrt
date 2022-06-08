from ctypes import byref, c_double, c_int, c_uint64, CDLL, POINTER, Structure
from glob import glob
from logging import getLogger
from pathlib import Path

from numpy import asarray, ones, zeros
from numpy.ctypeslib import ndpointer


info = getLogger(__name__).info
library = glob(str(Path(__file__).parent / "libgrtoptics*.so"))[0]
gas_optics = CDLL(library)
HOST = -1
GRTCODE_SUCCESS = 0
cm_to_m = 0.01  # [m cm-1].
avogadro = 6.0221409e23 # Avogadro's constant [mol-1].


def cpointer(var_type):
    """Creates a contiguous numpy.ctypeslib ndpointer.

    Args:
        var_type: ctypes variable type object.

    Returns:
        A contiguous ndpointer object.
    """
    return ndpointer(var_type, flags="C_CONTIGUOUS")


class SpectralLines(object):
    """HITRAN spectral line parameters.

    Attributes:
        d_air: Transition pressure shift factor.
        en: Lower state energy level.
        gamma_air: Air-broadened halfwidth.
        gamma_self: Self-broadened halfwidth.
        iso: Isotopologue id.
        n_air: Air-broadening temperature dependence power.
        s: Transition strength.
        v: Transition wavenumber.
    """
    parameters = ["d_air", "en", "gamma_air", "gamma_self", "iso", "n_air", "s", "v"]

    def __init__(self, transitions):
        """Converts from Transition objects to numpy arrays.

        Args:
            transitions: List of Transition objects.
        """
        for x in self.parameters:
            setattr(self, x, [])
        for transition in transitions:
            self.d_air.append(transition.delta_air)
            self.en.append(transition.elower)
            self.gamma_air.append(transition.gamma_air)
            self.gamma_self.append(transition.gamma_self)
            self.iso.append(transition.local_iso_id)
            self.n_air.append(transition.n_air)
            self.s.append(transition.sw)
            self.v.append(transition.nu)
        for x in self.parameters:
            setattr(self, x, asarray(getattr(self, x)))


class Gas(object):
    def __init__(self, formula, mass, transitions, total_parition_function=None, device="host"):
        """Initializes object.

        Args:
            formula: String chemical formula.
            mass: List of isotopologue masses.
            transitions: List of TransitionTable objects.
            total_partition_function: Not used.
            device: Device (host or GPU id) to run GRTcode on.
        """
        self.device = HOST if device.lower() == "host" else device
        self.spectal_lines = SpectralLines(transitions)
        self.mass = asarray(mass)
        self.mol_id = transitions[0].molecule_id
        self.num_iso = len(mass)

        # Do initial correction to the line strengths.
        q = total_partition_functions(self.mol_id, self.num_iso, asarray([296.,]))
        self.spectral_lines.s[:] /= correct_line_strengths(self.num_iso,
                                                           self.spectral_lines.iso,
                                                           ones(self.spectral_lines.v.size),
                                                           self.spectral_lines.v,
                                                           self.spectral_lines.en,
                                                           asarray([296.,]), q)[0, :]

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid,
                               remove_pedestal=False):
        """Calculates absorption coefficients for the gas using GRTCODE.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Wavenumber grid [cm-1].

        Returns:
            Absorption coefficients [m2].
        """
        info("Calculating GRTCODE optics for molecule {}.".format(self.mol_id))

        #Convert pressure from Pa to atm.
        p = asarray([pressure*9.86923e-6,])
        t = asarray([temperature,])
        vmr = asarray([volume_mixing_ratio,])
        pedestal = (0., 2.e4) if remove_pedestal and self.mol_id == 1 else None

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
        alpha = doppler_halfwidths(self.mass[self.iso[0] - 1], self.spectral_lines.v, t)

        #Calculate the molecule's absorption coefficients.
        bins = create_spectral_bins(p.size, grid[0], grid.size, grid[1] - grid[0], 1.5)
        k = absorption_coefficients(center, strength, gamma, alpha, bins, pedestal)
        destroy_spectral_bins(bins)
        return k*cm_to_m*cm_to_m


class SpectralBins(Structure):
    """Binding to a GRTcode SpectralBins c struct."""
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
    """Initializes a GRTcode SpectralBins c struct.

    Args:
        num_layers: Number of layers.
        w0: Spectral grid lower bound [cm-1].
        n: Spectral grid size.
        wres: Spectral grid resolution [cm-1].
        bin_width: Width of each spectral bin [cm-1].
        device: Device (host or GPU id) that GRTcode will run on.

    Returns:
        A SpectralBins object (binded to an initialized GRTcode SpectralBins c struct.
    """
    types = [POINTER(SpectralBins), c_int, c_double, c_uint64, c_double, c_double, c_int]
    gas_optics.create_spectral_bins.argtypes = types
    gas_optics.create_spectral_bins.restype = check_return_code
    bins = SpectralBins()
    gas_optics.create_spectral_bins(byref(bins), c_int(num_layers), c_double(w0),
                                    c_uint64(n), c_double(wres), c_double(bin_width),
                                    c_int(device))
    return bins


def destroy_spectral_bins(bins):
    """Frees memory held by a GRTcode SpectralBins c struct.

    Args:
        bins: A SpectralBins object.
    """
    gas_optics.destroy_spectral_bins.argtypes = [POINTER(SpectralBins),]
    gas_optics.destroy_spectral_bins.restype = check_return_code
    gas_optics.destroy_spectral_bins(byref(bins))


def pressure_shift_line_centers(center, shift, pressure):
    """Pressure shifts the transition wavenumbers.

    Args:
        center: Numpy array of transition wavenumbers [cm-1].
        shift: Numpy array of pressure-shift factors [cm-1 atm-1].
        pressure: Numpy array of pressures [atm].

    Returns:
        Numpy array of pressure-shifted transition wavenumbers [cm-1].
    """
    num_layers, num_lines = pressure.size, center.size
    shifted_center = zeros((num_layers, num_lines))
    gas_optics.calc_line_centers.argtypes = [c_uint64, c_int] + 4*[cpointer(c_double)]
    gas_optics.calc_line_centers.restype = check_return_code
    gas_optics.calc_line_centers(c_uint64(num_lines), c_int(num_layers),
                                 center, shift, pressure, shifted_center)
    return shifted_center


def total_partition_functions(mol_id, num_iso, temperature):
    """Calculates the total partition functions.

    Args:
        mol_id: HITRAN molecule id.
        num_iso: Number of isotopologues for the molecule.
        temperature: Numpy array of temperatures [K].

    Returns:
        Numpy array of total partition functions.
    """
    num_layers = temperature.size
    partition_function = zeros((num_layers, num_iso))
    gas_optics.calc_partition_functions.argtypes = 3*[c_int] + 2*[cpointer(c_double)]
    gas_optics.calc_partition_functions.restype = check_return_code
    gas_optics.calc_partition_functions(c_int(num_layers), c_int(mol_id), c_int(num_iso),
                                        temperature, partition_function)
    return partition_function


def correct_line_strengths(num_iso, iso, strength_t0, center, energy, temperature,
                           partition_function):
    """Temperature corrects the transition strengths.

    Args:
        num_iso: Number of isotopologues for the molecule.
        iso: Numpy array of HITRAN isotopologue ids.
        strength_t0: Numpy array of transition strengths [cm].
        center: Numpy array of transition wavenumbers [cm-1].
        energy: Numpy array of lower state energies [cm-1].
        temperature: Numpy array of temperatures [K].
        partition_function: Numpy array of total partition functions.

    Returns:
        Numpy array of temperature-corrected transition strengths [cm].
    """
    num_layers, num_lines = temperature.size, center.size
    strength = zeros((num_layers, num_lines))
    iso_id = asarray(iso, dtype=c_int)
    types = [c_uint64, c_int, c_int, cpointer(c_int)] + 6*[cpointer(c_double)]
    gas_optics.calc_line_strengths.argtypes = types
    gas_optics.calc_line_strengths.restype = check_return_code
    gas_optics.calc_line_strengths(c_uint64(num_lines), c_int(num_layers), c_int(num_iso),
                                   iso_id, strength_t0, center, energy, temperature,
                                   partition_function, strength)
    return strength


def lorentz_halfwidths(n, yair, yself, temperature, pressure, volume_mixing_ratio):
    """Calculates Lorentz halfwidths.

    Args:
        n: Numpy array of temperature dependence powers.
        yair: Numpy array of air-broadened halfwidths [cm-1 atm-1] at 296 K and 1 atm.
        yself: Numpy array of self-broadened halfwidths [cm-1 atm-1] at 296 K and 1 atm.
        temperature: Numpy array of temperatures [K].
        pressure: Numpy array of pressures [atm].
        volume_mixing_ratio: Numpy array of volume mixing ratios [mol mol-1].

    Returns:
        Numpy array of Lorentz halfwidths [cm-1].
    """
    num_layers, num_lines = temperature.size, yair.size
    gamma = zeros((num_layers, num_lines))
    ps = volume_mixing_ratio[:]*pressure[:]
    gas_optics.calc_lorentz_hw.argtypes = [c_uint64, c_int] + 7*[cpointer(c_double)]
    gas_optics.calc_lorentz_hw.restype = check_return_code
    gas_optics.calc_lorentz_hw(c_uint64(num_lines), c_int(num_layers), n,
                               yair, yself, temperature, pressure, ps, gamma)
    return gamma


def doppler_halfwidths(mass, center, temperature):
    """Calculates Doppler halfwidths.

    Args:
        mass: Molar mass of the molecule [g mol-1].
        center: Numpy array of transition wavenumbers [cm-1].
        temperature: Numpy array of temperatures [K].

    Returns:
        Numpy array of Doppler halfwidths [cm-1].
    """
    num_layers, num_lines = temperature.size, center.size
    alpha = zeros((num_layers, num_lines))
    gas_optics.calc_doppler_hw.argtypes = [c_uint64, c_int, c_double] + 3*[cpointer(c_double)]
    gas_optics.calc_doppler_hw.restype = check_return_code
    gas_optics.calc_doppler_hw(c_uint64(num_lines), c_int(num_layers),
                               c_double(mass/avogadro), center, temperature, alpha)
    return alpha


def sort_lines(center, strength, gamma, alpha):
    """Sorts transitions in order of increasing wavenumber.

    Args:
        center: Numpy array of transition wavenumbers [cm-1].
        strength: Numpy array of transitin strengths [cm].
        gamma: Numpy array of Lorentz halfwidths [cm-1].
        alpha: Numpy array of Doppler halfwidths [cm-1].
    """
    num_layers, num_lines = center.shape
    gas_optics.sort_lines.argtypes = [c_uint64, c_int] + 4*[cpointer(c_double)]
    gas_optics.sort_lines.restype = check_return_code
    gas_optics.sort_lines(c_uint64(num_lines), c_int(num_layers), center, strength,
                          gamma, alpha)


def absorption_coefficients(center, strength, gamma, alpha, bins, pedestal_bounds=None):
    """Calculates absorption coefficients.

    Args:
        center: Numpy array of transition wavenumbers [cm-1].
        strength: Numpy array of transitin strengths [cm].
        gamma: Numpy array of Lorentz halfwidths [cm-1].
        alpha: Numpy array of Doppler halfwidths [cm-1].
        bins: SpectralBins object.
        pedestal_bounds: Tuple describing wavenumber [cm-1] range where pedestal is removed.

    Returns:
        Numpy array of absorption coefficients [cm2].
    """
    num_layers, num_lines = center.shape
    absorption_coefficient = zeros((num_layers, bins.num_wpoints))
    n = ones(num_layers)
    lower = None if pedestal_bounds is None else byref(c_double(pedestal_bounds[0]))
    upper = None if pedestal_bounds is None else byref(c_double(pedestal_bounds[1]))
    types = [c_uint64, c_int] + 5*[cpointer(c_double)] + [SpectralBins,] + \
            [cpointer(c_double),] + 2*[POINTER(c_double)]
    gas_optics.calc_optical_depth_line_sample.argtypes = types
    gas_optics.calc_optical_depth_line_sample.restype = check_return_code
    gas_optics.calc_optical_depth_line_sample(c_uint64(num_lines), c_int(num_layers),
                                              center, strength, gamma, alpha, n,
                                              bins, absorption_coefficient, lower, upper)
    return absorption_coefficient


def water_vapor_continuum(temperature, pressure, volume_mixing_ratio, grid, c_foreign,
                          c_self, t_foreign, t_self):
    """Calculates water vapor continnum absorption coefficients.

    Args:
        temperature: Numpy array of temperatures [K].
        pressure: Numpy array of pressures [atm].
        volume_mixing_ratio: Numpy array of volume mixing ratios [mol mol-1].
        grid: Numpy array specifying the wavenumber grid [cm-1].
        c_foreign: Numpy array of foreign broadening coefficients.
        c_self: Numpy array of self broadening coefficients.
        t_foreign: Numpy array of foreign broadening temperature dependence factors.
        t_self: Numpy array of self broadening temperature dependence factors.

    Returns:
        Numpy array of absorption coefficients [cm2].
    """
    num_layers, num_wpoints = temperature.size, grid.size
    absorption_coefficient = zeros((num_layers, num_wpoints))
    ps = pressure[:]*volume_mixing_ratio[:]
    n = ones(num_layers)
    types = [c_uint64, c_int] + 9*[cpointer(c_double)]
    gas_optics.calc_water_vapor_ctm_optical_depth.argtypes = types
    gas_optics.calc_water_vapor_ctm_optical_depth.restype = check_return_code
    gas_optics.calc_water_vapor_ctm_optical_depth(c_uint64(num_wpoints), c_int(num_layers),
                                                  absorption_coefficient, c_self,
                                                  temperature, ps, n, t_self,
                                                  c_foreign, pressure, t_foreign)
    return absorption_coefficient


def ozone_continuum(temperature, grid, cross_section):
    """Calculates ozone continuum absorption coefficients.

    Args:
        temperature: Numpy array of temperatures [K].
        grid: Numpy array specifying the wavenumber grid [cm-1].
        cross_section: Numpy array of ozone cross sections [cm2].

    Returns:
        Numpy array of absorption coefficients [cm2].
    """
    num_layers, num_wpoints = temperature.size, grid.size
    absorption_coefficient = zeros((num_layers, num_wpoints))
    n = ones(num_layers)
    types = [c_uint64, c_int] + 3*[cpointer(c_double)]
    gas_optics.calc_ozone_ctm_optical_depth.argtypes = types
    gas_optics.calc_ozone_ctm_optical_depth.restype = check_return_code
    gas_optics.calc_ozone_ctm_optical_depth(c_uint64(num_wpoints), c_int(num_layers),
                                            cross_section, n, absorption_coefficient)
    return absorption_coefficient


def check_return_code(value):
    """Raises an exception if one of the GRTcode c functions fails.

    Args:
        value: Integer code returned from GRTcode c function.

    Raises:
        ValueError if GRTcode c function fails.

    Returns:
        Integer code returned from GRTcode c function.
    """
    if value != GRTCODE_SUCCESS:
        raise ValueError("GRTCODE c function returned error {}".format(value))
    return value
