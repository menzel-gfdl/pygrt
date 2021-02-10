from netCDF4 import Dataset
from numpy import arange

from pygrt.gas_optics import Gas


if __name__ == "__main__":
    gas = Gas("H2O")
    grid = arange(1., 3250., 0.01)
    pressure, temperature, vmr = 100*1000., 300., 0.006637074
    k = gas.absorption_coefficient(pressure, temperature, vmr, grid)
    with Dataset("spectra.nc", "w") as dataset:
        dataset.createDimension("wavenumber", grid.size)
        v = dataset.createVariable("wavenumber", float, ("wavenumber",))
        v.setncattr("units", "cm-1")
        v[:] = grid
        v = dataset.createVariable("absorption_coefficient", float, ("wavenumber",))
        v.setncattr("units", "m2")
        v[:] = k[:]
