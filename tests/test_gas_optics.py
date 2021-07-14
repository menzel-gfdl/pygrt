from netCDF4 import Dataset
from numpy import arange, asarray

from pygrt.gas_optics import Gas


if __name__ == "__main__":
    gas = Gas("H2O")
    grid = arange(1., 3000., 0.1)
    temperature = asarray([230.92, 269.01, 203.37, 288.99])
    pressure = asarray([10.0, 117.0, 11419.0, 98388.0])
    vmr = asarray([4.072945e-06, 5.244536e-06, 3.039952e-06, 0.006637074])
    with Dataset("spectra.nc", "w") as dataset:
        dataset.createDimension("wavenumber", grid.size)
        dataset.createDimension("pressure", pressure.size)
        for name, units, data in zip(["wavenumber", "pressure", "temperature", "vmr"],
                                     ["cm-1", "Pa", "K", "mol mol-1"],
                                     [grid, pressure, temperature, vmr]):
            dim = "wavenumber" if name == "wavenumber" else "pressure"
            v = dataset.createVariable(name, "f8", (dim,))
            v.setncattr("units", units)
            v[:] = data
        v = dataset.createVariable("absorption_coefficient", "f8",
                                   ("pressure", "wavenumber"))
        v.setncattr("units", "m2")
        for i in range(pressure.size):
            v[i, :] = gas.absorption_coefficient(temperature[i], pressure[i], vmr[i], grid)
