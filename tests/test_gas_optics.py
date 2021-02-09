from numpy import arange

from pygrt.gas_optics import Gas


if __name__ == "__main__":
    gas = Gas("H2O")
    grid = arange(1., 3250., 0.01)
    k = gas.absorption_coefficient(101000., 300., 0.1, grid)
