#!/usr/bin/env python3
"""Construct a spiner-compatible table from the pAGN opacity tables

The pAGN team put a lot of work into developing a set of opacity tables
  which are valid Rosseland Mean Opacity estimates from 5 to 10^8 Kelvin.
They are not evenly spaced, which is a requirement for singularity-opac.

At low temperatures, pAGN uses the Semenov+2003 opacities for ice and dust.
https://arxiv.org/pdf/astro-ph/0308344
At high temperatures, pAGN uses Badnell+2005 
https://academic.oup.com/mnras/article/360/2/458/964515?login=false

The Semenov opacities were provided (rather than in a csv format) in the body
  of Fortran 77 subroutines.
The pAGN developers created an F2py module to provide a Python interface,
  which lives only in the 'utils' directory and is not installed through
  PyPi when you `pip install pagn`.
The size of the grid pAGN used for the Semenov opacities was arbitrary, and
  one could modify the spacing by rerunning the module.

The Badnell opacities have an even spacing in logT below and above 10^6 K,
  and are interpolated directly from published tables.

The pAGN developers claim that "for the smoothest transition" between the
  two opacity models, the highest temperature recorded from Semenov is
  at 4,544 Kelvin and the lowest temperature recorded from Badnell is at
  7079 Kelvin.
I feel more comfortable with this gap knowing that they looked at it and
  decided this was best after trying different things.

In order to take advantage of their hard work, while using singularity-opac,
  we need to repackage it in singularity-opac's format, which means sampling
  on an evenly spaced grid.
The pAGN developers use a Bivariate Quadratic Spline 
  through scipy.interpolate's RectBivariateSpline object,
  assigning kx=2, ky=2, and s=0.
As a result, it passes through every point exactly,
  and the function and its first derivative are continuous,
  but its second derivative may change abruptly.
When we re-evaluate, we should use the same interpolator.
"""
######## Imports ########
#### Standard Library ####
import os
from importlib import resources as impresources
#### Third Party ####
import numpy as np
import h5py
import scipy
#### pAGN ####
import pagn
from pagn.opacities import get_combined_opacity

######## Setup ########
NAME    = "pagn_combined"
OUT     = f"{NAME}.sp5"
OPACID  = 1001

## Temperature grid
LOG10T_MIN  = np.log10(5)
LOG10T_MAX  = 8.
LOG10T_NUM  = 147 # 147 for 0.05 spacing or 244 for 0.03 spacing
LOG10T_GRID = np.linspace(LOG10T_MIN, LOG10T_MAX, LOG10T_NUM)

## Density grid
LOG10R_MIN    = -7.5
LOG10R_MAX    = 1.0
LOG10R_NUM    = 18
LOG10R_GRID = np.linspace(LOG10R_MIN, LOG10R_MAX, LOG10R_NUM)

#### Load the pAGN data ####
#print(os.listdir(impresources.files(pagn.opacities) / "opacity_tables" ))
opacity_dir = impresources.files(pagn.opacities) / "opacity_tables"
PAGN_LOGRT_COMBINED_DATA = np.loadtxt(opacity_dir / "logRT_combined_table.dat")
PAGN_LOGRT_COMBINED_T = np.loadtxt(opacity_dir / "logRT_combined_Taxis.dat")
PAGN_LOGRT_COMBINED_R = np.loadtxt(opacity_dir / "logRT_combined_Raxis.dat")
_R, _T = np.meshgrid(PAGN_LOGRT_COMBINED_R, PAGN_LOGRT_COMBINED_T)
CEIL = 9.999
CEIL_MASK = np.isclose(PAGN_LOGRT_COMBINED_DATA, CEIL)
APPLY_MASK = False

## Trim data
if APPLY_MASK:
    _R = _R[~CEIL_MASK]
    _T = _T[~CEIL_MASK]
    _DATA = PAGN_LOGRT_COMBINED_DATA[~CEIL_MASK]
else:
    _R = _R.flatten()
    _T = _T.flatten()
    _DATA = PAGN_LOGRT_COMBINED_DATA.flatten()
PAGN_INTERP = scipy.interpolate.RectBivariateSpline(
    PAGN_LOGRT_COMBINED_T,#_T,
    PAGN_LOGRT_COMBINED_R,#_R,
    PAGN_LOGRT_COMBINED_DATA,#_DATA,
    kx=2,
    ky=2,
    s=0,
)
#PAGN_INTERP = get_combined_opacity()

######## Main ########
def main():
    # Generate the meshgrid
    log10T, log10R = np.meshgrid(LOG10T_GRID, LOG10R_GRID, indexing="ij")
    # pAGN's interpolator goes in rho[g/cm^3] / (T/1_000_000 K)^3
    # singularity-opac wants CGS
    #log10Ropac = log10R - 3*log10T + 21
    # Evaluate the meshgrid
    log10Kappa = PAGN_INTERP(log10T, log10R, grid=False)
    # Get linear kappa in cgs (cm^2 / g from m^2 / kg)
    Kappa = 10**(log10Kappa + 1.0)
    with h5py.File(OUT, 'w') as F:
        # Create a top level group with the name of the model
        group = F.create_group(NAME)
        # Assign attributes
        group.attrs["opacid"] = np.int32(OPACID)
        group.attrs["opac_name"] = NAME
        # Instantiate the array to write to SPS
        kappa_sps = np.empty((Kappa.shape[1], Kappa.shape[0], 1))
        kappa_sps[:,:,0] = Kappa.T

        ## Write opacity group ##
        sub_name = f"{NAME}/Rosseland group opacity"
        sub = F.create_group(sub_name)
        sub.attrs["rank"] = np.int32(3)
        sub.attrs["dims"] = np.asarray(list(kappa_sps.shape),dtype=np.int32)[::-1]
        sub.attrs["index_types"] = np.array([2, 0, 0])
        sub.attrs["index_types_info"] = "Interpolated:0\nNamed:1\nIndexed:2"
        sub.create_dataset("data", data=kappa_sps)
        grid_group = sub.create_group("grids")
        grid_group.create_dataset("grid_[2]", data=np.array([
            LOG10T_GRID[0],
            LOG10T_GRID[-1],
            LOG10T_GRID[1] - LOG10T_GRID[0]
        ]))
        grid_group["grid_[2]"].attrs["npoints"] = np.int32(LOG10T_GRID.size)
        grid_group["grid_[2]"].attrs["range columns"] = "[0]:min [1]:max [2]:dx"
        grid_group.create_dataset("grid_[3]", data=np.array([
            LOG10R_GRID[0],
            LOG10R_GRID[-1],
            LOG10R_GRID[1] - LOG10R_GRID[0]
        ]))
        grid_group["grid_[3]"].attrs["npoints"] = np.int32(LOG10R_GRID.size)
        grid_group["grid_[3]"].attrs["range columns"] = "[0]:min [1]:max [2]:dx"

        ## Write bounds group ##
        sub_name = f"{NAME}/group bounds"
        sub = F.create_group(sub_name)
        sub.attrs["rank"] = np.int32(1)
        sub.attrs["dims"] = np.asarray([2,], dtype=np.int32)
        sub.attrs["index_types"] = np.array([2])
        sub.attrs["index_types_info"] = "Interpolated:0\nNamed:1\nIndexed:2"
        sub.create_dataset("data", data=np.asarray([0.0, np.inf]))
        # Create empty grid group for some reason
        grid_group = sub.create_group("grids")

######## Execution ########
if __name__ == "__main__":
    main()

