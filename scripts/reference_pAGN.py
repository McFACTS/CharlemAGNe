#!/usr/bin/env python3
"""Test against the existing implementation in pAGN
"""
######## Imports ########
#### Standard Library ####
#### Third Party ####
import numpy as np
import scipy
from astropy import units as u
#### pAGN ####
from pagn import Sirko

######## Setup ########
SG03_PAGN_KWARGS = {
    "Mbh"   : 1.0e8 * u.solMass.to(u.kg),
    "alpha" : 0.01,
    "le"    : 0.5,
    "eps"   : 0.1,
}
SG03_GRID_PTS = 1e4

######## Functions ########
def make_SG03():
    disk = Sirko.SirkoAGN(**SG03_PAGN_KWARGS)
    disk.solve_disk(SG03_GRID_PTS)
    return disk
    

######## Main ########
def main():
    SG03 = make_SG03()
    print(f"Mdot: {SG03.Mdot}")
    print(f"Mdot: {SG03.Mdot * (u.kg.to(u.solMass) * u.yr.to(u.s))} msun / yr")
    return

######## Execution ########
if __name__ == "__main__":
    main()

