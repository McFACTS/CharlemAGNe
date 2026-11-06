#!/usr/bin/env python3
######## Imports ########
#### Standard Library ####
from importlib import resources as impresources
import os
import time
#### Local ####
from charlemagne import data
from charlemagne._abi import abi_solve_sg03, abi_solve_tqm05
######## Setup ########
OPAC_PATH = str(impresources.files(data) / 'pagn_combined.sp5')
assert os.path.isfile(OPAC_PATH)
OPAC_NAME = os.path.split(OPAC_PATH)[-1].split(".")[0]
OPAC_ID = 1001
G_SI = 6.67e-11
C_SI = 299792458
MSUN = 1.988e30

######## Tests ########
def test_sg03_solve():
    tic = time.perf_counter()
    profile_dict = abi_solve_sg03(OPAC_NAME, OPAC_PATH, OPAC_ID)
    toc = time.perf_counter()
    print(f"SG03 time: {toc-tic:.6f}s")

def test_tqm05_solve():
    tic = time.perf_counter()
    profile_dict = abi_solve_tqm05(OPAC_NAME, OPAC_PATH, OPAC_ID)
    toc = time.perf_counter()
    print(f"TQM05 time: {toc-tic:.6f}s")

def test_rmax():
    profile_dict = abi_solve_sg03(OPAC_NAME, OPAC_PATH, OPAC_ID)
    r_cyl = profile_dict["r_cyl"]
    rmin_default, rmax_default = r_cyl.min(), r_cyl.max()
    rmin_default = rmin_default / (1.e8 * MSUN * G_SI / C_SI**2)
    rmax_default = rmax_default / (1.e8 * MSUN * G_SI / C_SI**2)
    assert abs((rmin_default - 6.) / 6.) < 1.e-1
    assert abs((rmax_default - 1.e7) / 1.e7) < 1.e-1
    print(rmin_default)
    print(rmax_default)
    profile_dict2 = abi_solve_sg03(
        OPAC_NAME, OPAC_PATH, OPAC_ID,
        rmax=rmax_default*2
    )
    r_cyl2 = profile_dict2["r_cyl"]
    rmin2, rmax2 = r_cyl2.min(), r_cyl2.max()
    rmin2 = rmin2 / (1.e8 * MSUN * G_SI / C_SI**2)
    rmax2 = rmax2 / (1.e8 * MSUN * G_SI / C_SI**2)
    print(rmin2)
    print(rmax2)
    assert abs((rmin2 - 6.) / 6.) < 1.e-1
    assert abs((rmax2 - 2.e7) / 2.e7) < 1.e-1

######## Execution ########
if __name__ == "__main__":
    test_sg03_solve()
    test_tqm05_solve()
    test_rmax()
