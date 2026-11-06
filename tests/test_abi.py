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

######## Execution ########
if __name__ == "__main__":
    test_sg03_solve()
    test_tqm05_solve()
