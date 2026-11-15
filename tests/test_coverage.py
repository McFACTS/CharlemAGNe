#!/usr/bin/env python3
######## Imports ########
#### Standard Library ####
import collections
import contextlib
import itertools
import os
import time
#### Third Party ####
import numpy as np
#### Local ####
from charlemagne.disk import SG03, TQM05, SG03_DEFAULTS, TQM05_DEFAULTS

######## Setup ########
# Taken from <https://stackoverflow.com/a/9098295/4761692>
def named_product(**items):
    Options = collections.namedtuple('Options', items.keys())
    return itertools.starmap(Options, itertools.product(*items.values()))

ALL_SMBH_MASS = np.geomspace(1.e5, 1.e10, 16)
ALL_TQM_MASS = np.geomspace(1.e6, 1.e10, 13)
ALL_EDDINGTON = np.linspace(0.1, 1.0, 4)
NUM_SIM = ALL_SMBH_MASS.size * ALL_EDDINGTON.size
PROFILE_FINITE = ["omega", "Teff", "Mdot", "T", "cs", "rho", "kappa", "tauV", "h", "dr_cyl", "dA", "dV", "mass_disk", "Mdot_sf"]

######## Tests ########
def test_sg03_coverage():
    """Really just see if it crashes"""
    print(f"Running {NUM_SIM} SG03 disks!")
    tic = time.perf_counter()
    for config in named_product(
        smbh_mass = ALL_SMBH_MASS,
        eddington = ALL_EDDINGTON,
    ):
        # Construct a disk
        with open(os.devnull, 'w') as devnull:
          with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
            sg03 = SG03(
                smbh_mass=config.smbh_mass,
                eddington_ratio=config.eddington,
            )
        # Check that there's no NaNs or infs
        for q in PROFILE_FINITE:
            try:
                assert np.all(np.isfinite(sg03[q]))
            except AssertionError as exc:
                print(f"smbh_mass: {config.smbh_mass}; lum: {config.eddington}")
                raise
    toc = time.perf_counter()
    print(f"SG03 time: {toc-tic:.6f}s")

def test_tqm05_coverage():
    print(f"Running {NUM_SIM} TQM05 disks!")
    tic = time.perf_counter()
    for config in named_product(
            smbh_mass = ALL_TQM_MASS,
            eddington = ALL_EDDINGTON,
        ):
        try:
          with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
                tqm05 = TQM05(
                    smbh_mass=config.smbh_mass,
                    #eddington_ratio=config.eddington,
                )
        except Exception as exc:
            print(f"Crash: smbh_mass: {config.smbh_mass}; lum: {config.eddington}")
            raise
        print(f"smbh_mass: {config.smbh_mass}; lum: {config.eddington}")
        # Check that there's no NaNs or infs
        for q in PROFILE_FINITE:
            try:
                assert np.all(np.isfinite(tqm05[q]))
            except AssertionError as exc:
                print(f"smbh_mass: {config.smbh_mass}; lum: {config.eddington}")
                raise
    toc = time.perf_counter()
    print(f"TQM05 time: {toc-tic:.6f}s")

######## Execution ########
if __name__ == "__main__":
    test_sg03_coverage()
    test_tqm05_coverage()
