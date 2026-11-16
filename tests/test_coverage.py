#!/usr/bin/env python3
######## Imports ########
#### Standard Library ####
import collections
import itertools
from contextlib import redirect_stdout, redirect_stderr
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
ALL_EDDINGTON = np.geomspace(0.01, 10.0, 13)
NUM_SIM = ALL_SMBH_MASS.size * ALL_EDDINGTON.size
PROFILE_FINITE = [
    "omega", "Teff", "Mdot", "T", "cs", "rho",
    "kappa", "tauV", "h", "dr_cyl", "dA", "dV", "mass_disk", "Mdot_sf"
]
PASS_THRESHOLD = 10.0
NPTS = int(1e4)

######## Tests ########
def test_sg03_coverage():
    """Really just see if it crashes"""
    print(f"Running {NUM_SIM} SG03 disks!")
    tic = time.perf_counter()
    for config in named_product(
        smbh_mass = ALL_SMBH_MASS,
        eddington = ALL_EDDINGTON,
    ):
        with open(os.devnull, 'w') as devnull:
          with redirect_stdout(devnull), redirect_stderr(devnull):
            # Construct a disk
            sg03 = SG03(
                smbh_mass=config.smbh_mass,
                eddington_ratio=config.eddington,
                npts=NPTS,
            )
        print(
            f"SG03 -- smbh_mass: {config.smbh_mass:.4e}; "
            f"lum: {config.eddington:.4e}; "
            f"max_residual: {sg03.max_residual:.3f}"
        )
        # Check that there's no NaNs or infs
        for q in PROFILE_FINITE:
            try:
                assert np.all(np.isfinite(sg03[q]))
            except AssertionError as exc:
                raise
    toc = time.perf_counter()
    print(f"SG03 time: {toc-tic:.6f}s")

def test_tqm05_coverage():
    print(f"Running {NUM_SIM} TQM05 disks!")
    n_fail = 0
    tic = time.perf_counter()
    for config in named_product(
            smbh_mass = ALL_SMBH_MASS,
            eddington = ALL_EDDINGTON,
        ):
        try:
          with open(os.devnull, 'w') as devnull:
            with redirect_stdout(devnull), redirect_stderr(devnull):
                tqm05 = TQM05(
                    smbh_mass=config.smbh_mass,
                    eddington_ratio=config.eddington,
                    npts=NPTS,
                )
        except Exception as exc:
            print(f"Crash: smbh_mass: {config.smbh_mass}; lum: {config.eddington}")
            raise
        print(
            f"TQM05 -- smbh_mass: {config.smbh_mass:.4e}; "
            f"lum: {config.eddington:.4e}; "
            f"max_residual: {tqm05.max_residual:.3f}"
        )
        if tqm05.max_residual > PASS_THRESHOLD:
            n_fail += 1
        # Check that there's no NaNs or infs
        for q in PROFILE_FINITE:
            try:
                assert np.all(np.isfinite(tqm05[q]))
            except AssertionError as exc:
                print(f"smbh_mass: {config.smbh_mass}; lum: {config.eddington}")
                raise
    toc = time.perf_counter()
    print(f"TQM05 time: {toc-tic:.6f}s ({n_fail} sims failed)")

######## Execution ########
if __name__ == "__main__":
    test_sg03_coverage()
    test_tqm05_coverage()
