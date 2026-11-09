#!/usr/bin/env python3
######## Imports ########
#### Standard Library ####
from importlib import resources as impresources
import os
import time
#### Third Party ####
from matplotlib import pyplot as plt
from pagn import Sirko, Thompson
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

SG03_DEFAULTS = {
    "rmin" : 12.0,
    "rmax" : 2.e7,
    "npts" : int(1e4),
    "smbh_mass" : 1.e8,
    "eddington_ratio" : 0.5,
    "alpha" : 0.01,
    "X_mass_fraction" : 0.7,
    "radiative_efficiency" : 0.1,
    "viscosity_switch" : False,
}

SG03_pAGN_ALIAS = {
    "smbh_mass" : "Mbh",
    "eddington_ratio" : "le",
    "alpha" : "alpha",
    "X_mass_fraction" : "X",
    "radiative_efficiency" : "eps",
    "viscosity_switch" : "b",
}

PLOT_LABELS = {
    "rho"   : r"Density ($\rho$)",
    "T"     : r"Temperature ($T$)",
    "kappa" : r"Opacity ($\kappa$)",
    "tauV"  : r"Optical Depth ($\tau_{V}$)",
    "h"     : r"Aspect Ratio ($h$)",
}

######## Tests ########
def make_disks(params):
    pagn_params = {}
    for key in params:
        if key in SG03_pAGN_ALIAS:
            pagn_params[SG03_pAGN_ALIAS[key]] = params[key]
    pagn_params["Mbh"] *= MSUN
    tic = time.perf_counter()
    pagn_disk = Sirko.SirkoAGN(**pagn_params)
    pagn_disk.solve_disk(N=params["npts"])
    toc = time.perf_counter()
    print(f"pAGN time: {toc-tic:.6f}s")
    tic = time.perf_counter()
    charlemagne_dict = abi_solve_sg03(OPAC_NAME, OPAC_PATH, OPAC_ID, **params)
    toc = time.perf_counter()
    print(f"CharlemAGNe time: {toc-tic:.6f}s")

    # Get profile things
    x_pagn = pagn_disk.R / (pagn_disk.Rs / 2)
    pagn_surf_dens = 2. * pagn_disk.h * pagn_disk.rho
    pagn_kappa = 2. * pagn_surf_dens * pagn_disk.h
    x_char = charlemagne_dict["r_cyl"] / (1.e8 * MSUN * G_SI / C_SI**2)
    print(list(charlemagne_dict.keys()))

    # Make plots
    plt.style.use('bmh')
    for item, label in PLOT_LABELS.items():
        # Density
        fig, ax = plt.subplots()
        fig.suptitle(label)
        ax.loglog(x_pagn, getattr(pagn_disk, item), label="pAGN")
        ax.loglog(x_char, charlemagne_dict[item], label="CharlemAGNe", linestyle="dashed")
        ax.legend()
        plt.savefig(f"SG03_{item}.jpg")
        plt.close()

######## Execution ########
if __name__ == "__main__":
    make_disks(SG03_DEFAULTS)
