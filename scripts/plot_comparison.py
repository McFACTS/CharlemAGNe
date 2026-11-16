#!/usr/bin/env python3
######## Imports ########
#### Standard Library ####
import os
import time
#### Third Party ####
from matplotlib import pyplot as plt
from pagn import Sirko, Thompson
#### Local ####
from charlemagne.disk import SG03, TQM05

######## Setup ########
G_SI = 6.67e-11
C_SI = 299792458
MSUN = 1.988e30

SG03_DEFAULTS = {
    "rmin"                  : 12.0,     # 12.0 Rg
    "rmax"                  : 2.e7,     # 2.e7 Rg
    "npts"                  : int(1e4),
    "smbh_mass"             : 1.e8,     # MSUN
    "eddington_ratio"       : 0.5,
    "alpha"                 : 0.01,
    "X_mass_fraction"       : 0.7,
    "radiative_efficiency"  : 0.1,
    "viscosity_switch"      : False,
}

TQM05_DEFAULTS = {
    "rmin"                  : 6.0,          # 6.0 Rg
    "rmax"                  : 4.178e+7,     # 200 pc for 1.e8 smbh in Rg
    "npts"                  : int(1e4),
    "smbh_mass"             : 1.e8,         # MSUN
    "eddington_ratio"       : 1.23271,
    "velocity_dispersion"   : 0.0,
    "radiative_efficiency"  : 1.e-3,
    "transport_efficiency"  : 0.2,
    "supernova_pressure"    : 1.0,
    "X_mass_fraction"       : 0.7,
    "radius_inner"          : 6.0,          # 6.0 Rg
    "radius_outer"          : 4.178e+7,     # 200 pc for 1.e8 smbh in Rg
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
    charles_disk = SG03(**params)
    toc = time.perf_counter()
    print(f"CharlemAGNe time: {toc-tic:.6f}s")

    # Get profile things
    x_pagn = pagn_disk.R / (pagn_disk.Rs / 2)
    pagn_surf_dens = 2. * pagn_disk.h * pagn_disk.rho
    pagn_kappa = 2. * pagn_surf_dens * pagn_disk.h
    x_char = charles_disk["r_cyl"] / (1.e8 * MSUN * G_SI / C_SI**2)
    print(charles_disk.fields)

    # Make plots
    plt.style.use('bmh')
    for item, label in PLOT_LABELS.items():
        # Density
        fig, ax = plt.subplots()
        fig.suptitle(label)
        ax.loglog(x_pagn, getattr(pagn_disk, item), label="pAGN")
        ax.loglog(x_char, charles_disk[item], label="CharlemAGNe", linestyle="dashed")
        ax.legend()
        plt.savefig(f"SG03_{item}.jpg")
        plt.close()

######## Execution ########
if __name__ == "__main__":
    make_disks(SG03_DEFAULTS)
