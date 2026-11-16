"""A one-dimensional disk wrapper in Python"""
######## Imports ########
#### Standard Library ####
from abc import ABC, abstractmethod
from importlib import resources as impresources
import os
#### Local ####
from charlemagne import data
from charlemagne._abi import abi_solve_sg03, abi_solve_tqm05

######## Setup ########
OPAC_PATH = str(impresources.files(data) / 'pagn_combined.sp5')
assert os.path.isfile(OPAC_PATH)
OPAC_NAME = os.path.split(OPAC_PATH)[-1].split(".")[0]
OPAC_ID = 1001

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


class Disk1D(ABC):
    """This object is the pure Python representation of the disk module"""

    @abstractmethod
    def _solve_disk(self):
        raise NotImplementedError

    def __getattr__(self, key):
        """Yeah, I'm lazy and don't want to define a method for all these
          different params

        Parameters
        ----------
        key : str
            The thing we want from params
        """
        # Check if it's in params
        if key in self.params:
            return self.params[key]
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute: {key}"
        )

    def __getitem__(self, key):
        """Okay. This one is for profiles from the solve method
        
        Parameters
        ----------
        key : str
            The name of the profile we want to load (E.g. r_cyl, kappa, etc...)

        Returns
        -------
        np.ndarray
            The profile values
        """
        return self.profile[key]

    @property
    def fields(self):
        """Return a list of available profile fields

        Returns
        -------
        list
            The profile dictionary keys available
        """
        return list(self.profile.keys())

    @property
    def max_residual(self):
        return self.profile["max_residual"].max()

class SG03(Disk1D):
    """The Standard Sirko-Goodman disk with constant mdot throughout the disk
    """
    def __init__(
        self,
        rmin                  = SG03_DEFAULTS["rmin"],
        rmax                  = SG03_DEFAULTS["rmax"],
        npts                  = SG03_DEFAULTS["npts"],
        smbh_mass             = SG03_DEFAULTS["smbh_mass"],
        eddington_ratio       = SG03_DEFAULTS["eddington_ratio"],
        alpha                 = SG03_DEFAULTS["alpha"],
        X_mass_fraction       = SG03_DEFAULTS["X_mass_fraction"],
        radiative_efficiency  = SG03_DEFAULTS["radiative_efficiency"],
        viscosity_switch      = SG03_DEFAULTS["viscosity_switch"],
    ):
        """Initialize the SG03 Disk object

        Parameters
        ----------
        """

        # Define the parameters dictionary
        self.params = {
            "rmin"                  : rmin,
            "rmax"                  : rmax,
            "npts"                  : npts,
            "smbh_mass"             : smbh_mass,
            "eddington_ratio"       : eddington_ratio,
            "alpha"                 : alpha,
            "X_mass_fraction"       : X_mass_fraction,
            "radiative_efficiency"  : radiative_efficiency,
            "viscosity_switch"      : viscosity_switch,
        }
        # Solve the disk
        self.profile = self._solve_disk()

    def _solve_disk(self):
        """Call the C++ library and solve the disk"""
        return abi_solve_sg03(
            OPAC_NAME, OPAC_PATH, OPAC_ID,
            rmin                    = self.rmin,
            rmax                    = self.rmax,
            npts                    = self.npts,
            smbh_mass               = self.smbh_mass,
            eddington_ratio         = self.eddington_ratio,
            alpha                   = self.alpha,
            X_mass_fraction         = self.X_mass_fraction,
            radiative_efficiency    = self.radiative_efficiency,
            viscosity_switch        = self.viscosity_switch,
        )

class TQM05(Disk1D):
    """The TQM05 disk integrating inward
    """
    def __init__(
        self,
        rmin                    = TQM05_DEFAULTS["rmin"],
        rmax                    = TQM05_DEFAULTS["rmax"],
        npts                    = TQM05_DEFAULTS["npts"],
        smbh_mass               = TQM05_DEFAULTS["smbh_mass"],
        eddington_ratio         = TQM05_DEFAULTS["eddington_ratio"],
        velocity_dispersion     = TQM05_DEFAULTS["velocity_dispersion"],
        radiative_efficiency    = TQM05_DEFAULTS["radiative_efficiency"],
        transport_efficiency    = TQM05_DEFAULTS["transport_efficiency"],
        supernova_pressure      = TQM05_DEFAULTS["supernova_pressure"],
        X_mass_fraction         = TQM05_DEFAULTS["X_mass_fraction"],
        radius_inner            = TQM05_DEFAULTS["radius_inner"],
        radius_outer            = TQM05_DEFAULTS["radius_outer"],
    ):
        """Initialize the TQM05 Disk object

        Parameters
        ----------
        """
        # Define the parameters dictionary
        self.params = {
            "rmin"                  : rmin,
            "rmax"                  : rmax,
            "npts"                  : npts,
            "smbh_mass"             : smbh_mass,
            "eddington_ratio"       : eddington_ratio,
            "velocity_dispersion"   : velocity_dispersion,
            "radiative_efficiency"  : radiative_efficiency,
            "transport_efficiency"  : transport_efficiency,
            "supernova_pressure"    : supernova_pressure,
            "X_mass_fraction"       : X_mass_fraction,
            "radius_inner"          : radius_inner,
            "radius_outer"          : radius_outer,
        }
        # Solve the disk
        self.profile = self._solve_disk()

    def _solve_disk(self):
        """Call the C++ library and solve the disk"""
        return abi_solve_tqm05(
            OPAC_NAME, OPAC_PATH, OPAC_ID,
            rmin                    = self.rmin,
            rmax                    = self.rmax,
            npts                    = self.npts,
            smbh_mass               = self.smbh_mass,
            eddington_ratio         = self.eddington_ratio,
            velocity_dispersion     = self.velocity_dispersion,
            radiative_efficiency    = self.radiative_efficiency,
            transport_efficiency    = self.transport_efficiency,
            supernova_pressure      = self.supernova_pressure,
            X_mass_fraction         = self.X_mass_fraction,
            radius_inner            = self.radius_inner,
            radius_outer            = self.radius_outer,
        )
