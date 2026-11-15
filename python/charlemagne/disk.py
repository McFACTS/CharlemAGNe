"""A one-dimensional disk wrapper in Python"""
######## Imports ########
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

class Disk1D(object):
    """This object is the pure Python representation of the disk module"""

class SG03(Disk1D):
    """The Standard Sirko-Goodman disk with constant mdot throughout the disk
    """
    def __init__(
        self,
        rmin                    = 6.0,
        rmax                    = 1.e07,
        npts                    = int(1e4),
        smbh_mass               = 1.0e8,
        eddington_ratio         = 0.5,
        alpha                   = 0.01,
        X_mass_fraction         = 0.7,
        radiative_efficiency    = 0.1,
        viscosity_switch        = False,
    ):
        """Initialize the SG03 Disk object

        Parameters
        ----------
        """

        # Define the parameters dictionary
        self.params = {
        rmin                    : rmin,
        rmax                    : rmax,
        npts                    : npts,
        smbh_mass               : smbh_mass,
        eddington_ratio         : eddington_ratio,
        alpha                   : alpha,
        X_mass_fraction         : X_mass_fraction,
        radiative_efficiency    : radiative_efficiency,
        viscosity_switch        : viscosity_switch,
    }

    def _solve_disk(self):
        """Call the C++ library and solve the disk"""
        return abi_solve_sg03(
            OPAC_NAME, OPAC_PATH, OPAC_ID,
            **self.params
        )
