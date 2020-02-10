"""Makes reactions in chemical kinetic model all irreversible.
"""
# Standard libraries
import copy
import math
import warnings
from itertools import repeat

try:
    import numpy as np
except ImportError:
    print('Error: NumPy must be installed.')
    raise

import pint

# Local imports
from . import chem_utilities as chem
from . import mech_interpret as mech
from .irrev_mech import write_mech


def remove_plog(
    mech_name, therm_name=None, pressure='1.0 atm', output_file="mech.txt"
    ):
    """Remove PLOG reactions from Chemkin-style model based on desired pressure.

    This utility is very simple: given an input pressure, all PLOG reactions
    are converted to pressure-independent reactions, using the Arrhenius parameters
    from the nearest pressure.

    Parameters
    ----------
    mech_name: str
        Kinetic model filename (e.g. 'mech.dat')
    therm_name: str, optional
        Thermodynamic database filename (e.g. 'therm.dat') or None
    pressure: str, optional
        Pressure with units for hard-coding PLOG reactions. 
        If no units specified, atm assumed.
    output_file: str, optional
        Filename for new kinetic model

    Returns
    -------
    None

    """
    ureg = pint.UnitRegistry()
    pressure = ureg.Quantity(pressure)
    try:
        pressure.ito(ureg.pascal)
    except pint.DimensionalityError:
        warnings.warn(
            'No units specified, or units incompatible with pressure. ',
            'Assuming atm.'
            )
        pressure = (pressure.magnitude * ureg.atm).to(ureg.pascal)

    # interpret reaction mechanism file
    [elems, specs, reacs] = mech.read_mech(mech_name, therm_name)

    for idx, rxn in enumerate(reacs):
        if rxn.plog:
            pars_closest = rxn.plog_par[0]
            for pars in rxn.plog_par:
                if np.abs(pressure - pars[0]) < np.abs(pressure - pars_closest):
                    pars_closest = pars

            rxn.A, rxn.B, rxn.E = pars_closest[1], pars_closest[2], pars_closest[3]
            rxn.plog = False
            rxn.plog_par = None

            reacs[idx] = rxn

    # write new reaction list to new file
    write_mech(output_file, elems, specs, reacs)
