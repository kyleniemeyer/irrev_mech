"""Makes reactions in chemical kinetic model all irreversible.
"""
# Standard libraries
import warnings

import numpy as np
from pint import DimensionalityError

# Local imports
from .chem_utilities import units, Q_, ReacInfo
from .mech_interpret import read_mech
from .irrev_mech import write_mech


def convert_plog(rxn, pressure):
    """Convert PLOG reaction to elementary Arrhenius reaction based on nearest pressure.

    Parameters
    ----------
    rxn : ReacInfo
        PLOG reaction to be converted
    pressure: pint.Quantity
        Target pressure

    Returns
    -------
    rxn : ReacInfo
        Converted reaction

    """
    assert rxn.plog

    pars_closest = rxn.plog_par[0]
    # internally PLOG pressure given in Pa
    for pars in rxn.plog_par:
        if np.abs(pressure - pars[0]) < np.abs(pressure - pars_closest[0]):
            pars_closest = pars

    rxn = ReacInfo(
        rxn.rev, rxn.reac, rxn.reac_nu, rxn.prod, rxn.prod_nu, 
        pars_closest[1], pars_closest[2], pars_closest[3]
        )

    return rxn


def remove_plog_reactions(
    mech_name, therm_name=None, pressure='1.0 atm', output_file='mech.out'
    ):
    """Remove PLOG reactions from Chemkin-style model based on desired pressure.

    This utility is very simple: given an input pressure, all PLOG reactions
    are converted to pressure-independent reactions, using the Arrhenius parameters
    from the nearest pressure.

    Parameters
    ----------
    mech_name : str
        Kinetic model filename (e.g. 'mech.dat')
    therm_name : str, optional
        Thermodynamic database filename (e.g. 'therm.dat') or None
    pressure : str, optional
        Pressure with units for hard-coding PLOG reactions. 
        If no units specified, atm assumed.
    output_file : str, optional
        Filename for new kinetic model

    Returns
    -------
    None

    """
    pressure = Q_(pressure)
    try:
        pressure.ito('pascal')
    except DimensionalityError:
        warnings.warn(
            'No units specified, or units incompatible with pressure. ',
            'Assuming atm.'
            )
        pressure = (pressure.magnitude * units.atm).to('pascal')

    # interpret reaction mechanism file
    [elems, specs, reacs] = read_mech(mech_name, therm_name)

    for idx, rxn in enumerate(reacs):
        if rxn.plog:
            rxn = convert_plog(rxn, pressure)
            reacs[idx] = rxn

    # write new reaction list to new file
    write_mech(output_file, elems, specs, reacs)
