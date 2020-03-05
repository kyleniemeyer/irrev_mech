"""Makes reactions in chemical kinetic model all irreversible.
"""
# Standard libraries
import logging

import numpy as np
from pint import DimensionalityError

# Local imports
from .chem_utilities import units, Q_, ReacInfo
from .mech_interpret import read_mech
from .irrev_mech import write_mech


def find_nearest_parameters(pressure, parameter_list):
    """Finds closest set of rate parameters based on pressure.

    Parameters
    ----------
    pressure: pint.Quantity
        Target pressure
    parameter_list: float
        List of pressure and associated Arrhenius parameters

    Returns
    -------
    pars_closest : list of float
        Arrhenius parameters

    """
    pars_closest = parameter_list[0]
    # internally PLOG pressure given in Pa
    for pars in parameter_list:
        if np.abs(pressure - pars[0]) < np.abs(pressure - pars_closest[0]):
            pars_closest = pars
            
    return pars_closest[1:]


def interpolate_parameters(pressure, parameter_list):
    """Interpolates rate parameters based on pressure.

    Parameters
    ----------
    pressure: pint.Quantity
        Target pressure
    parameter_list: float
        List of pressure and associated Arrhenius parameters

    Returns
    -------
    parameters : list of float
        Pressure and Arrhenius parameters

    """
    if pressure < parameter_list[0][0]:
        logging.info(
            'Warning: pressure below lowest given in PLOG reaction. '
            'Using nearest set of values.'
            )
        parameters = parameter_list[0][1:]
    elif pressure > parameter_list[-1][0]:
        logging.info(
            'Warning: pressure above highest given in PLOG reaction. '
            'Using nearest set of values.'
            )
        parameters = parameter_list[-1][1:]
    else:
        pars_under = None
        pars_over = None
        for idx, pars in enumerate(parameter_list):
            if pressure >= pars[0] and pressure < parameter_list[idx + 1][0]:
                pars_under = parameter_list[idx]
                pars_over = parameter_list[idx + 1]
                break

        pressure_log = (
            np.log(pressure / pars_under[0]) / np.log(pars_over[0] / pars_under[0])
            ).magnitude
        parameters = [0.0, 0.0, 0.0]
        parameters[0] = (
            np.power(pars_under[1], 1.0 - pressure_log) * 
            np.power(pars_over[1], pressure_log)
            )
        parameters[1] = pars_under[2] * (1.0 - pressure_log) + pars_over[2] * pressure_log
        parameters[2] = pars_under[3] * (1.0 - pressure_log) + pars_over[3] * pressure_log
            
    return parameters


def convert_plog(rxn, pressure, method='interpolate'):
    """Convert PLOG reaction to elementary Arrhenius reaction based on nearest pressure.

    Parameters
    ----------
    rxn : ReacInfo
        PLOG reaction to be converted
    pressure: pint.Quantity
        Target pressure
    method: str
        Method for choosing new rate parameters, either "interpolate" or "nearest"

    Returns
    -------
    rxn : ReacInfo
        Converted reaction

    """
    assert rxn.plog
    assert method in ['interpolate', 'nearest']

    if method == 'nearest':
        rate_parameters = find_nearest_parameters(pressure, rxn.plog_par)
    elif method == 'interpolate':
        rate_parameters = interpolate_parameters(pressure, rxn.plog_par)

    rxn = ReacInfo(
        rxn.rev, rxn.reac, rxn.reac_nu, rxn.prod, rxn.prod_nu, 
        rate_parameters[0], rate_parameters[1], rate_parameters[2]
        )

    return rxn


def remove_plog_reactions(
    mech_name, therm_name=None, pressure='1.0 atm', output_file='mech.out',
    method='interpolate'
    ):
    """Remove PLOG reactions from Chemkin-style model based on desired pressure.

    This utility converts all PLOG reactions to pressure-independent reactions, 
    using Arrhenius parameters either taken from those given at the nearest pressure
    or by interpolating between those given at the surrounding pressures.

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
    method : str, optional
        Method for de-plogging, either "interpolate" or "nearest"

    Returns
    -------
    None

    """
    pressure = Q_(pressure)
    try:
        pressure.ito('pascal')
    except DimensionalityError:
        logging.info(
            'No units specified, or units incompatible with pressure. '
            'Assuming atm.'
            )
        pressure = (pressure.magnitude * units.atm).to('pascal')

    assert method in ['interpolate', 'nearest']

    # interpret reaction mechanism file
    [elems, specs, reacs] = read_mech(mech_name, therm_name)

    for idx, rxn in enumerate(reacs):
        if rxn.plog:
            rxn = convert_plog(rxn, pressure)
            reacs[idx] = rxn

    # write new reaction list to new file
    write_mech(output_file, elems, specs, reacs)
