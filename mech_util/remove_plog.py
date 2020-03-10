"""Makes reactions in chemical kinetic model all irreversible.
"""
# Standard libraries
import logging
import warnings

import numpy as np
from pint import DimensionalityError
from scipy.optimize import least_squares

# Local imports
from .chem_utilities import units, Q_, ReacInfo, Arrhenius
from .mech_interpret import read_mech
from .write_mech import write_mech


def calc_rate_coeff(pars, T):
    """Calculate Arrhenius reaction rate coefficient.

    Parameters
    ----------
    pars: Arrhenius
        Arrhenius parameters for reaction rate coefficient
    T: float
        Temperature [K]
    
    Returns
    -------
    float
        Reaction rate coefficient
    
    """
    return pars.pre_factor * np.exp(
        pars.temp_exponent * np.log(T) - (pars.act_energy.to('kelvin').magnitude / T)
        )


def find_nearest_parameters(pressure, parameter_list):
    """Finds closest set(s) of rate parameters based on pressure.

    If one set of coefficients are given per pressure, then this returns
    the closest set. If there are duplicate coefficients for the nearest 
    pressure, then both sets are returned.

    Parameters
    ----------
    pressure: pint.Quantity
        Target pressure
    parameter_list: float
        List of pressure and associated Arrhenius parameters

    Returns
    -------
    pars_closest : list of Arrhenius
        List of one or more Arrhenius parameters sets

    """
    # find closest, while also checking for duplicates
    diff = [np.abs(pressure - pars[0]) for pars in parameter_list]
    idx_closest = diff.index(min(diff))
    dups = [
        idx for idx, d in enumerate(diff) 
        if np.allclose(d, diff[idx_closest]) and idx != idx_closest
        ]

    pars_closest = [parameter_list[idx_closest][1]] + [parameter_list[d][1] for d in dups]

    return pars_closest


def residuals(pars, rhs, T):
    """Residual for calculating pressure-log rate coefficients.

    Parameters
    ----------
    pars: list of float
        List of multiple sets of Arrhenius factors
    rhs: float
        Calculated right-hand side of pressure-log expression at temperature ``T``
    T: float
        Temperature [K]
    
    Returns
    -------
    float
    
    """
    return rhs - np.log(np.sum([
        calc_rate_coeff(
            Arrhenius(pre_factor=A, temp_exponent=b, act_energy=E*units('kelvin')), T
            ) 
        for A, b, E in zip(pars[0::3], pars[1::3], pars[2::3])
        ]))


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
    parameters : list of Arrhenius
        List of one or more Arrhenius parameters sets

    """
    # find parameters given at closest pressure, and any duplicates
    diff = [np.abs(pressure - pars[0]) for pars in parameter_list]
    idx_closest = diff.index(min(diff))
    dups = [
        idx for idx, d in enumerate(diff) 
        if np.allclose(d, diff[idx_closest]) and idx != idx_closest
        ]
    idx_closest = [idx_closest] + dups
    pressure_closest = parameter_list[idx_closest[0]][0]
    
    if np.allclose(
        pressure.to('pascal').magnitude, 
        pressure_closest.to('pascal').magnitude
        ):
        parameters = [parameter_list[idx][1] for idx in idx_closest]

    elif pressure <= min([pars[0] for pars in parameter_list]):
        logging.info(
            'Warning: pressure below lowest given in PLOG reaction. '
            'Using nearest set of values.'
            )
        parameters = [parameter_list[idx][1] for idx in idx_closest]

    elif pressure >= max([pars[0] for pars in parameter_list]):
        logging.info(
            'Warning: pressure above highest given in PLOG reaction. '
            'Using nearest set of values.'
            )
        parameters = [parameter_list[idx][1] for idx in idx_closest]

    else:
        # find parameters at next closest pressure
        diff_next = [d for d in diff if not np.allclose(d, diff[idx_closest[0]])]
        idx_next_closest = diff.index(min(diff_next))
        dups_next = [
            idx for idx, d in enumerate(diff) 
            if np.allclose(d, diff[idx_next_closest]) and idx != idx_next_closest
            ]
        idx_next_closest = [idx_next_closest] + dups_next
        
        if pressure >= pressure_closest:
            pressure_under = parameter_list[idx_closest[0]][0]
            pars_under = [parameter_list[idx][1] for idx in idx_closest]
            pressure_over = parameter_list[idx_next_closest[0]][0]
            pars_over = [parameter_list[idx][1] for idx in idx_next_closest]
        else:
            pressure_under = parameter_list[idx_next_closest[0]][0]
            pars_under = [parameter_list[idx][1] for idx in idx_next_closest]
            pressure_over = parameter_list[idx_closest[0]][0]
            pars_over = [parameter_list[idx][1] for idx in idx_closest]

        if len(idx_closest) > 1 or len(idx_next_closest) > 1:
            # When there are multiple sets of parameters at each pressure,
            # need to use nonlinear optimization scheme to find new parameters

            pressure_log = (
                np.log(pressure / pressure_under) / 
                np.log(pressure_over / pressure_under)
                ).magnitude

            # calculate right-hand side of expression based on densely sampled temperatures
            temps = np.linspace(300.0, 5000.0, 1000)
            rhs_vals = np.zeros(len(temps))
            for idx, temp in enumerate(temps):
                rhs_vals[idx] = (
                    (1.0 - pressure_log) * np.log(np.sum([
                        calc_rate_coeff(pars, temp) for pars in pars_under
                        ])) + 
                    pressure_log * np.log(np.sum([
                        calc_rate_coeff(pars, temp) for pars in pars_over
                        ]))
                    )
            
            # number of unknowns will be based on the minimum number given at the nearest
            # two pressures (i.e., if only one pressure given on either side, then just three
            # unknowns even if two pressures given at the other side)
            vals_lower = []
            for par in pars_under:
                vals_lower.append(par.pre_factor)
                vals_lower.append(par.temp_exponent)
                vals_lower.append(par.act_energy.to('kelvin').magnitude)
            vals_higher = []
            for par in pars_over:
                vals_higher.append(par.pre_factor)
                vals_higher.append(par.temp_exponent)
                vals_higher.append(par.act_energy.to('kelvin').magnitude)
            if len(vals_lower) <= len(vals_higher):
                initial_vals = vals_lower[:]
            else:
                initial_vals = vals_higher[:]

            lower_bounds = []
            upper_bounds = []
            for idx in range(len(initial_vals) // 3):
                lower_bounds += [0.0, -10.0, -1.0e5]
                upper_bounds += [np.inf, 10.0, 1.0e5]

            # Start with low number of max function evals, increase if needed.
            warnings.filterwarnings('error')
            for mx in [1000, 5000, 10000, 20000, 40000]:
                try:
                    val_lsq = least_squares(
                        residuals, initial_vals, args=(rhs_vals, temps), 
                        method='trf', max_nfev=mx, 
                        bounds=(lower_bounds, upper_bounds),
                        )
                    break
                except RuntimeWarning:
                    continue
            else:
                logging.info(
                    'Warning: minimization failed to converge at pressure ' +
                    f'{pressure}'
                    )
            parameters = [
                Arrhenius(pre_factor=A, temp_exponent=b, act_energy=E*units('kelvin'))
                for A, b, E in zip(val_lsq.x[0::3], val_lsq.x[1::3], val_lsq.x[2::3])
                ]

        else:
            pars_under = pars_under[0]
            pars_over = pars_over[0]

            pressure_log = (
                np.log(pressure / pressure_under) / np.log(pressure_over / pressure_under)
                ).magnitude
            pre_exponential = (
                np.power(pars_under.pre_factor, 1.0 - pressure_log) * 
                np.power(pars_over.pre_factor, pressure_log)
                )
            temp_exponent = (
                pars_under.temp_exponent * (1.0 - pressure_log) + 
                pars_over.temp_exponent * pressure_log
                )
            activation_energy = (
                pars_under.act_energy * (1.0 - pressure_log) + 
                pars_over.act_energy * pressure_log
                )
            parameters = [Arrhenius(pre_exponential, temp_exponent, activation_energy)]
            
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
    reactions : list of ReacInfo
        List of one or more converted reaction(s)

    """
    assert rxn.plog
    assert method in ['interpolate', 'nearest']

    if method == 'nearest':
        rate_parameters = find_nearest_parameters(pressure, rxn.plog_par)
    elif method == 'interpolate':
        rate_parameters = interpolate_parameters(pressure, rxn.plog_par)

    # return one or more reactions
    reactions = []
    for rate in rate_parameters:
        reactions.append(ReacInfo(
            rxn.rev, rxn.reac, rxn.reac_nu, rxn.prod, rxn.prod_nu, 
            rate.pre_factor, rate.temp_exponent, rate.act_energy
            ))

    return reactions


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

    new_reacs = []
    for rxn in reacs:
        if rxn.plog:
            new_rxns = convert_plog(rxn, pressure, method)
            for r in new_rxns:
                new_reacs.append(r)
        else:
            new_reacs.append(rxn)

    # write new reaction list to new file
    write_mech(output_file, elems, specs, new_reacs)
