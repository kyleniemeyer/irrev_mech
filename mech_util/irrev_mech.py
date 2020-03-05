"""Makes reactions in chemical kinetic model all irreversible.
"""

# Standard libraries
import copy
import warnings
import logging
from multiprocessing import Pool
from itertools import repeat

try:
    import numpy as np
except ImportError:
    print('Error: NumPy must be installed.')
    raise
try:
    from scipy.optimize import leastsq
except ImportError:
    print('Error: SciPy must be installed.')
    raise

# Local imports
from .chem_utilities import units, Q_, GAS_CONSTANT, calc_spec_smh
from .mech_interpret import read_mech
from .write_mech import write_mech


def calc_rate_coeff(p, T):
    """Calculate Arrhenius reaction rate coefficient."""
    A, b, E = p
    k = A * np.exp(b * np.log(T) - (E / T))
    return k


def calc_rev_rate_coeff(T, p_Arr, specs, rxn):
    """Calculate reverse Arrhenius rate coefficient."""

    # calculate forward reaction rate (ignoring pressure dependence, since it
    # is accounted for in both directions by the specific formulation)
    k_fwd = calc_rate_coeff(p_Arr, T)

    # equilibrium constant
    # first get entropy minus enthalpy for all species
    spec_smh = calc_spec_smh(T, specs)

    Kp = 0.0
    # products
    for sp in rxn.prod:
        isp = next(i for i in range(len(specs)) if specs[i].name == sp)
        Kp += rxn.prod_nu[rxn.prod.index(sp)] * spec_smh[isp]

    # reactants
    for sp in rxn.reac:
        isp = next(i for i in range(len(specs)) if specs[i].name == sp)
        Kp -= rxn.reac_nu[rxn.reac.index(sp)] * spec_smh[isp]

    Kc = np.exp(Kp) * (
        Q_(1.0, 'atm').to('pascal').magnitude / (GAS_CONSTANT.magnitude * T)
        ) ** (sum(rxn.prod_nu) - sum(rxn.reac_nu))
    k_rev = k_fwd / Kc

    return k_rev


def residuals(p, y, x):
    """Residual for calculating rate coefficient."""
    err = y - calc_rate_coeff(p, x)
    return err


def calc_rev_Arrhenius(specs, rxn, rxn_id, Tfit, coeffs):
    """Calculate reverse Arrhenius coefficients for a particular reaction.

    Using three temperatures, fit reverse Arrhenius coefficients by
    calcluating forward and reverse reaction rates.

    Input
    rxn_id: reaction index
    rxn:  reaction object
    Tfit: tuple of three temperatures
    """

    # various constants for ease of calculation
    T1 = Tfit[0]
    T2 = Tfit[1]
    T3 = Tfit[2]

    A = coeffs[0]
    b = coeffs[1]
    E = coeffs[2]

    x1 = np.log(T1)
    x2 = np.log(T2)
    x3 = np.log(T3)

    p_Arr = A, b, E.to(units.kelvin).magnitude

    # calculate reverse reaction rates for each temperature

    # T1
    kr1 = calc_rev_rate_coeff(T1, p_Arr, specs, rxn)

    # T2
    kr2 = calc_rev_rate_coeff(T2, p_Arr, specs, rxn)

    # T3
    kr3 = calc_rev_rate_coeff(T3, p_Arr, specs, rxn)

    a1 = np.log(kr1)
    a2 = np.log(kr2)
    a3 = np.log(kr3)

    den = x1 * T1 * (T3 - T2) + x2 * T2 * (T1 - T3) + x3 * T3 * (T2 - T1)
    br = (a1 * T1 * (T3 - T2) + a2 * T2 * (T1 - T3) +
          a3 * T3 * (T2 - T1)
          ) / den
    Er = T1 * T2 * T3 * (a1 * (x2 - x3) + a2 * (x3 - x1) +
                         a3 * (x1 - x2)
                         ) / den
    Ar = (a1 * T1 * (x2 * T2 - x3 * T3) + a2 * T2 * (x3 * T3 - x1 * T1) +
          a3 * T3 * (x1 * T1 - x2 * T2)
          ) / den
    Ar = np.exp(Ar)

    # Now perform nonlinear least-squares minimization using these
    # values as the initial guesses.
    x = np.linspace(T1, T3, 1000)
    y = np.zeros(len(x))
    for idx, val in enumerate(x):
        y[idx] = calc_rev_rate_coeff(val, p_Arr, specs, rxn)

    # Start with low number of max function evals, increase if needed.
    warnings.filterwarnings('error')
    for mx in [1000, 5000, 10000, 20000, 40000]:
        try:
            val_lsq = leastsq(residuals, [Ar, br, Er], args=(y, x), maxfev=mx)
            break
        except RuntimeWarning:
            continue
    else:
        logging.info('Warning: minimization failed to converge for reaction ' +
              str(rxn_id) + '.'
              )

    parameters = val_lsq[0]
    return [parameters[0], parameters[1], parameters[2] * units.kelvin]


def process_reaction(arg):
    """
    Worker process for the multiprocessing. The single argument is expected to be a tuple
    with elements:
        0) Instance of ReacInfo
        1) List of species in the mechanism
        2) Index of the current reaction
        3) Tuple of fitting temperatures
    """
    rxn = arg[0]
    specs = arg[1]
    idx = arg[2]
    Tfit = arg[3]
    if not rxn.rev:
        return (idx, rxn, None)
    # create 2 irreversible reactions from reversible
    rxn.rev = False
    irrev_rxn = copy.deepcopy(rxn)

    # switch reactants and products
    irrev_rxn.reac = copy.copy(rxn.prod)
    irrev_rxn.reac_nu = copy.copy(rxn.prod_nu)
    irrev_rxn.prod = copy.copy(rxn.reac)
    irrev_rxn.prod_nu = copy.copy(rxn.reac_nu)

    # Calculate explicit reverse Arrhenius parameters for reaction
    if not rxn.rev_par:
        coeffs = [rxn.A, rxn.b, rxn.E]
        rev_par = calc_rev_Arrhenius(specs, rxn, idx, Tfit, coeffs)
    else:
        rev_par = rxn.rev_par

    irrev_rxn.A = rev_par[0]
    irrev_rxn.b = rev_par[1]
    irrev_rxn.E = rev_par[2]

    if rxn.pdep:
        # get reverse high-/low-pressure limit coeffs
        if rxn.low:
            coeffs = rxn.low
        elif rxn.high:
            coeffs = rxn.high
        
        rev_par = calc_rev_Arrhenius(specs, rxn, idx, Tfit, coeffs)

        if rxn.low:
            irrev_rxn.low = copy.copy(rev_par)
        elif rxn.high:
            irrev_rxn.high = copy.copy(rev_par)

    elif rxn.plog:
        # Pressure-log reaction
        irrev_rxn.plog = True
        irrev_rxn.plog_par = []
        for par in rxn.plog_par:
            rev_par = calc_rev_Arrhenius(
                specs, rxn, idx, Tfit, [par[1], par[2], par[3]]
                )
            plog_par = [par[0], rev_par[0], rev_par[1], rev_par[2]]
            irrev_rxn.plog_par.append(plog_par)

    elif rxn.cheb:
        irrev_rxn.cheb = True
        raise NotImplementedError('CHEB reactions not yet supported')

    rxn.rev_par = []
    irrev_rxn.rev_par = []
    return (idx, rxn, irrev_rxn)


def convert_mech_irrev(
    mech_name, therm_name=None, temp_range=[300., 5000.], 
    output_file="mech_irrev.txt", n_procs=None
    ):
    """Convert Chemkin-style mechanism with reversible reactions.

    Input
    mech_name: string with reaction mechanism filename (e.g. 'mech.dat')
    therm_name: string with thermodynamic database filename (e.g. 'therm.dat')
                or None if info in mech_name
    """

    # interpret reaction mechanism file
    [elems, specs, reacs] = read_mech(mech_name, therm_name)

    # tuple holding fit temperatures
    #Tfit = 300.0, 1000.0, 5000.0
    Tmid = temp_range[0] + 0.5*(temp_range[1] - temp_range[0])
    Tfit = temp_range[0], Tmid, temp_range[1]

    # Check for any Chebyshev reactions; not currently supported
    if any([rxn for rxn in reacs if rxn.cheb]):
        raise NotImplementedError('CHEB reactions not yet supported')

    if n_procs and n_procs > 1:
        with Pool(processes=n_procs) as pool:
            result = pool.map(
                process_reaction, 
                zip(reacs, repeat(specs), 
                    [reacs.index(rxn) for rxn in reacs], repeat(Tfit)
                    )
                )
    else:
        result = []
        for idx, rxn in enumerate(reacs):
            input = tuple([rxn, specs, idx, Tfit])
            result.append(process_reaction(input))


    reacs = []
    for idx, rxn, irrev_rxn in sorted(result, key=lambda tup: tup[0]):
        # now recreate reaction list
        reacs.append(rxn)
        if irrev_rxn:
            reacs.append(irrev_rxn)

    # Need to reevaluate duplicate reactions. Some marked as duplicate when
    # reversible may no longer be.
    dup_reacs = [rxn for rxn in reacs if rxn.dup]
    for rxn in dup_reacs:
        fnd_dup = False

        for rxn2 in dup_reacs:
            if rxn2 == rxn: continue

            # Compare lists of reactants and products;
            # shouldn't need to also compare stoich coefficients
            if (sorted(rxn.reac) == sorted(rxn2.reac) and
                sorted(rxn.prod) == sorted(rxn2.prod)):
                fnd_dup = True
                break

        # If no duplicates,
        if not fnd_dup:
            reacs[reacs.index(rxn)].dup = False

    # write new reaction list to new file
    write_mech(output_file, elems, specs, reacs)

    return
