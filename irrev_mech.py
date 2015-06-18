#! /usr/bin/env python

# Python 2 compatibility
from __future__ import division
from __future__ import print_function

# Standard libraries
import copy
import math
from argparse import ArgumentParser
import numpy as np
from scipy.optimize import leastsq
import warnings

# Local imports
import chem_utilities as chem
import mech_interpret as mech


def calc_rate_coeff(p, T):
    """Calculate Arrhenius reaction rate coefficient."""
    A, b, E = p
    k = A * np.exp(b * np.log(T) - (E / T))
    return k


def calc_rev_rate_coeff(T, p_Arr, specs, reac):
    """Calculate reverse Arrhenius rate coefficient."""

    # calculate forward reaction rate (ignoring pressure dependence, since it
    # is accounted for in both directions by the specific formulation)
    k_fwd = calc_rate_coeff(p_Arr, T)

    # equilibrium constant
    # first get entropy minus enthalpy for all species
    spec_smh = chem.calc_spec_smh(T, specs)

    Kp = 0.0
    # products
    for sp in reac.prod:
        isp = next(i for i in xrange(len(specs)) if specs[i].name == sp)
        Kp += reac.prod_nu[reac.prod.index(sp)] * spec_smh[isp]

    # reactants
    for sp in reac.reac:
        isp = next(i for i in xrange(len(specs)) if specs[i].name == sp)
        Kp -= reac.reac_nu[reac.reac.index(sp)] * spec_smh[isp]

    Kp = math.exp(Kp)
    Kc = Kp * (chem.PA / (chem.RU * T)) ** (sum(reac.prod_nu) - sum(reac.reac_nu))
    k_rev = k_fwd / Kc

    return k_rev


def residuals(p, y, x):
    """Residual for calculating rate coefficient."""
    A, b, E = p
    err = y - calc_rate_coeff(p, x)
    return err


def calc_rev_Arrhenius(specs, reac, reac_id, Tfit, units, coeffs):
    """Calculate reverse Arrhenius coefficients for a particular reaction.

    Using three temperatures, fit reverse Arrhenius coefficients by
    calcluating forward and reverse reaction rates.

    Input
    reac_id: reaction index
    reac:  reaction object
    Tfit: tuple of three temperatures
    """

    # various constants for ease of calculation
    T1 = Tfit[0]
    T2 = Tfit[1]
    T3 = Tfit[2]

    A = coeffs[0]
    b = coeffs[1]

    x1 = math.log(T1)
    x2 = math.log(T2)
    x3 = math.log(T3)

    # use correct gas constant based on units
    efac = 1.0
    if 'kelvin' not in units:
        if 'kcal/mole' in units:
            efac = 4184.0 / chem.RU_JOUL
        elif 'cal/mole' in units:
            efac = 4.184 / chem.RU_JOUL
        elif 'kjoule' in units:
            efac = 1000.0 / chem.RU_JOUL
        elif 'joules' in units:
            efac = 1.00 / chem.RU_JOUL
        elif 'evolt' in units:
            efac = 11595.0
        else:
            # default is cal/mole
            efac = 4.184 / chem.RU_JOUL
    E = coeffs[2] * efac

    p_Arr = A, b, E

    # calculate reverse reaction rates for each temperature

    # T1
    kr1 = calc_rev_rate_coeff(T1, p_Arr, specs, reac)

    # T2
    kr2 = calc_rev_rate_coeff(T2, p_Arr, specs, reac)

    # T3
    kr3 = calc_rev_rate_coeff(T3, p_Arr, specs, reac)

    a1 = math.log(kr1)
    a2 = math.log(kr2)
    a3 = math.log(kr3)

    den = x1 * T1 * (T3 - T2) + x2 * T2 * (T1 - T3) + x3 * T3 * (T2 - T1)
    br = (a1 * T1 * (T3 - T2) + a2 * T2 * (T1 - T3) + a3 * T3 * (T2 - T1)) / den
    Er = T1 * T2 * T3 * (a1 * (x2 - x3) + a2 * (x3 - x1) + a3 * (x1 - x2)) / den
    Ar = (a1 * T1 * (x2 * T2 - x3 * T3) + a2 * T2 * (x3 * T3 - x1 * T1) + a3 * T3 * (x1 * T1 - x2 * T2)) / den
    Ar = math.exp(Ar)

    # Now perform nonlinear least-squares minimization using these
    # values as the initial guesses.
    x = np.linspace(T1, T3, 1000)
    y = np.zeros(len(x))
    for idx, val in enumerate(x):
        y[idx] = calc_rev_rate_coeff(val, p_Arr, specs, reac)

    # Start with low number of max function evals, increase if needed.
    warnings.filterwarnings('error')
    for mx in [1000, 5000, 10000, 20000]:
        try:
            val_lsq = leastsq(residuals, [Ar, br, Er], args=(y, x), maxfev=mx)
            break
        except RuntimeWarning:
            continue
    else:
        print('Warning: minimization failed to converge for reaction ' +
              reac_id + '.'
              )

    # correct gas constant based on units
    Er /= efac

    return [Ar, br, Er]


def write_mech(filename, elems, specs, reacs, units):
    """Write Chemkin-format mechanism.

    Input
    """
    file = open(filename, 'w')

    # elements
    file.write('elements\n')

    for e in elems:
        # write atomic weight if necessary
        if e in mech.elem_new:
            file.write(e + ' /' + str(mech.elem_wt[e.lower()]) + '/ \n')
        else:
            file.write(e + '\n')

    file.write('end\n\n')

    # species
    file.write('species\n')

    for sp in specs:
        file.write(sp.name + '\n')

    file.write('end\n\n')

    # reactions
    file.write('reactions                           ' + units + '\n')

    for rxn in reacs:
        line = ''

        # reactants
        for sp in rxn.reac:
            isp = rxn.reac.index(sp)
            # print stoich coefficient if other than one
            if rxn.reac_nu[isp] != 1:
                line += str(rxn.reac_nu[isp]) + sp
            else:
                line += sp

            if (len(rxn.reac) - 1) > isp:
                line += '+'

        # third body in reactants
        if rxn.pdep:
            if rxn.thd:
                line += '(+m)'
            else:
                line += '(+{:s})'.format(rxn.pdep_sp)
        elif rxn.thd:
            line += '+m'

        if rxn.rev:
            line += '='
        else:
            line += '=>'

        # products
        for sp in rxn.prod:
            isp = rxn.prod.index(sp)
            # print stoich coefficient if other than one
            if rxn.prod_nu[isp] != 1:
                line += str(rxn.prod_nu[isp]) + sp
            else:
                line += sp

            if (len(rxn.prod) - 1) > isp:
                line += '+'

        # third body in products
        if rxn.pdep:
            if rxn.thd:
                line += '(+m)'
            else:
                line += '(+{:s})'.format(rxn.pdep_sp)
        elif rxn.thd:
            line += '+m'

        # now add Arrhenius coefficients to the same line
        line += ' {:.4e} {:.4e} {:.4e}'.format(rxn.A, rxn.b, rxn.E)

        line += '\n'
        file.write(line)


        # line for reverse Arrhenius parameters, if any
        if rxn.rev:
            line = '  rev/ {:.4e}  {:.4e}  {:.4e} /\n'.format(rxn.rev_par[0], rxn.rev_par[1], rxn.rev_par[2])
            file.write(line)


        # write Lindemann low- or high-pressure limit Arrhenius parameters
        if rxn.pdep:
            if rxn.low:
                line = '  low /{:.4e}  {:.4e}  {:.4e} /\n'.format(rxn.low[0], rxn.low[1], rxn.low[2])
            else:
                line = '  high /{:.4e}  {:.4e}  {:.4e} /\n'.format(rxn.high[0], rxn.high[1], rxn.high[2])

            file.write(line)

        # write Troe parameters if any
        if rxn.troe:
            troe = rxn.troe_par
            if len(troe) == 3:
                line = '  troe/ {:.4e} {:.4e} {:.4e} /\n'.format(troe[0], troe[1], troe[2])
            else:
                line = '  troe/ {:.4e} {:.4e} {:.4e} {:.4e} /\n'.format(troe[0], troe[1], troe[2], troe[3])
            file.write(line)


        # write SRI parameters if any
        if rxn.sri:
            sri = rxn.sri_par
            if len(sri) == 3:
                line = '  sri/ {:.4e} {:.4e} {:.4e} /\n'.format(sri[0], sri[1], sri[2])
            else:
                line = '  sri/ {:.4e} {:.4e} {:.4e} {:.4e} {:.4e} /\n'.format(sri[0], sri[1], sri[2], sri[3], sri[4])


        # third-body efficiencies
        if rxn.thd_body:
            line = '  '
            for thd_body in rxn.thd_body:
                thd_eff = '{:.2f}'.format(thd_body[1])
                line += thd_body[0] + '/' + thd_eff + '/ '

                # move to next line if long
                if len(line) >= 60 and rxn.thd_body.index(thd_body) is not (len(rxn.thd_body) - 1):
                    line += '\n'
                    file.write(line)
                    line = '  '

            line += '\n'
            file.write(line)

        # duplicate reaction flag
        if rxn.dup:
            file.write('  DUPLICATE\n')

    file.write('end')

    file.close()
    return


def convert_mech_irrev(mech_name, therm_name, temp_range):
    """Convert Chemkin-style mechanism with reversible reactions.

    Input
    mech_name: string with reaction mechanism filename (e.g. 'mech.dat')
    therm_name: string with thermodynamic database filename (e.g. 'therm.dat') or None if info in mech_name
    """

    # interpret reaction mechanism file
    [elems, specs, reacs, units] = mech.read_mech(mech_name)

    # interpret thermodynamic database file (if it exists)
    if therm_name:
        file = open(therm_name, 'r')
        mech.read_thermo(file, elems, specs)
        file.close()
    else:
        # copy therm data into new file
        mech_file = open(mech_name, 'r')
        file = open('therm_irrev.txt', 'w')

        flag = False
        for line in mech_file:
            if line[0:4].lower() == 'ther':
                file.write('thermo\n')
                flag = True
                continue

            if flag: file.write(line)

            if flag and line[0:3].lower() == 'end': break

        file.close()

    # tuple holding fit temperatures
    #Tfit = 300.0, 1000.0, 5000.0
    Tmid = temp_range[0] + 0.5*(temp_range[1] - temp_range[0])
    Tfit = temp_range[0], Tmid, temp_range[1]

    # now loop through reversible reactions
    for rxn in [rxn for rxn in reacs[:] if rxn.rev]:

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
            rev_par = calc_rev_Arrhenius(specs, rxn, reacs.index(rxn), Tfit, units, coeffs)
        else:
            rev_par = rxn.rev_par

        irrev_rxn.A = rev_par[0]
        irrev_rxn.b = rev_par[1]
        irrev_rxn.E = rev_par[2]

        # get reverse high-/low-pressure limit coeffs
        if rxn.pdep:

            if rxn.low:
                coeffs = rxn.low
            elif rxn.high:
                coeffs = rxn.high

            rev_par = calc_rev_Arrhenius(specs, rxn, reacs.index(rxn), Tfit, units, coeffs)

            if rxn.low:
                irrev_rxn.low = copy.copy(rev_par)
            elif rxn.high:
                irrev_rxn.high = copy.copy(rev_par)

        rxn.rev_par = []
        irrev_rxn.rev_par = []

        # now insert into reaction list
        reacs.insert(reacs.index(rxn) + 1, irrev_rxn)

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

    mod_mech = 'mech_irrev.txt'

    # write new reaction list to new file
    write_mech(mod_mech, elems, specs, reacs, units)

    return


if __name__ == "__main__":

    # command line arguments
    parser = ArgumentParser(description = 'Generates reaction mechanism with '
                                          'only irreversible reactions.')
    parser.add_argument('-m', '--mech',
                        type = str,
                        required = True,
                        help = 'Input mechanism filename (e.g., mech.dat).')
    parser.add_argument('-t', '--thermo',
                        type = str,
                        default = None,
                        help = 'Thermodynamic database filename (e.g., '
                        'therm.dat), or nothing if in mechanism.')
    parser.add_argument('-r', '--range',
                        type = float, nargs=2,
                        default = [300.0, 5000.0],
                        help = 'Temperature range for fit in Kelvin (e.g., 300 5000).')

    args = parser.parse_args()
    convert_mech_irrev(args.mech, args.thermo, args.range)
