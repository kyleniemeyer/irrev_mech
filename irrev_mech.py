#! /usr/bin/env python

import math
from chem_utilities import *
from mech_interpret import *

def calc_rev_Arrhenius(specs, reac, Tfit, units, coeffs):
    """Calculate reverse Arrhenius coefficients.
    
    Using three temperatures, fit reverse Arrhenius coefficients by
    calcluating forward and reverse reaction rates.
    
    Input
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
            #E = reac.E * 1000.0 / RUC
            efac = 4184.0 / RU_JOUL
        elif 'cal/mole' in units:
            #E = reac.E / RUC
            efac = 4.184 / RU_JOUL
        elif 'kjoule' in units:
            efac = 1000.0 / RU_JOUL
        elif 'joules' in units:
            # Ru = 8.3144621 J/(mol * K)
            #E = reac.E / 8.3144621
            efac = 1.00 / RU_JOUL
        elif 'evolt' in units:
            # Ru = 5.189e19 eV/(mol * K)
            efac = 11595.
        else:
            # default is cal/mole
            efac = 4.184 / RU_JOUL
    E = coeffs[2] * efac
    
    
    # calculate reverse reaction rates for each temperature
    
    # T1
    
    # calculate forward reaction rate (ignoring pressure dependence, since it
    # is accounted for in both directions by the specific formulation
    k1 = A * math.exp(b * x1 - (E / T1))
    
    # equilibrium constant
    # first get entropy minus enthalpy for all species
    spec_smh = calc_spec_smh(T1, specs)
    
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
    Kc = Kp * (PA / (RU * T1)) ** (sum(reac.prod_nu) - sum(reac.reac_nu))
    kr1 = k1 / Kc
    
    # T2
    
    # calculate forward reaction rate (ignoring pressure dependence, since it
    # is accounted for in both directions by the specific formulation
    k2 = A * math.exp(b * x2 - (E / T2))
    
    # equilibrium constant
    # first get entropy minus enthalpy for all species
    spec_smh = calc_spec_smh(T2, specs)
    
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
    Kc = Kp * (PA / (RU * T2)) ** (sum(reac.prod_nu) - sum(reac.reac_nu))
    kr2 = k2 / Kc
    
    # T3
    
    # calculate forward reaction rate (ignoring pressure dependence, since it
    # is accounted for in both directions by the specific formulation
    k3 = A * math.exp(b * x3 - (E / T3))
    
    # equilibrium constant
    # first get entropy minus enthalpy for all species
    spec_smh = calc_spec_smh(T3, specs)
    
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
    Kc = Kp * (PA / (RU * T3)) ** (sum(reac.prod_nu) - sum(reac.reac_nu))
    kr3 = k3 / Kc
    
    a1 = math.log(kr1)
    a2 = math.log(kr2)
    a3 = math.log(kr3)
    
    den = x1 * T1 * (T3 - T2) + x2 * T2 * (T1 - T3) + x3 * T3 * (T2 - T1)
    br = (a1 * T1 * (T3 - T2) + a2 * T2 * (T1 - T3) + a3 * T3 * (T2 - T1)) / den
    Er = T1 * T2 * T3 * (a1 * (x2 - x3) + a2 * (x3 - x1) + a3 * (x1 - x2)) / den
    Ar = (a1 * T1 * (x2 * T2 - x3 * T3) + a2 * T2 * (x3 * T3 - x1 * T1) + a3 * T3 * (x1 * T1 - x2 * T2)) / den
    Ar = math.exp(Ar)
    
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
        if e in elem_mw_new:
            file.write(e + ' /' + str(elem_mw_new[e.lower()]) + '/ \n')
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
                line += str(rxn.reac_nu[isp]) + ' ' + sp
            else:
                line += sp
            
            if (len(rxn.reac) - 1) > isp:
                line += ' + '
        
        # third body in reactants
        if rxn.thd:
            line += ' + m'
        elif rxn.pdep:
            line += ' (+ {:s})'.format(rxn.pdep_sp)
        
        if rxn.rev:
            line += ' = '
        else:
            line += ' => '
        
        # products
        for sp in rxn.prod:
            isp = rxn.prod.index(sp)
            # print stoich coefficient if other than one
            if rxn.prod_nu[isp] != 1:
                line += str(rxn.prod_nu[isp]) + ' ' + sp
            else:
                line += sp
            
            if (len(rxn.prod) - 1) > isp:
                line += ' + '
        
        # third body in products
        if rxn.thd:
            line += ' + m'
        elif rxn.pdep:
            line += ' (+ {:s})'.format(rxn.pdep_sp)
        
        # now add Arrhenius coefficients to the same line
        line += '    {:.4e}  {:.4e}  {:.4e}'.format(rxn.A, rxn.b, rxn.E)
        
        line += '\n'
        file.write(line)
        
        
        # line for reverse Arrhenius parameters, if any
        if rxn.rev:
            line = '    rev /  {:.4e}  {:.4e}  {:.4e}  / \n'.format(rxn.rev_par[0], rxn.rev_par[1], rxn.rev_par[2])
            file.write(line)
        
        
        # write Lindemann low- or high-pressure limit Arrhenius parameters
        if rxn.pdep:
            if rxn.low:
                line = '    low / {:.4e}  {:.4e}  {:.4e} / \n'.format(rxn.low[0], rxn.low[1], rxn.low[2])
            else:
                line = '    high / {:.4e}  {:.4e}  {:.4e} / \n'.format(rxn.high[0], rxn.high[1], rxn.high[2])
            
            file.write(line)
        
        # write Troe parameters if any
        if rxn.troe:
            troe = rxn.troe_par
            if len(troe) == 3:
                line = '    troe / {:.4e} {:.4e} {:.4e} / \n'.format(troe[0], troe[1], troe[2])
            else:
                line = '    troe / {:.4e} {:.4e} {:.4e} {:.4e} / \n'.format(troe[0], troe[1], troe[2], troe[3])
            file.write(line)
        
        
        # write SRI parameters if any
        if rxn.sri:
            sri = rxn.sri_par
            if len(sri) == 3:
                line = '    sri / {:.4e} {:.4e} {:.4e} / \n'.format(sri[0], sri[1], sri[2])
            else:
                line = '    sri / {:.4e} {:.4e} {:.4e} {:.4e} {:.4e} / \n'.format(sri[0], sri[1], sri[2], sri[3], sri[4])
        
        
        # third-body efficiencies
        if rxn.thd_body:
            line = '    '
            for thd_body in rxn.thd_body:
                thd_eff = '{:.2f}'.format(thd_body[1])
                line += thd_body[0] + ' / ' + thd_eff + ' / '
                
                # move to next line if long
                if len(line) >= 60:
                    line += '\n'
                    file.write(line)
                    line = '    '
            
            line += '\n'
            file.write(line)
        
        # duplicate reaction flag
        if rxn.dup:
            file.write('    DUPLICATE\n')
    
    file.write('end')
    
    file.close()
    return


def convert_mech_irrev(mech_name, therm_name = None):
    """Convert Chemkin-style mechanism with reversible reactions.
    
    Input
    mech_name: string with reaction mechanism filename (e.g. 'mech.dat')
    therm_name: string with thermodynamic database filename (e.g. 'therm.dat') or nothing if info in mech_name
    """
    import copy
    
    elems = []
    specs = []
    reacs = []
    
    # interpret reaction mechanism file
    [num_e, num_s, num_r, units] = read_mech(mech_name, elems, specs, reacs)
    
    # interpret thermodynamic database file (if it exists)
    if therm_name:
        file = open(therm_name, 'r')
        read_thermo(file, elems, specs)
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
    Tfit = 1000.0, 1750.0, 2500.0
    
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
        
        # reaction doesn't have explicit reverse Arrhenius parameters, so calculate
        if not rxn.rev_par:
            coeffs = [rxn.A, rxn.b, rxn.E]
            rev_par = calc_rev_Arrhenius(specs, rxn, Tfit, units, coeffs)
        
        irrev_rxn.A = rev_par[0]
        irrev_rxn.b = rev_par[1]
        irrev_rxn.E = rev_par[2]
        
        # get reverse high-/low-pressure limit coeffs
        if rxn.pdep:
            
            if rxn.low:
                coeffs = rxn.low
            elif rxn.high:
                coeffs = rxn.high
            
            rev_par = calc_rev_Arrhenius(specs, rxn, Tfit, units, coeffs)
            
            if rxn.low:
                irrev_rxn.low = copy.copy(rev_par)
            elif rxn.high:
                irrev_rxn.high = copy.copy(rev_par)
        
        rxn.rev_par = []
        irrev_rxn.rev_par = []
        
        # now insert into reaction list
        reacs.insert(reacs.index(rxn) + 1, irrev_rxn)
    
    mod_mech = 'mech_irrev.txt'
    
    # write new reaction list to new file
    write_mech(mod_mech, elems, specs, reacs, units)
    
    return


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) == 2:
        convert_mech_irrev(sys.argv[1])
    elif len(sys.argv) == 3:
        convert_mech_irrev(sys.argv[1], sys.argv[2])
    else:
        print 'Incorrect number of arguments'
