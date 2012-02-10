#! /usr/bin/env python

import math

# universal gas constants, cgs units
RU = 8.314510e7 # erg/(mole * K)
RUC = 1.9858775 # cal/(mole * K)
RU_JOUL = 8.314510e0
# pressure of one standard atmosphere, dynes/cm^2
PA = 1.01325e6

# element atomic weights
elem_mw = dict( [ 
    ('h',   1.00797),  ('he',   4.00260), ('li',   6.93900), ('be',   9.01220), 
    ('b',  10.81100),  ('c',   12.01115), ('n',   14.00670), ('o',   15.99940), 
    ('f',  18.99840),  ('ne',  20.18300), ('na',  22.98980), ('mg',  24.31200),
    ('al', 26.98150),  ('si',  28.08600), ('p',   30.97380), ('s',   32.06400), 
    ('cl', 35.45300),  ('ar',  39.94800), ('k',   39.10200), ('ca',  40.08000), 
    ('sc', 44.95600),  ('ti',  47.90000), ('v',   50.94200), ('cr',  51.99600), 
    ('mn', 54.93800),  ('fe',  55.84700), ('co',  58.93320), ('ni',  58.71000), 
    ('cu', 63.54000),  ('zn',  65.37000), ('ga',  69.72000), ('ge',  72.59000), 
    ('as', 74.92160),  ('se',  78.96000), ('br',  79.90090), ('kr',  83.80000), 
    ('rb', 85.47000),  ('sr',  87.62000), ('y',   88.90500), ('zr',  91.22000),
    ('nb', 92.90600),  ('mo',  95.94000), ('tc',  99.00000), ('ru', 101.07000),
    ('rh', 102.90500), ('pd', 106.40000), ('ag', 107.87000), ('cd', 112.40000),
    ('in', 114.82000), ('sn', 118.69000), ('sb', 121.75000), ('te', 127.60000),
    ('i',  126.90440), ('xe', 131.30000), ('cs', 132.90500), ('ba', 137.34000),
    ('la', 138.91000), ('ce', 140.12000), ('pr', 140.90700), ('nd', 144.24000),
    ('pm', 145.00000), ('sm', 150.35000), ('eu', 151.96000), ('gd', 157.25000),
    ('tb', 158.92400), ('dy', 162.50000), ('ho', 164.93000), ('er', 167.26000),
    ('tm', 168.93400), ('yb', 173.04000), ('lu', 174.99700), ('hf', 178.49000),
    ('ta', 180.94800), ('w',  183.85000), ('re', 186.20000), ('os', 190.20000),
    ('ir', 192.20000), ('pt', 195.09000), ('au', 196.96700), ('hg', 200.59000),
    ('tl', 204.37000), ('pb', 207.19000), ('bi', 208.98000), ('po', 210.00000),
    ('at', 210.00000), ('rn', 222.00000), ('fr', 223.00000), ('ra', 226.00000),
    ('ac', 227.00000), ('th', 232.03800), ('pa', 231.00000), ('u',  238.03000),
    ('np', 237.00000), ('pu', 242.00000), ('am', 243.00000), ('cm', 247.00000),
    ('bk', 249.00000), ('cf', 251.00000), ('es', 254.00000), ('fm', 253.00000),
    ('d',    2.01410), ('e', 5.48578e-4) ] )

# dict for any new element definitions
elem_mw_new = dict()

def read_str_num(string):
    """Pull list of floats from a string"""
        
    # separate string into space-delimited strings of numbers
    num_str = string.split()
    
    nums = []
    for n in num_str:
        nums.append( float(n) )
    
    return nums


def split_str(seq, length):
    """Separate a string seq into length-sized pieces"""
    return [seq[i : i + length] for i in range(0, len(seq), length)]


class ReacInfo:
    """Class for reaction information"""
    
    def __init__(self, rev, reactants, reac_nu, products, prod_nu, A, b, E):
        self.reac = reactants
        self.reac_nu = reac_nu
        self.prod = products
        self.prod_nu = prod_nu
        
        # Arrhenius coefficients
        self.A = A
        self.b = b
        self.E = E
        
        # reversible reaction properties
        self.rev = rev
        self.rev_par = [] # reverse A, b, E
        
        # duplicate reaction
        self.dup = False
        
        # third-body efficiencies
        self.thd = False
        self.thd_body = [] # in pairs with species and efficiency
        
        # pressure dependence
        self.pdep = False
        self.pdep_sp = ''
        self.low = []
        self.high = []
        
        self.troe = False
        self.troe_par = []
        
        self.sri = False
        self.sri_par = []


class SpecInfo:
    """Class for species information"""
    
    def __init__(self, name):
        self.name = name
        
        # elemental composition
        self.elem = []
        # molecular weight
        self.mw = 0.0
        # high-temp range thermodynamic coefficients
        self.hi = [0.0 for j in range(7)]
        # low-temp range thermodynamic coefficients
        self.lo = [0.0 for j in range(7)]
        # temperature range for thermodynamic coefficients
        self.Trange = [0.0, 0.0, 0.0]


def read_mech(filename, elems, specs, reacs):
    """Read and interpret mechanism file for elements, species, and reactions.
    
    Doesn't support element names with digits.
    
    Input
    filename:  reaction mechanism filename (e.g. 'mech.dat')
    """
    
    file = open(filename, 'r')
    
    num_e = 0
    num_s = 0
    num_r = 0
    
    units = ''
    key = ''
    
    # start line reading loop
    while True:
        line = file.readline()
        
        # end of file
        if not line: break
        
        # skip blank or commented lines
        if line == '\n' or line == '\r\n' or line[0:1] == '!': continue
        
        # don't convert everything, since thermo needs to match (for Chemkin)
        ## convert to lowercase
        #line = line.lower()
        
        # remove any comments from end of line
        ind = line.find('!')
        if ind > 0: line = line[0:ind]
        
        # now determine key
        if line[0:4].lower() == 'elem':
            key = 'elem'
            
            # check for any entries on this line
            line_split = line.split()
            if len(line_split) > 1:
                ind = line.index( line_split[1] )
                line = line[ind:]
            else:
                continue
            
        elif line[0:4].lower() == 'spec':
            key = 'spec'
            
            # check for any entries on this line
            line_split = line.split()
            if len(line_split) > 1:
                ind = line.index( line_split[1] )
                line = line[ind:]
            else:
                continue
            
        elif line[0:4].lower() == 'reac':
            key = 'reac'
            
            # get Arrhenius coefficient units
            line_split = line.split()
            if len(line_split) > 1:
                units = line[ line.index(line_split[1]) : ].strip()
            else:
                # default units
                units = 'cal/mole'
            
            continue
            
        elif line[0:3].lower() == 'end':
            continue
        
        line = line.strip()
        
        if key == 'elem':
            # if any atomic weight declarations, replace / with spaces
            line = line.replace('/', ' ')
            
            line_split = line.split()
            e_last = ''
            for e in line_split:
                if e.isalpha():
                    if e[0:3] == 'end': continue
                    if e not in elems:
                        elems.append(e)
                        num_e += 1
                    e_last = e
                else:
                    # check either new element or updating existing atomic weight
                    if e_last in elem_mw:
                        elem_mw[e_last.lower()] = float(e)
                    
                    # in both cases add to 2nd dict to keep track
                    elem_mw_new[e_last.lower()] = float(e)
            
        elif key == 'spec':
            line_split = line.split()
            for s in line_split:
                if s[0:3] == 'end': continue
                if not next((sp for sp in specs if sp.name == s), None):
                    specs.append( SpecInfo(s) )
                    num_s += 1
            
        elif key == 'reac':
            # determine if reaction or auxiliary info line
            
            if '=' in line:
                # new reaction
                num_r += 1
                
                # get Arrhenius coefficients
                line_split = line.split()
                n = len(line_split)
                reac_A = float( line_split[n - 3] )
                reac_b = float( line_split[n - 2] )
                reac_E = float( line_split[n - 1] )
                
                ind = line.index( line_split[n - 3] )
                line = line[0:ind].strip()
                
                if '<=>' in line:
                    ind = line.index('<=>')
                    reac_rev = True
                    reac_str = line[0:ind].strip()
                    prod_str = line[ind + 3:].strip()
                elif '=>' in line:
                    ind = line.index('=>')
                    reac_rev = False
                    reac_str = line[0:ind].strip()
                    prod_str = line[ind + 2:].strip()
                else:
                    ind = line.index('=')
                    reac_rev = True
                    reac_str = line[0:ind].strip()
                    prod_str = line[ind + 1:].strip()
                
                thd = False
                pdep = False
                pdep_sp = ''
                
                reac_spec = []
                reac_nu = []
                prod_spec = []
                prod_nu = []
                
                # reactants
                
                # look for third-body species
                if '(' in reac_str: 
                    pdep = True
                    ind1 = reac_str.find('(')
                    ind2 = reac_str.find(')')
                    
                    # either 'm' or a specific species
                    pdep_sp = reac_str[ind1 + 1 : ind2].replace('+', ' ')
                    pdep_sp = pdep_sp.strip()
                    
                    # now remove from string
                    reac_str = reac_str[0:ind1] + reac_str[ind2 + 1:]
                
                reac_list = reac_str.split('+')
                
                for sp in reac_list:
                    
                    sp = sp.strip()
                    
                    # look for coefficient
                    if sp[0:1].isdigit(): 
                        # starts with number (coefficient)
                        
                        # search for first letter
                        for i in range( len(sp) ):
                            if sp[i : i + 1].isalpha(): break
                        
                        nu = sp[0:i]
                        if '.' in nu:
                            # float
                            nu = float(nu)
                        else:
                            # integer
                            nu = int(nu)
                        
                        sp = sp[i:].strip()
                    else:
                        # no coefficient given
                        nu = 1
                    
                    # check for third body
                    if sp == 'm' or sp == 'M':
                        thd = True
                        continue
                    
                    # check if species already in reaction
                    if sp not in reac_spec:
                        # new reactant
                        reac_spec.append(sp)
                        reac_nu.append(nu)
                    else:
                        # existing reactant
                        i = reac_spec.index(sp)
                        reac_nu[i] += nu
                
                # products
                
                # look for third-body species
                if '(' in prod_str: 
                    pdep = True
                    ind1 = prod_str.find('(')
                    ind2 = prod_str.find(')')
                    
                    # either 'm' or a specific species
                    pdep_sp = prod_str[ind1 + 1 : ind2].replace('+', ' ')
                    pdep_sp = pdep_sp.strip()
                    
                    # now remove from string
                    prod_str = prod_str[0:ind1] + prod_str[ind2 + 1:]
                
                prod_list = prod_str.split('+')
                
                for sp in prod_list:
                    
                    sp = sp.strip()
                    
                    # look for coefficient
                    if sp[0:1].isdigit(): 
                        # starts with number (coefficient)
                        
                        # search for first letter
                        for i in range( len(sp) ):
                            if sp[i : i + 1].isalpha(): break
                        
                        nu = sp[0:i]
                        if '.' in nu:
                            # float
                            nu = float(nu)
                        else:
                            # integer
                            nu = int(nu)
                        
                        sp = sp[i:].strip()
                    else:
                        # no coefficient given
                        nu = 1
                    
                    # check for third body
                    if sp == 'm' or sp == 'M':
                        thd = True
                        continue
                    
                    # check if species already in reaction
                    if sp not in prod_spec:
                        # new product
                        prod_spec.append(sp)
                        prod_nu.append(nu)
                    else:
                        # existing product
                        i = prod_spec.index(sp)
                        prod_nu[i] += nu
                
                # add reaction to list
                reacs.append( ReacInfo(reac_rev, reac_spec, reac_nu, prod_spec, prod_nu, reac_A, reac_b, reac_E) )
                reacs[num_r - 1].thd = thd
                reacs[num_r - 1].pdep = pdep
                if pdep: reacs[num_r - 1].pdep_sp = pdep_sp
                
            else:
                # auxiliary reaction info
                
                aux = line[0:3].lower()
                if aux == 'dup':
                    reacs[num_r - 1].dup = True
                    
                elif aux == 'rev':
                    line = line.replace('/', ' ')
                    line_split = line.split()
                    reacs[num_r - 1].rev_par.append( float( line_split[1] ) )
                    reacs[num_r - 1].rev_par.append( float( line_split[2] ) )
                    reacs[num_r - 1].rev_par.append( float( line_split[3] ) )
                    
                elif aux == 'low':
                    line = line.replace('/', ' ')
                    line_split = line.split()
                    reacs[num_r - 1].low.append( float( line_split[1] ) )
                    reacs[num_r - 1].low.append( float( line_split[2] ) )
                    reacs[num_r - 1].low.append( float( line_split[3] ) )
                    
                elif aux == 'hig':
                    line = line.replace('/', ' ')
                    line_split = line.split()
                    reacs[num_r - 1].high.append( float( line_split[1] ) )
                    reacs[num_r - 1].high.append( float( line_split[2] ) )
                    reacs[num_r - 1].high.append( float( line_split[3] ) )
                    
                elif aux == 'tro':
                    line = line.replace('/', ' ')
                    line_split = line.split()
                    reacs[num_r - 1].troe = True
                    reacs[num_r - 1].troe_par.append( float( line_split[1] ) )
                    reacs[num_r - 1].troe_par.append( float( line_split[2] ) )
                    reacs[num_r - 1].troe_par.append( float( line_split[3] ) )
                    
                    # optional fourth parameter
                    if len(line_split) > 4:
                        reacs[num_r - 1].troe_par.append( float( line_split[4] ) )
                    
                elif aux == 'sri':
                    line = line.replace('/', ' ')
                    line_split = line.split()
                    reacs[num_r - 1].sri = True
                    reacs[num_r - 1].sri_par.append( float( line_split[1] ) )
                    reacs[num_r - 1].sri_par.append( float( line_split[2] ) )
                    reacs[num_r - 1].sri_par.append( float( line_split[3] ) )
                    
                    # optional fourth and fifth parameters
                    if len(line_split) > 4:
                        reacs[num_r - 1].sri_par.append( float( line_split[4] ) )
                        reacs[num_r - 1].sri_par.append( float( line_split[5] ) )
                else:
                    # enhanced third body efficiencies
                    line = line.replace('/', ' ')
                    
                    line_split = line.split()
                    for i in range(0, len(line_split), 2):
                        reacs[num_r - 1].thd_body.append( [line_split[i], float(line_split[i + 1])] )
    
    return (num_e, num_s, num_r, units)


def read_thermo(filename, elems, specs):
    """Read and interpret thermodynamic database for species data.
    
    Reads the file therm.dat and returns the species thermodynamic coefficients
    as well as the species-specific temperature range values (if given)
    
    Input
    filename:  thermo database filename (e.g. 'therm.dat')
    elems: list of element names
    specs: list of species names (SpecInfo class)
    """
    
    # lines in thermo file are 80 characters long
    file = open(filename, 'r')
    
    # loop through intro lines
    while True:
        line = file.readline()
    
        # skip blank or commented lines
        if line == '\n' or line == '\r\n' or line[0:1] == '!': continue
    
        # skip 'thermo' at beginning
        if line[0:6].lower() == 'thermo': break
    
    # next line has common temperature ranges
    line = file.readline()
    
    T_ranges = read_str_num(line)
    
    # now start reading species thermo info
    while True:
        # first line of species info
        line = file.readline()
        
        # don't convert to lowercase, needs to match thermo for Chemkin
        #line = line.lower()
        
        # break if end of file
        if line is None: break
        if line[0:3] == 'end': break
        # skip blank/commented line
        if line == '\n' or line == '\r\n' or line[0:1] == '!': continue
        
        # species name, columns 0:18
        spec = line[0:18].strip()
        
        # apparently in some cases notes are in the columns of shorter species names
        # so make sure no spaces
        if spec.find(' ') > 0:
            spec = spec[0 : spec.find(' ')]
        
        # now need to determine if this species is in mechanism
        if next((sp for sp in specs if sp.name == spec), None):
            sp_ind = next(i for i in xrange(len(specs)) if specs[i].name == spec)
        else:
            # not in mechanism, read next three lines and continue
            line = file.readline()
            line = file.readline()
            line = file.readline()
            continue
        
        # set species to the one matched
        spec = specs[sp_ind]
        
        # now get element composition of species, columns 24:44
        # each piece of data is 5 characters long (2 for element, 3 for #)
        elem_str = split_str(line[24:44], 5)
        
        for e_str in elem_str:
            e = e_str[0:2].strip()
            # skip if blank
            if e == '': continue
            # may need to convert to float first, in case of e.g. "1."
            e_num = float( e_str[2:].strip() )
            e_num = int(e_num)
            
            spec.elem.append([e, e_num])
            
            # calculate molecular weight
            spec.mw += e_num * elem_mw[e.lower()]
        
        # temperatures for species
        T_spec = read_str_num(line[45:73])
        T_low  = T_spec[0]
        T_high = T_spec[1]
        if len(T_spec) == 3: T_com = T_spec[2]
        else: T_com = T_common[1]
        
        spec.Trange = [T_low, T_com, T_high]
        
        # second species line
        line = file.readline()
        coeffs = split_str(line[0:75], 15)
        spec.hi[0] = float( coeffs[0] )
        spec.hi[1] = float( coeffs[1] )
        spec.hi[2] = float( coeffs[2] )
        spec.hi[3] = float( coeffs[3] )
        spec.hi[4] = float( coeffs[4] )
        
        # third species line
        line = file.readline()
        coeffs = split_str(line[0:75], 15)
        spec.hi[5] = float( coeffs[0] )
        spec.hi[6] = float( coeffs[1] )
        spec.lo[0] = float( coeffs[2] )
        spec.lo[1] = float( coeffs[3] )
        spec.lo[2] = float( coeffs[4] )
        
        # fourth species line
        line = file.readline()
        coeffs = split_str(line[0:75], 15)
        spec.lo[3] = float( coeffs[0] )
        spec.lo[4] = float( coeffs[1] )
        spec.lo[5] = float( coeffs[2] )
        spec.lo[6] = float( coeffs[3] )
        
        # stop reading if all species in mechanism accounted for
        if not next((sp for sp in specs if sp.mw == 0.0), None): break
    
    file.close()
    return


def calc_spec_smh(T, specs):
    """Calculate standard-state entropies minus enthalpies for all species
    
    Input
    T: temperature
    specs: list of species (SpecInfo class)
    """
    
    spec_smh = []
    
    Tlog = math.log(T)
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T

    Thalf = T / 2.0
    T2 = T2 / 6.0
    T3 = T3 / 12.0
    T4 = T4 / 20.0

    for sp in specs:
        if T <= sp.Trange[1]:
            smh = sp.lo[0] * (Tlog - 1.0) + sp.lo[1] * Thalf + sp.lo[2] * T2 + \
                  sp.lo[3] * T3 + sp.lo[4] * T4 - (sp.lo[5] / T) + sp.lo[6]
        else:
            smh = sp.hi[0] * (Tlog - 1.0) + sp.hi[1] * Thalf + sp.hi[2] * T2 + \
                  sp.hi[3] * T3 + sp.hi[4] * T4 - (sp.hi[5] / T) + sp.hi[6]
        
        spec_smh.append(smh)
    
    return(spec_smh)


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


def convert_mech_irrev(mech_name, therm_name):
    """Convert Chemkin-style mechanism with reversible reactions.
    
    
    """
    import copy
    
    elems = []
    specs = []
    reacs = []
    
    # interpret reaction mechanism file
    [num_e, num_s, num_r, units] = read_mech(mech_name, elems, specs, reacs)
    
    # interpret thermodynamic database file
    read_thermo(therm_name, elems, specs)
    
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
    
    mod_mech = 'mech_irrev.dat'
    
    # write new reaction list to new file
    write_mech(mod_mech, elems, specs, reacs, units)
    
    return


if __name__ == "__main__":
    import sys
    convert_mech_irrev(sys.argv[1], sys.argv[2])
