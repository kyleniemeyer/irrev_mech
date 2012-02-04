#! /usr/bin/env python

import math
import collections

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


def read_mech(filename, elem_list, spec_list, ):
    """Read and interpret mechanism file for elements, species, and reactions.
    
    Input
    filename:  reaction mechanism filename (e.g. 'mech.dat'
    """
    
    file = open(filename, 'r')
    
    return (num_e, num_s, num_r)


def read_thermo(filename, elem_list, spec_list, spec_elem, spec_mw, spec_hi, spec_lo):
    """Read and interpret thermodynamic database for species data.
    
    Reads the file therm.dat and returns the species thermodynamic coefficients
    as well as the species-specific temperature range values (if given)
    
    Input
    filename:  thermo database filename (e.g. 'therm.dat')
    elem_list: list of element names
    spec_list: list of species names
    spec_elem: list of elemental composition of each species (empty)
    spec_mw:   list of species molecular weights (empty)
    spec_hi:   list of species thermo coefficients (high temp)
    spec_lo:   list of species thermo coefficients (low temp)
    """
    # make local copy of species list
    species = list(spec_list)
    
    # lines in thermo file are 80 characters long
    file = open(filename, 'r')
    
    # loop through intro lines
    while True:
        line = file.readline()
    
        # skip blank or commented lines
        if line == '\n' or line[0:1] == '!': continue
    
        # skip 'thermo' at beginning
        if line[0:6].lower() == 'thermo': break
    
    # next line has common temperature ranges
    line = file.readline()
    
    T_ranges = read_str_num(line)
    
    # now start reading species thermo info
    while True:
        # first line of species info
        line = file.readline()
        
        # break if end of file
        if line[0:3].lower() == 'end': break
        if line is None: break
        
        # species name, columns 0:18
        spec = line[0:18].strip()
        
        # apparently in some cases notes are in the columns of shorter species names
        # so make sure no spaces
        if spec.find(' ') > 0:
            spec = spec[0 : spec.find(' ')]
        
        # now need to determine if this species is in mechanism
        spec_match = [x for x in species if (x in spec and len(x) == len(spec) )]
        
        # if not, read next three lines and continue
        if spec_match == []:
            line = file.readline()
            line = file.readline()
            line = file.readline()
            continue
        
        # get species index from list of species
        s = spec_list.index(spec)
        
        # now get element composition of species, columns 24:44
        # each piece of data is 5 characters long (2 for element, 3 for #)
        elem_str = split_str(line[24:44], 5)
        
        for e_str in elem_str:
            e = e_str[0:2].strip()
            # skip if blank
            if e == '': continue
            e_num = int( e_str[2:].strip() )
            
            e_ind = elem_list.index(e)
            spec_elem[s][e_ind] = e_num
            
            # calculate molecular weight
            spec_mw[s] += e_num * elem_mw[e]
        
        # temperatures for species
        T_spec = read_str_num(line[45:73])
        T_low  = T_spec[0]
        T_high = T_spec[1]
        if len(T_spec) == 3: T_com = T_spec[2]
        else: T_com = T_common[1]
        
        # second species line
        line = file.readline()
        coeffs = split_str(line[0:75], 15)
        spec_hi[s][0] = float( coeffs[0] )
        spec_hi[s][1] = float( coeffs[1] )
        spec_hi[s][2] = float( coeffs[2] )
        spec_hi[s][3] = float( coeffs[3] )
        spec_hi[s][4] = float( coeffs[4] )
        
        # third species line
        line = file.readline()
        coeffs = split_str(line[0:75], 15)
        spec_hi[s][5] = float( coeffs[0] )
        spec_hi[s][6] = float( coeffs[1] )
        spec_lo[s][0] = float( coeffs[2] )
        spec_lo[s][1] = float( coeffs[3] )
        spec_lo[s][2] = float( coeffs[4] )
        
        # fourth species line
        line = file.readline()
        coeffs = split_str(line[0:75], 15)
        spec_lo[s][3] = float( coeffs[0] )
        spec_lo[s][4] = float( coeffs[1] )
        spec_lo[s][5] = float( coeffs[2] )
        spec_lo[s][6] = float( coeffs[3] )
        
        # remove processed species from list
        species.remove(spec)
        
        # leave loop if all species in mechanism accounted for
        if species == []: break
    
    file.close()
    return


def calc_spec_enthalpy(T, spec_mw, spec_lo, spec_hi, spec_T):
    """Return species standard-state molar enthalpy
    
    Given species thermodynamic coefficients and species-specific temperature
    ranges (if given), calculate standard-state molar enthalpy for all species
    
    Input
    T: temperature
    spec_mw: list of species molecular weights
    spec_lo: list of species thermo coefficients (lower temp range)
    spec_hi: list of species thermo coefficients (higher temp range)
    spec_T: list of species-specific temperature range limits
    """
    
    spec_enthalpy = []
    
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T
    T5 = T4 * T
    
    T2 = T2 / 2.0
    T3 = T3 / 3.0
    T4 = T4 / 4.0
    T5 = T5 / 5.0
    
    for i in range( len(spec_mw) ):
        if T <= spec_T[i]:
            h =  spec_lo[i][0] * T  + spec_lo[i][1] * T2 + spec_lo[i][2] * T3 + \
                 spec_lo[i][3] * T4 + spec_lo[i][4] * T5 + spec_lo[i][5]
            h = h * Ru / spec_mw[i]
            
        else:
            h =  spec_hi[i][0] * T  + spec_hi[i][1] * T2 + spec_hi[i][2] * T3 + \
                 spec_hi[i][3] * T4 + spec_hi[i][4] * T5 + spec_hi[i][5]
            h = h * Ru / spec_mw[i]
        
        spec_enthalpy.append(h)
    
    return spec_enthalpy


def calc_spec_entropy(T, spec_coeff, spec_coeff_T):
    """Return species standard-state molar entropy
    
    Given species thermodynamic coefficients and species-specific temperature
    ranges (if given), calculate standard-state molar entropy for all species
    
    Input
    T: temperature
    spec_coeff: list of species thermo coefficients
    spec_coeff_T: list of species-specific temperature range limits
    """
    
    spec_entropy = []
    
    Tl = math.log(T)
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T

    T2 = T2 / 2.0
    T3 = T3 / 3.0
    T4 = T4 / 4.0

    for i in range( len(spec_mw) ):
        if T <= spec_T[i]:
            s =  spec_lo[i][0] * Tl + spec_lo[i][1] * T  + spec_lo[i][2] * T2 + \
                 spec_lo[i][3] * T3 + spec_lo[i][4] * T4 + spec_lo[i][6]
            s = s * Ru / spec_mw[i]

        else:
            s =  spec_hi[i][0] * Tl + spec_hi[i][1] * T  + spec_hi[i][2] * T2 + \
                 spec_hi[i][3] * T3 + spec_hi[i][4] * T4 + spec_hi[i][6]
            s = s * Ru / spec_mw[i]

        spec_entropy.append(s)
    
    return spec_entropy


def convert_mech_irrev(mech_name, therm_name):
    """Convert Chemkin-style mechanism with reversible reactions.
    
    
    """
    
    elem_list = []
    spec_list = []
    
    # interpret reaction mechanism file
    #[num_e, num_s, num_r] = read_mech(mech_name, elem_list, spec_list)
    
    elem_list = ['h', 'c', 'o', 'n', 'ar']
    spec_list = ['h', 'h2', 'o', 'o2', 'oh', 'h2o', 'n2', 'ho2', 'h2o2', 'ar']
    num_e = len(elem_list)
    num_s = len(spec_list)
    therm_name = 'therm.dat'
    
    # initialize empty lists for element numbers, thermo coefficients, temperature ranges
    spec_elem = [ [0 for j in range(num_e)] for i in range(num_s) ]
    spec_mw = [0.] * num_s
    spec_hi = [ [0.0 for j in range(7)] for i in range(num_s) ]
    spec_lo = [ [0.0 for j in range(7)] for i in range(num_s) ]
    
    # interpret thermodynamic database file
    read_thermo(therm_name, elem_list, spec_list, spec_elem, spec_mw, spec_hi, spec_lo)
    
    return


if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])
