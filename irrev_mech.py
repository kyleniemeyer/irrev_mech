#! /usr/bin/env python

import math

def read_string_num (string):
    """Pull list of numbers from a string"""
        
    # separate string into space-delimited strings of numbers
    num_str = string.split()
    
    nums = []
    for n in num_str:
        nums.append( float(n) )
    
    return nums

def read_thermo (elem_list, spec_list):
    """Read and interpret thermodynamic database for species
    
    Reads the file therm.dat and returns the species thermodynamic coefficients
    as well as the species-specific temperature range values (if given)
    
    Input
    elem_list: list of element names
    spec_list: list of species names
    """
    # make local copy of species list
    species = spec_list
    
    # initialize empty lists for element numbers, thermo coefficients, temperature ranges
    spec_elem = [[0] * len(elem_list)] * len(spec_list)
    spec_mw = [0.] * len(spec_list)
    spec_lo = [[0.] * 7]] * len(spec_list)
    spec_hi = [[0.] * 7]] * len(spec_list)
    
    # lines in thermo file are 80 characters long
    file = open('therm.dat', 'r')
    
    # loop through intro lines
    for line in file:
        
        # skip blank or commented lines
        if line == '\n' or line[0:1] == '!': continue
        
        # skip 'thermo' at beginning
        if line[0:6].lower() == 'thermo': break
    
    # next line has common temperature ranges
    line = file.readline()
    
    T_common = read_string_num(line)
    
    # now start reading species thermo info
    while len(species) > 0:
        # first line of species info
        line = file.readline()
        
        # break if end of file
        if line[0:3].lower() != 'end': break
        
        # species name, columns 0:18
        spec = line[0:18].strip()
        
        # now need to determine if this species is in mechanism
        spec_match = [x for x in species if (x in spec and len(x) == len(spec) )]
        
        # if not, read next three lines and continue
        if spec_match == []:
            line = file.readline()
            line = file.readline()
            line = file.readline()
            continue
        
        
        # get species index from list of species
        s = spec_list.index( spec )
        
        
    
    # remove processed species from list
    species.remove(spec)



def calc_spec_enthalpy (T, spec_mw, spec_lo, spec_hi, spec_T):
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
            h =  spec_lo[i][0] * T  + spec_lo[i][1] * T2 + spec_lo[i][2] * T3 + 
                 spec_lo[i][3] * T4 + spec_lo[i][4] * T5 + spec_lo[i][5]
            h = h * Ru / spec_mw[i]
            
        else:
            h =  spec_hi[i][0] * T  + spec_hi[i][1] * T2 + spec_hi[i][2] * T3 + 
                 spec_hi[i][3] * T4 + spec_hi[i][4] * T5 + spec_hi[i][5]
            h = h * Ru / spec_mw[i]
        
        spec_enthalpy.append(h)
    
    return spec_enthalpy


def calc_spec_entropy (T, spec_coeff, spec_coeff_T):
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
            s =  spec_lo[i][0] * Tl + spec_lo[i][1] * T  + spec_lo[i][2] * T2 + 
                 spec_lo[i][3] * T3 + spec_lo[i][4] * T4 + spec_lo[i][6]
            s = s * Ru / spec_mw[i]

        else:
            s =  spec_hi[i][0] * Tl + spec_hi[i][1] * T  + spec_hi[i][2] * T2 + 
                 spec_hi[i][3] * T3 + spec_hi[i][4] * T4 + spec_hi[i][6]
            s = s * Ru / spec_mw[i]

        spec_entropy.append(s)
    
    return spec_entropy

