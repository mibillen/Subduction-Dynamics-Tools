#!/opt/anaconda3/bin/python3
# Haoyuna: I create this file in order to keep track the rheology I use in 2D subduction in aspect

import numpy as np
import matplotlib.pyplot as plt

#Physical Constants
R=8.314 #J/mol*K

class RHEOLOGY_PRM():
    """
    class for rheologies
    """
    def __init__(self):
        '''
        Initiation, initiate rheology parameters
        '''
        # dislocation creep in Hirth & Kohlstedt 2003
        self.HK03_disl = \
            {
                "A": 90,
                "p": 0.0,
                "r": 1.2,
                "n": 3.5,
                "E": 480e3,
                "V": 11e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }

        # diffusion creep in Hirth & Kohlstedt 2003
        self.HK03_diff = \
            {
                "A" : 1.0e6,
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 335e3,
                "V" : 4e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }
        
        # dislocation creep in Arredondo & Billen 2017
        self.AB17_disl = \
            {
                "A": 2.57e-20,
                "p": 0.0,
                "r": 1.2,
                "n": 3.5,
                "E": 496e3,
                "V": 11e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }
        
        # diffusion creep in Arredondo & Billen 2017
        self.AB17_diff = \
            {
                "A" : 2.85e6,
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 317e3,
                "V" : 4e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }
        
        # modify dislocation creep in Hirth & Kohlstedt 2003
        self.HK03v1_disl = \
            {
                "A": 0.9,
                "p": 0.0,
                "r": 1.2,
                "n": 3.5,
                "E": 480e3,
                "V": 11e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }

        # diffusion creep in Hirth & Kohlstedt 2003
        self.HK03v1_diff = \
            {
                "A" : 1.0e6,
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 335e3,
                "V" : 4e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }
        
        # modified creep laws from Hirth & Kohlstedt 2003
        # for detail, refer to magali's explain_update_modHK03_rheology.pdf file
        self.HK03_wet_mod_diff = \
            {
                # "A" : 10**6.9,  # MPa^(-n-r)*um**p/s
                "A" : 7.1768e6,  # MPa^(-n-r)*um**p/s
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 375e3,
                "V" : 23e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet": 1.0  # I use this to mark this is a wet rheology, so I need to account for V and E for water later.
            }

        self.HK03_wet_mod_disl = \
            {
                "A" : 10**2.65,
                "p" : 0.0,
                "r" : 1.0,
                "n" : 3.5,
                "E" : 520e3,
                "V" : 24e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet" : 1.0
            }
        
        # modified creep laws from Hirth & Kohlstedt 2003
        # I bring the values to the limit of the range
        # for detail, refer to magali's explain_update_modHK03_rheology.pdf file
        self.HK03_wet_mod1_diff = \
            {
                # "A" : 10**6.9,  # MPa^(-n-r)*um**p/s
                "A" : 7.1768e6,  # MPa^(-n-r)*um**p/s
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 375e3 - 25e3,
                "V" : 23e-6 -5.5e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet": 1.0  # I use this to mark this is a wet rheology, so I need to account for V and E for water later.
            }

        self.HK03_wet_mod1_disl = \
            {
                "A" : 10**2.65,
                "p" : 0.0,
                "r" : 1.0,
                "n" : 3.5,
                "E" : 520e3 + 40e3,
                "V" : 24e-6 + 4e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet" : 1.0
            }
        
        
        self.water = \
            {
                "A" : 87.75,             # H/(10^6*Si)/MPa
                "E" : 50e3,                     # J/mol +/-2e3
                "V" : 10.6e-6                     # m^3/mol+/-1
            }


def GetRheology(rheology):
    '''
    read rheology parameters, and account for effects of water if it is a wet rheology
    '''
    RheologyPrm = RHEOLOGY_PRM()
    diffusion_creep = getattr(RheologyPrm, rheology + "_diff")
    dislocation_creep = getattr(RheologyPrm, rheology + "_disl")
    try:
        _ = diffusion_creep['wet']
    except KeyError:
        pass
    else:
        ### effects of water accounted, see Magali's file explain_update_modHK03_rheology eq(5)
        water_creep = getattr(RheologyPrm, "water")
        diffusion_creep['A'] = diffusion_creep['A'] / (water_creep['A'] ** diffusion_creep['r'])
        diffusion_creep['V'] = diffusion_creep['V'] - water_creep['V'] * diffusion_creep['r']
        diffusion_creep['E'] = diffusion_creep['E'] - water_creep['E'] * diffusion_creep['r']
        dislocation_creep['A'] = dislocation_creep['A'] / (water_creep['A'] ** dislocation_creep['r'])
        dislocation_creep['V'] = dislocation_creep['V'] - water_creep['V'] * dislocation_creep['r']
        dislocation_creep['E'] = dislocation_creep['E'] - water_creep['E'] * dislocation_creep['r']
    return diffusion_creep, dislocation_creep


def CreepRheology(creep_type, strain_rate, P, T, d=None, Coh=None, **kwargs):
    """
    def CreepRheology(creep_type, strain_rate, P, T, d, Coh):

    Calculate viscosity by flow law in form of (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - Coh: H / 10^6 Si
     - Return value: Pa*s
    Pay attention to pass in the right value, this custom is inherited
    """
    A = creep_type['A']
    p = creep_type['p']
    r = creep_type['r']
    n = creep_type['n']
    E = creep_type['E']
    V = creep_type['V']
    # compute value of F(pre factor)
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 1 / (2**((n-1)/n)*3**((n+2)/2/n))
    else:
        F = 1.0

    if d is None:
        d = creep_type['d']
    if Coh is None:
        Coh = creep_type['Coh']
    # calculate B
    B = A * d**(-p) * Coh**r
    eta = 1/2.0 * F * (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T)) * 1e6

    return eta