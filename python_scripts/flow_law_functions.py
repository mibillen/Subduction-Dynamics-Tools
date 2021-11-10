#!/opt/anaconda3/bin/python3
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import python_scripts.TwoD_haoyuan_2021.rheology_haoyuan_2021 as rheology_haoyuan_2021

output_dir = "output"
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

#Physical Constants
R=8.314 #J/mol*K

# Input parameter choices as:
# T = temperature               (K)
# P = pressure                (Pa)
# d = grainsize               (microns) -> converted to meters inside.
# sigd = differential stress  (MPa)
# coh = water content         (ppm)
# fh2o = MPa

# Edev, Vdev:  "mid", "min" or "max" value from reported error range.
# mod:  "orig" or "new" refers to using the orignal values of V (and A) from HK03 or new values 
# following Ohuchi 2015.

# Most flow laws require the water fugacity rather than the water content, but we usually
# define water as C_OH, so need to convert this to pass to flow-laws.
# ** Exception is the original HK03 "Constant C_OH" flow laws which already did this conversion.

# examples:
#    plot peierls rheology with dislocation creep, Mk 10 and hansen 19 are included
# 	python -m python_scripts.flow_law_functions plot_hansen19_dislocation
 
def convert_OH_fugacity(T,P,coh):

    # T [K]             Temperature
    # P [MPa]            Pressure
    # coh [ppm-H/si]    OH content (note can convert C_h2o = Coh/16.25, confirmed by Ohuchi)
    
    # Per e-mail from Ohuchi, use data from the table in Keppler & Bolfan-Casanova, 2006
    # to related C_H2O to f_H2O, checked against Kohlstedt et al 1996 data. 
    # I used Kohlstedt et al., 1996 data to show that CH2O = COH/16.25
    # substituting this leads to multiply Ah2o by 16.25Then change units
    # bars to MPa... the resulting value (87.75 ppm-H/Si/MPAa) is very close to the 
    # Ah2o in Zhao (90 # ppm-H/Si/MPa +/- 10) but we ignore the iron content dependence.  
       
    Ah2o = 0.54*(10)*16.25             # (ppm/bar)*(MPa/bar)*COH/CH2O
    Eh2o = 50e3                     # J/mol +/-2e3
    Vh2o = 10.6e-6                     # m^3/mol+/-1
    fH2O = coh/(Ah2o*np.exp(-(Eh2o + P*Vh2o)/(R*T)))       #MPa
    
    return fH2O

# Use choices of water (wet/dry/constant), mod (orig,modified) and deviation (mid, min, max)
# to set the values of A, E and V to use in the strain rate or viscosity calculation
# functions.

# ** I think this information might be better kept in a dictionary?? ** 
# wet or dry or con
# uses "mid" reported values, 
# and define dE or dV, which is +/- below
       
def get_diff_HK_params(water,mod,Edev,Vdev):
    p = 3 
    
    if water == 'dry':
        r = 0
        if mod == 'orig':
            A = 1.5e9/(1e6)**p # MPa^(-n-r) * m^p * s^-1 (converted from microns to meters)
            E = 375e3      # J/mol 
            V = 6e-6       # m^3/mol, Mid value form range in Table 1 of HK 03.
            dE = 50e3      # error on activation energy
            dV = 4e-6      # error on activation volume 
        else: 
            A = 10**(-10.40) #MPa^(-n-r) * m^p * s^-1  Hansen et al (2011), Table A2 logAdif=7.6 MPa, microns, K; grain size difference; convert microns to meters 
            E = 375e3       # J/mol     
            V = 8.2e-6      # m^3/mol       Nishihara et al. (2014) for forsterite   8.2 +/- 0.9 cm^3.mol
            dE = 50e3          # error on activation energy from HK03
            dV = 0.9e-6      # error on activation volume from Nishihara         
    elif water == 'wet':      # HK03 wet diffusion, (NOT constant COH)
        r = 1
        if mod == 'orig':
            A = 10**(-10.6) # MPa^(-n-r) * m^p * s^-1 # HK03: 2.5e7 converted from microns to meters        
            E = 375e3         # J/mol                                    
            V = 10e-6          # m^/mol, from HK03   
            dE = 75          # error on activation energy, from HK03
            dV = 10e-6      # error on activation volume, from HK03  
        else: 
            A = (10**(-10.6))/3.5 #MPa^(-n-r) * m^p * s^-1 # 2.5e7 converted from microns to meters; bell03 water correction        
            E = 375e3  # J/mol                                    
            V = 23e-6  # m^/mol, Ohuchi et al. (2012)   
            dE = 75e3    # error on activation energy from HK03
            dV = 4.5e-6  # Ohuchi et al, 2012
    elif water == 'con': 
        print('Using Diffusion Constant-COH')     
        r = 1        
        A = 1.0e6/(1e6)**p # Pa^(-n-r) * m^p * s^-1 (converted from microns to meters)      
        E = 335e3         # J/mol                                    
        V = 4e-6          # m^/mol, from HK03   
        dE = 75e3          # error on activation energy, from HK03
        dV = 4e-6          # use error for from HK03  
        if mod != 'orig':
            print("Error no modified versions for Constant-COH HK03, using originals")
    else:
        print("water must be dry, wet or con")
                    
    if Edev == 'min':  
        E = E - dE  
    elif Edev == 'max':
        E = E + dE  

    if Vdev == 'min':  
        V = V - dV  
    elif Vdev == 'max':
        V = V + dV
        
    return(A,E,V,p,r)
    
    
def get_disl_HK_params(water,mod,Edev,Vdev):
    p = 0   # no grain size dependence
    
    if water == 'dry':
        r = 0
        n = 3.5
        dn = 0.3
        if mod == 'orig':
            A = 1.1e5 # MPa^(-n-r) * s^-1 (=10^5.04)
            E = 530e3      # J/mol 
            V = 18e-6   # m^3/mol, Mid value from range in Table 2 of HK 03.
            dE = 4e3      # error on activation energy
            dV = 4.0e-6      # error on activation volume 
        else: 
            A = 1.1e5 # MPa^(-n-r) * s^-1 (=10^5.04)
            E = 530e3      # J/mol 
            V = 17e-6   # m^3/mol, from Kawazoe etal., 2009            
            dE = 4e3      # error on activation energy
            dV = 2.5e-6      # error on activation volume     
    elif water == 'wet':      # HK03 wet diffusion, (NOT constant COH)
        r = 1.2
        dr = 0.4
        n = 3.5
        dn = 0.3
        if mod == 'orig':
            A = 1600         # MPa^(-n-r) # HK03 (=10^3.2)
            E = 520e3         # J/mol                                    
            V = 22e-6          # m^/mol, from HK03   
            dE = 40e3          # error on activation energy, from HK03
            dV = 11e-6      # error on activation volume, from HK03  
        else: 
            A = 1600/3.5 #MPa^(-n-r) s^-1 # HK03 modified for Bell_etal 03: (=10^2.66)
            E = 520e3  # J/mol                                    
            V = 24e-6  # m^/mol, Ohuchi et al. (2012)   
            dE = 40e3    # error on activation energy from HK03
            dV = 3e-6  # NEED TO CHECK Ohuchi et al, 2012
    elif water == 'con': 
        print('Using Dislocation Constant-COH')     
        r = 1.2 
        dr = 0.4
        n = 3.5
        dn = 0.3       
        A = 90 # MPa^(-n-r) * s^-1     
        E = 480e3         # J/mol                                    
        V = 11e-6          # m^/mol, from HK03   
        dE = 40e3          # error on activation energy, from HK03
        dV = 5e-6          # use error for from HK03  
        if mod != 'orig':
            print("Error no modified versions for Constant-COH HK03, using originals")
    else:
        print("water must be dry, wet or con")
                    
    if Edev == 'min':  
        E = E - dE  
    elif Edev == 'max':
        E = E + dE  

    if Vdev == 'min':  
        V = V - dV  
    elif Vdev == 'max':
        V = V + dV
        
    return(A,E,V,n,r)

# Preferred dislocation creep parameters from Kawazoe 
def get_disl_KAW09_params(water,mod,Edev,Vdev):
    p = 0   # no grain size dependence
    
    if water == 'dry':  # same as modified HK03
        r = 0
        n = 3.5
        dn = 0.3
        A = 1.1e5 # MPa^(-n-r) * s^-1Scaling of edot/A * np.exp(Delta_F/(R*T)) term in the sinh function  (=10^5.04)
        E = 530e3      # J/mol 
        V = 17e-6   # m^3/mol, from Kawazoe etal., 2009            
        dE = 4e3      # error on activation energy
        dV = 2.5e-6      # error on activation volume    
    elif water == 'wet':      # HK03 wet diffusion, (NOT constant COH)
        r = 1.2
        dr = 0.4
        n = 3.0
        dn = 0.1
        A = 10**2.9 #MPa^(-n-r) s^-1 # HK03 modified for Bell_etal 03: (=10^2.66)
        E = 470e3  # J/mol                                    
        V = 24e-6  # m^/mol,  
        dE = 40e3  #
        dV = 3e-6  # 
    else:
        print("water must be dry or wet")
                    
    if Edev == 'min':  
        E = E - dE  
    elif Edev == 'max':
        E = E + dE  

    if Vdev == 'min':  
        V = V - dV  
    elif Vdev == 'max':
        V = V + dV
    return(A,E,V,n,r)    
    
# Olivine diffusion creep: strain-rate
def edot_diff_HK(T,P,d,sigd,coh,water,mod,Edev,Vdev):   

    A, E, V, p, r = get_diff_HK_params(water,mod,Edev,Vdev)    

    if water == 'con': # use constant COH equation/values from HK03
        edot = A*(sigd)**n*(d/1e6)**(-p)*((coh)**r)*np.exp(-(E+P*V)/(R*T)) #s^-1
    else:
        fh2o = convert_OH_fugacity(T,P,coh) # Convert coh to water fugacity
        edot = A*(sigd)**n*(d/1e6)**(-p)*((fh2o)**r)*np.exp(-(E+P*V)/(R*T)) #s^-1

    return edot
    
# Olivine dislocation creep: strain-rate
def edot_disc_HK(T,P,sigd,coh,water,mod,Edev,Vdev):   
    
    A, E, V, n, r = get_disc_HK_params(water,mod,Edev,Vdev)    

    if water == 'con': # use constant COH equation/values from HK03
        edot = A*(sigd)**n*((coh)**r)*np.exp(-(E+P*V)/(n*R*T)) #s^-1
    else:
        fh2o = convert_OH_fugacity(T,P,coh) # Convert coh to water fugacity
        edot = A*(sigd)**n*((fh2o)**r)*np.exp(-(E+P*V)/(n*R*T)) #s^-1

    return edot
    
# Olivine dislocation creep: strain-rate
def edot_disc_KAW09(T,P,sigd,coh,water,mod,Edev,Vdev):   
    
    A, E, V, n, r = get_disc_KAW09_params(water,mod,Edev,Vdev)    

    fh2o = convert_OH_fugacity(T,P,coh) # Convert coh to water fugacity
    edot = A*(sigd)**n*((fh2o)**r)*np.exp(-(E+P*V)/(n*R*T)) #s^-1

    return edot

# Olivine diffusion creep: viscosity  
# Note that A in this case needs to be convert from MPa to Pa to give visc in Pa-s
# instead of MPa-s 
def visc_diff_HK(T,P,d,coh,water,mod,Edev,Vdev):   

    A, E, V, p, r = get_diff_HK_params(water,mod,Edev,Vdev)    
    dm = d/1e6  # convert from microns to meters
    Am = A/1e6  # convert from MPa^-1 to Pa^-1
        
    if water == 'con': # use constant COH equation/values from HK03        
        visc = 0.5*(dm)**p/(Am*coh**r)*np.exp((E + P*V)/(R*T)) # Pa s
    else:     
        fh2o = convert_OH_fugacity(T,P,coh)    # Convert coh to water fugacity
        visc = 0.5*(dm)**p/(Am*fh2o**r)*np.exp((E + P*V)/(R*T)) # Pa s

    return visc
    
# Olivine dislocation creep: viscosity  
# Note that A in this case needs to be convert from MPa to Pa to give visc in Pa-s
# instead of MPa-s 
def visc_disl_HK(T,P,edot,coh,water,mod,Edev,Vdev):   
    
    A, E, V, n, r = get_disl_HK_params(water,mod,Edev,Vdev)    
    Am = A/(1e6)**n  # convert from MPa^-1 to Pa^-1
        
    if water == 'con': # use constant COH equation/values from HK03        
        visc = 0.5*edot**(1/n - 1)/(Am*coh**r)**(1/n)*np.exp((E + P*V)/(n*R*T)) # Pa s
    else:     
        fh2o = convert_OH_fugacity(T,P,coh)    # Convert coh to water fugacity
        visc =  0.5*edot**(1/n - 1)/(Am*fh2o**r)**(1/n)*np.exp((E + P*V)/(n*R*T)) # Pa s

    return visc

# Olivine dislocation creep: viscosity  
# Note that A in this case needs to be convert from MPa to Pa to give visc in Pa-s
# instead of MPa-s 
# For Kawazoe et al., 2009 preferred values    
def visc_disl_KAW09(T,P,edot,coh,water,mod,Edev,Vdev):   
    
    A, E, V, n, r = get_disl_KAW09_params(water,mod,Edev,Vdev)    
    Am = A/(1e6)**n  # convert from MPa^-1 to Pa^-1
        
    fh2o = convert_OH_fugacity(T,P,coh)    # Convert coh to water fugacity
    visc =  0.5*edot**(1/n - 1)/(Am*fh2o**r)**(1/n)*np.exp((E + P*V)/(n*R*T)) # Pa s

    return visc

# Peierls creep flow law 
# flv: flow law version
# MK10: for Mei and Kohlstedt, 2010     (gam = 0.17)
# gam: fitting parameter = sig_ref/sigp
#
def peierls_approx_visc(flv,gam,P,T,edot):

    mpa = 1e6  # MPa to Pa

    if flv == "MK10":
        # Mei et al., JGR 2010
        q = 1.0
        p = 0.5
        n = 2.0
        sigp0 = 5.9e9                    # Pa (+/- 0.2e9 Pa)
        A = 1.4e-7/np.power(mpa,n)     # s^-1 Pa^-2
        E = 320e3                      # J/mol (+/-50e3 J/mol)
    
    # Pressure dependence of Peierls stress
    # From Kawazoe et al. PEPI 2009 and parameters from Liu et al., GRL 2005 
    G0 = 77.4*1e9 # GPa  
    Gp = 1.61 # GPa/GPa 
    sigp = sigp0*(1 + (Gp/G0)*P)
    
    s = (E/(R*T))*p*q*((1-gam**p)**(q-1))*(gam**p)
    x = 1/(s+n)
    visc = (0.5*gam*sigp*edot**(x-1)) / ((A*(gam*sigp)**n)**x) * np.exp( (E*(1-gam**p)**q)/(R*T*(s+n)) ) 
    
    return visc


# Peierls creep flow law 
# flv: flow law version
# This form of rheology is from Hansen etal 2019
# Every inputs need to be in U/I
# P dependence is not implemented yet
def peierls_hansen19_visc(T,edot,d, **kwargs):    
    R = 8.314
    A = kwargs.get('A', 10**20.7)
    Delta_F = kwargs.get('Delta_F', 550e3)  # change in F
    sigma_l = kwargs.get('sigma_l', 3.1e9)  # internal resistance of the lattice to dislocation motion
    K = kwargs.get('K', 3.2e9*(1e-6)**0.5)  # material-dependent constant
    gamma = kwargs.get('gamma', 75.0)  # rate constant
    sigma_b_max = kwargs.get('sigma_b_max', 1.8e9)  # maximum back stress
    # try reading in sigma_l
    try:
        sigma_b = kwargs['sigma_b']
    except KeyError:
        ep = kwargs.get('ep', 0.0)  # strain
        sigma_b = sigma_b_max * (1 - exp(-gamma*ep))
    # derive Sigma
    Sigma = sigma_l + K * d**(-0.5)
    # derive sigma(differential stress)
    sigma = R * T * Sigma / Delta_F * np.arcsinh(edot/A * np.exp(Delta_F/(R*T))) + sigma_b
    # derive viscosity
    visc = sigma / (2*edot)
    return visc


# plot peierls rheology with dislocation creep
# Mk 10 and hansen 19 are included
# MK 10 has a free parameter gam
# Hansen 19 has a free parameter sigma_b, which is actually strain-dependent
def plot_peierls_hansen19_visc(**kwargs):
    def plot_temperature_grain_size(ax, P, T, d, edots, dislocation_creep, gam):
	# define this function to plot some temperature and grain size
        eta_disls = rheology_haoyuan_2021.CreepRheology(dislocation_creep, edots, P, T, use_effective_strain_rate=True)
        stress_disls = 2 * eta_disls * edots
        ax.loglog(edots, stress_disls, 'b--', label='dislocation creep')  # plot
        # compute viscosity from Hansoen 1
        sigma_b = 0.0
        eta_peierls_0 = peierls_hansen19_visc(T,edots,d, sigma_b=sigma_b)  # sigma_b = 0, corresponding to epsl = 0.0
        stress_peierls_0 = 2 * eta_peierls_0 * edots
        ax.loglog(edots, stress_peierls_0, 'r-', label='Hansen 19, sigma_b = %.1e Pa' % sigma_b)
        sigma_b = 1e9
        eta_peierls_1 = peierls_hansen19_visc(T,edots,d, sigma_b=sigma_b)  # sigma_b = 1 GPa, middle value
        stress_peierls_1 = 2 * eta_peierls_1 * edots
        ax.loglog(edots, stress_peierls_1, 'r--', label='Hansen 19, sigma_b = %.1e Pa' % sigma_b)
        sigma_b = 1.8e9
        eta_peierls_2 = peierls_hansen19_visc(T,edots,d, sigma_b=sigma_b)  # sigma_b = 1.8 GPa, largest value
        stress_peierls_2 = 2 * eta_peierls_2 * edots
        ax.loglog(edots, stress_peierls_2, 'r-.', label='Hansen 19, sigma_b = %.1e Pa' % sigma_b)
        # compute viscosity from Mei & K 10
        eta_MKs = peierls_approx_visc('MK10',gam,P,T,edots)
        stress_MKs = 2 * eta_MKs * edots
        ax.loglog(edots, stress_MKs, 'g-', label='M&K 10, gam = %.2f' % gam)
        # plot
        ax.set_xlabel('Strain rate [s^-1]')
        ax.set_ylabel('Stress [Pa]')
        ax.set_title('T = %.1f K, d = %.1e' % (T, d))
        ax.grid()
        ax.legend()

    # compute dislocation creep rheology
    R = 8.314
    P = 1e9
    edots = 10**np.linspace(-18, -11, 100)
    # get rheology for modified wet rheology
    # with a little deviation from the central values(E, V)
    dEdiff = -25e3  # -50 - 50e3
    dEdisl = 40.0e3  # -40 - 40e3
    dVdiff = -5.5e-6  # -5.5 - 5.5e-6
    dVdisl = 4.0e-6  # -4.0 - 4.0e-6
    rheology = kwargs.get('rheology', 'HK03_wet_mod')
    sigma_b = kwargs.get('sigma_b', 0.0)
    gam = 0.17
    diffusion_creep, dislocation_creep = rheology_haoyuan_2021.GetRheology(rheology)
    diffusion_creep['E'] += dEdiff
    dislocation_creep['E'] += dEdisl
    diffusion_creep['V'] += dVdiff
    dislocation_creep['V'] += dVdisl
    fig = plt.figure(figsize=(15, 15))  # initiate plot
    
    axs = fig.subplots(3,3)
    # 700 K
    d = 1e-2
    T = 700
    plot_temperature_grain_size(axs[0,0], P, T, d, edots, dislocation_creep, gam)

    # 900 K
    T = 900
    plot_temperature_grain_size(axs[0,1], P, T, d, edots, dislocation_creep, gam)
    
    # 1100 K
    T = 1100
    plot_temperature_grain_size(axs[0,2], P, T, d, edots, dislocation_creep, gam)
    
    # save figure
    fig.tight_layout()
    pdffile = os.path.join(output_dir, 'hansen19_dislocation_creep.pdf')
    fig.savefig(pdffile, bbox_inches='tight')
    

def main():
    '''
    main function of this module
    Inputs:
        sys.arg[1](str):
            commend
        sys.arg[2, :](str):
            options
    '''
    _commend = sys.argv[1]
    # parse options
    parser = argparse.ArgumentParser(description='Parse parameters')
    parser.add_argument('-i', '--inputs', type=str,
                        default='',
                        help='Some inputs')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'plot_hansen19_dislocation':
        # plot peierls rheology with dislocation creep
	# Mk 10 and hansen 19 are included	    
        plot_peierls_hansen19_visc(sigma_b=1e9)

# run script
if __name__ == '__main__':
    main()