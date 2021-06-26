# -*- coding: utf-8 -*-
#!/opt/anaconda3/bin/python3
r"""Test for peierls rheology

This outputs:

  - test results to value of variable "test_dir"

This depends on:

  - source files written from the value of "source_dir"

Examples of usage:

  - default usage:

        python -m pytest 

descriptions:
    every function is a separate test, combined usage with pytest module
""" 

# makes  a deformation mechanism map in P-T space for constant grain size and strainrate

import os
import shutil
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
# import module from python scripts folder
import python_scripts.flow_law_functions as flf

test_dir = "test_output"
if not os.path.isdir(test_dir):
    os.mkdir(test_dir)
# source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_plot_run_time')

def test_plot_peierls():
    # Unit conversions
    mpa = 1e6  # MPa to Pa
    mum = 1e6  # microns/meter for grain size 
    R = 8.314       # J/mol*K
    km2m = 1e3 # km to meters
    sec2yrs = 60*60*24*365.25 # sec per year
    
    # Depth in meters
    zmin = 0
    zmax = 660
    
    depth = np.array([zmin, zmax])
    
    # Temperature range 
    mad = 0.5715  # deg/km, approximating adiabat to match T(CMB) in Aspect
    Tmin = 200
    Tmax = 2000
    dT = 50 # C
    
    Pmin = 0
    Pmax = 22 # GPa
    dP = 1 # GPa
    
    T1 = np.arange(Tmin,Tmax+dT,dT)
    P1 = np.arange(Pmin,Pmax+dP,dP)*1e9  # Pa
    
    T, P = np.meshgrid(T1,P1)
    
    # This is just here so I could look at what the pressure dependence gives for a change
    # in the peierls stress: goes from 5.9 to 8.6 GPa
    # From Kawazoe et al PEPI 2009
    G0 = 77.4*1e9 # GPa  
    Gp = 1.61 # GPa/GPa 
    sigp0 = 5.9e9
    sigp = sigp0*(1 + (Gp/G0)*P1)
    	
    # constant grain size and strain rate
    d = 10e3 	  # microns (= 10 mm = 1.0 cm)
    #edot = 1e-15  # background mantle strain-Rate
    edot = 10**((3/1000)*T - 19)  # haoyuan: are they dependent on temperature?
    edot1 = 10**((3/1000)*T1 - 19)
    
    # Get viscosity
    water = 'wet'
    flver = 'mod'
    coh = 1000  
    
    # Diffusion Creep
    A, Emid1, Vmid1, p, r = flf.get_diff_HK_params(water,flver,'mid','mid')    
    dm = d/1e6  # convert from microns to meters
    Am = A/1e6  # convert from MPa^-1 to Pa^-1
    
    dE = -50e3      # dE = +/- 75, error on activation energy, from HK03
    dV = -4.5e-6  	# dV = +/-4.5e-6  # Ohuchi et al, 2012	
    
    Edf = Emid1 + dE
    Vdf = Vmid1 + dV	
    fh2o = flf.convert_OH_fugacity(T,P,coh)	# Convert coh to water fugacity
    etadf = 0.5*(dm)**p/(Am*fh2o**r)*np.exp((Edf + P*Vdf)/(R*T)) # Pa s
    
    # Dislocation Creep
    A, Emid2, Vmid2, n, r = flf.get_disl_HK_params(water,flver,'mid','mid')    
    Am = A/(1e6)**n  # convert from MPa^-1 to Pa^-1
    
    dE = -25e3    # dE = +/-40e3  error on activation energy from HK03
    dV = 3e-6  	 # dV = +/-3e-6 Ohuchi et al, 2012
    lab1 = 'df: ' + str(Edf/1e3) + '/' + str(Vdf*1e6)
    
    Eds = Emid2 + dE
    Vds = Vmid2 + dV
    lab2 = 'ds: ' + str(Eds/1e3) + '/' + str(Vds*1e6)	
    					
    fh2o = flf.convert_OH_fugacity(T,P,coh)	# Convert coh to water fugacity
    etads =  0.5*edot**(1/n - 1)/(Am*fh2o**r)**(1/n)*np.exp((Eds + P*Vds)/(n*R*T)) # Pa s
    	
    # Peierls Creep
    etap = flf.peierls_approx_visc('MK10',0.17,P,T,edot)
    #etap1 = flf.peierls_approx_visc('MK10',0.17,P1,T1,edot)
    
    
    # Composite viscosity
    etacomp1 = etadf*etads*etap/(etadf*etads + etadf*etap + etads*etap)
    etacomp = etadf*etads/(etadf+ etads)
    sig = 2*etacomp*edot
    p = np.where(sig>100e6)
    #etacomp[p] = etacomp1[p]
    etacomp = etadf*etads*etap/(etadf*etads + etadf*etap + etads*etap)

    # subfig1: diffusion creep viscosity	
    vmin1 = 18
    vmax1 = 25		
    fig = plt.figure(figsize=(8.5,11.0))
    ax1 = fig.add_subplot(3,2,1)	
    c1 = ax1.pcolor(T,P/1e9,np.log10(etadf),cmap='viridis')
    c1.cmap.set_over('k')
    c1.set_clim(vmin1,vmax1)
    ax1.set_ylabel('Pressure (GPa)')
    ax1.set_title('Diffusion',fontsize=10)	
    fig.colorbar(c1,ax=ax1)

    # subfig2: dislocation creep viscosity, the strain rate is T-dependent
    ax2 = fig.add_subplot(3,2,2)	
    c2 = ax2.pcolor(T,P/1e9,np.log10(etads),cmap='viridis')
    c2.cmap.set_over('k')
    c2.set_clim(vmin1,vmax1)
    ax2.set_title('Dislocation',fontsize=10)
    fig.colorbar(c2,ax=ax2)
 
    # subfig3: peierls creep, strain rate is T-dependent
    ax3 = fig.add_subplot(3,2,3)	
    c3 = ax3.pcolor(T,P/1e9,np.log10(etap),cmap='viridis')
    ax3.contour(T,P/1e9,np.log10(etap),colors='white')
    c3.cmap.set_over('k')
    c3.set_clim(vmin1,vmax1)
    ax3.set_xlabel('Temperature (C)')
    ax3.set_title('Peierls',fontsize=10)	
    fig.colorbar(c3,ax=ax3)

    # subfig 4: composite viscosity 
    ax4 = fig.add_subplot(3,2,4)	
    c4 = ax4.pcolor(T,P/1e9,np.log10(etacomp),cmap='viridis')
    c4.cmap.set_over('k')
    c4.set_clim(vmin1,vmax1)
    ax4.contour(T,P/1e9,np.log10(etacomp),colors='white')
    ax4.set_xlabel('Temperature (C)')
    ax4.set_title('Composite',fontsize=10)
    fig.colorbar(c4,ax=ax4)

    # subfig5: stress, with the strain rate dependent on T and the computed composite viscosity
    ax5 = fig.add_subplot(3,2,5)	
    c5 = ax5.pcolor(T,P/1e9,sig/1e9,cmap='viridis')
    c5.cmap.set_over('k')
    c5.set_clim(0,0.1)
    ax5.set_xlabel('Temperature (C)')
    ax5.set_title('Stress df-ds only',fontsize=10)
    fig.colorbar(c5,ax=ax5)
    
    skp = 1
    if skp == 0:
    	ax6 = fig.add_subplot(3,2,6)	
    	c6 = ax6.pcolor(T,P/1e9,np.log10(edot),cmap='viridis')
    	ax6.set_xlabel('Temperature (C)')
    	ax6.set_title('Strain-rate (T)',fontsize=10)
    	fig.colorbar(c6,ax=ax6)
    
    plt.tight_layout()
    pdffile = os.path.join(test_dir, 'peierls_defmech.pdf')
    if os.path.isfile(pdffile):
      os.remove(pdffile)  # remove old file
    fig.savefig(pdffile,bbox_inches='tight')
    assert(os.path.isfile(pdffile))


def test_plot_peierls1():
    # Unit conversions
    mpa = 1e6  # MPa to Pa
    mum = 1e6  # microns/meter for grain size 
    R = 8.314       # J/mol*K
    km2m = 1e3 # km to meters
    sec2yrs = 60*60*24*365.25 # sec per year
    
    # Depth in meters
    zmin = 0
    zmax = 660
    
    depth = np.array([zmin, zmax])
    
    # Temperature range 
    mad = 0.5715  # deg/km, approximating adiabat to match T(CMB) in Aspect
    Tmin = 200
    Tmax = 2000
    dT = 50 # C
    
    Pmin = 0
    Pmax = 22 # GPa
    dP = 1 # GPa
    
    T1 = np.arange(Tmin,Tmax+dT,dT)
    P1 = np.arange(Pmin,Pmax+dP,dP)*1e9  # Pa
    
    T, P = np.meshgrid(T1,P1)
    
    # This is just here so I could look at what the pressure dependence gives for a change
    # in the peierls stress: goes from 5.9 to 8.6 GPa
    # From Kawazoe et al PEPI 2009
    G0 = 77.4*1e9 # GPa  
    Gp = 1.61 # GPa/GPa 
    sigp0 = 5.9e9
    sigp = sigp0*(1 + (Gp/G0)*P1)
    	
    # constant grain size and strain rate
    d = 10e3 	  # microns (= 10 mm = 1.0 cm)
    edot = 1e-15  # background mantle strain-Rate
    edot1 = 1e-15  # background mantle strain-Rate
    
    # Get viscosity
    water = 'wet'
    flver = 'mod'
    coh = 1000  
    
    # Diffusion Creep
    A, Emid1, Vmid1, p, r = flf.get_diff_HK_params(water,flver,'mid','mid')    
    dm = d/1e6  # convert from microns to meters
    Am = A/1e6  # convert from MPa^-1 to Pa^-1
    
    dE = -50e3      # dE = +/- 75, error on activation energy, from HK03
    dV = -4.5e-6  	# dV = +/-4.5e-6  # Ohuchi et al, 2012	
    
    Edf = Emid1 + dE
    Vdf = Vmid1 + dV	
    fh2o = flf.convert_OH_fugacity(T,P,coh)	# Convert coh to water fugacity
    etadf = 0.5*(dm)**p/(Am*fh2o**r)*np.exp((Edf + P*Vdf)/(R*T)) # Pa s
    
    # Dislocation Creep
    A, Emid2, Vmid2, n, r = flf.get_disl_HK_params(water,flver,'mid','mid')    
    Am = A/(1e6)**n  # convert from MPa^-1 to Pa^-1
    
    dE = -25e3    # dE = +/-40e3  error on activation energy from HK03
    dV = 3e-6  	 # dV = +/-3e-6 Ohuchi et al, 2012
    lab1 = 'df: ' + str(Edf/1e3) + '/' + str(Vdf*1e6)
    
    Eds = Emid2 + dE
    Vds = Vmid2 + dV
    lab2 = 'ds: ' + str(Eds/1e3) + '/' + str(Vds*1e6)	
    					
    fh2o = flf.convert_OH_fugacity(T,P,coh)	# Convert coh to water fugacity
    etads =  0.5*edot**(1/n - 1)/(Am*fh2o**r)**(1/n)*np.exp((Eds + P*Vds)/(n*R*T)) # Pa s
    	
    # Peierls Creep
    etap = flf.peierls_approx_visc('MK10',0.17,P,T,edot)
    #etap1 = flf.peierls_approx_visc('MK10',0.17,P1,T1,edot)
    
    
    # Composite viscosity
    etacomp1 = etadf*etads*etap/(etadf*etads + etadf*etap + etads*etap)
    etacomp = etadf*etads/(etadf+ etads)
    sig = 2*etacomp*edot
    p = np.where(sig>100e6)
    #etacomp[p] = etacomp1[p]
    etacomp = etadf*etads*etap/(etadf*etads + etadf*etap + etads*etap)

    # subfig1: diffusion creep viscosity	
    vmin1 = 18
    vmax1 = 25		
    fig = plt.figure(figsize=(8.5,11.0))
    ax1 = fig.add_subplot(3,2,1)	
    c1 = ax1.pcolor(T,P/1e9,np.log10(etadf),cmap='viridis')
    c1.cmap.set_over('k')
    c1.set_clim(vmin1,vmax1)
    ax1.set_ylabel('Pressure (GPa)')
    ax1.set_title('Diffusion',fontsize=10)	
    fig.colorbar(c1,ax=ax1)

    # subfig2: dislocation creep viscosity, the strain rate is T-dependent
    ax2 = fig.add_subplot(3,2,2)	
    c2 = ax2.pcolor(T,P/1e9,np.log10(etads),cmap='viridis')
    c2.cmap.set_over('k')
    c2.set_clim(vmin1,vmax1)
    ax2.set_title('Dislocation',fontsize=10)
    fig.colorbar(c2,ax=ax2)
 
    # subfig3: peierls creep, strain rate is T-dependent
    ax3 = fig.add_subplot(3,2,3)	
    c3 = ax3.pcolor(T,P/1e9,np.log10(etap), cmap='viridis')
    # ax3.contour(T,P/1e9,np.log10(etap), [21.0, 22.0], colors='white')
    ax3.contour(T,P/1e9,np.log10(etap), [23.0], colors='white')
    c3.cmap.set_over('k')
    c3.set_clim(vmin1,vmax1)
    ax3.set_xlabel('Temperature (C)')
    ax3.set_title('Peierls',fontsize=10)	
    fig.colorbar(c3,ax=ax3)

    # subfig 4: composite viscosity 
    ax4 = fig.add_subplot(3,2,4)	
    c4 = ax4.pcolor(T,P/1e9,np.log10(etacomp),cmap='viridis')
    c4.cmap.set_over('k')
    c4.set_clim(vmin1,vmax1)
    ax4.contour(T,P/1e9,np.log10(etacomp), [21.0, 22.0],colors='white')
    ax4.set_xlabel('Temperature (C)')
    ax4.set_title('Composite',fontsize=10)
    fig.colorbar(c4,ax=ax4)

    # subfig5: stress, with the strain rate dependent on T and the computed composite viscosity
    ax5 = fig.add_subplot(3,2,5)	
    c5 = ax5.pcolor(T,P/1e9,sig/1e9,cmap='viridis')
    c5.cmap.set_over('k')
    c5.set_clim(0,0.1)
    ax5.set_xlabel('Temperature (C)')
    ax5.set_title('Stress df-ds only',fontsize=10)
    fig.colorbar(c5,ax=ax5)
    
    skp = 1
    if skp == 0:
    	ax6 = fig.add_subplot(3,2,6)	
    	c6 = ax6.pcolor(T,P/1e9,np.log10(edot),cmap='viridis')
    	ax6.set_xlabel('Temperature (C)')
    	ax6.set_title('Strain-rate (T)',fontsize=10)
    	fig.colorbar(c6,ax=ax6)
    
    plt.tight_layout()
    pdffile = os.path.join(test_dir, 'peierls_defmech_1e-15.pdf')
    if os.path.isfile(pdffile):
      os.remove(pdffile)  # remove old file
    fig.savefig(pdffile,bbox_inches='tight')
    assert(os.path.isfile(pdffile))


def test_plot_peierls_hansen19():
    # test using implementation from Hansen etal 2019
    # this function reproduces fig 8 from the paper
    P = 0.0  # no dependence on P
    T = 298  # room temperature
    edot = 1e-5  # guessed
    ds = 10**np.linspace(-2, 5, 100) * 1e-6
    viscs = flf.peierls_hansen19_visc(P,T,edot,ds, sigma_b = 1.8e9)
    stresses = 2 * viscs * edot
    std = np.array([[-0.9861750909741067, 0.9732704402515723],
     [-0.8179781662144259, 0.9150943396226414],
     [-0.6329799541448438, 0.8522012578616351],
     [-0.38082917902442137, 0.7735849056603773],
     [-0.19598872552112967, 0.7201257861635219],
     [-0.01120085820660055, 0.6698113207547169],
     [0.22385940556572148, 0.6132075471698113],
     [0.475589491176037, 0.5597484276729559],
     [0.6768368355735053, 0.5251572327044025],
     [0.8947014156202009, 0.4968553459119496],
     [1.0622935992091032, 0.4748427672955975],
     [1.246555604635999, 0.4559748427672956],
     [1.4977598283586788, 0.4339622641509433],
     [1.8492985002418951, 0.4119496855345912],
     [2.1170673734250434, 0.39937106918238996],
     [2.4015060684461824, 0.3899371069182389],
     [2.736059401358826, 0.3836477987421385],
     [2.987053280326454, 0.37421383647798745],
     [3.2880566248080596, 0.37421383647798745],
     [3.62255737153194, 0.371069182389937],
     [3.890168486148797, 0.3679245283018868]])
    # plot fig 8
    fig = plt.figure(figsize=(5, 10))
    ax1 = fig.add_subplot(2,1,1)
    ax1.loglog(ds/1e-6, stresses/1e9, 'b', label='test, edot=%.4e' % edot)
    ax1.loglog(10.0**std[:, 0], 10.0**std[:, 1], 'k--', label='standard')
    ax1.grid()
    ax1.set_xlim([1e-2, 1e5])
    ax1.set_ylim([1.0, 30.0])
    ax1.set_xlabel('Length scale(um)')
    ax1.set_ylabel('Differential stress(GPa)')
    ax1.legend()
    ax1.set_title('fig 8')
    pdffile = os.path.join(test_dir, 'peierls_hansen19.pdf')
    # fig 10,c (temperature dependence)
    P = 0
    Ts = np.linspace(200.0, 1000.0, 100)
    d = 0.01  # 1 cm
    viscs = flf.peierls_hansen19_visc(P,Ts,edot,d, sigma_b=0.0)  # sigma_b = 0, corresponding to epsl = 0.0
    stresses = 2 * viscs * edot
    std = np.array([[344.146218960934, 2.2154696132596685],
    [576.6130514142911, 1.541436464088398],
    [786.6619299040804, 0.9558011049723758],
    [1045.8894292609425, 0.22651933701657434]])
    # plot fig 10, c
    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(Ts, stresses/1e9, 'b', label='test, edot=%.4e' % edot)
    ax2.plot(std[:, 0], std[:, 1], 'k--', label='standard')
    ax2.set_xlim([0, 1200.0])
    ax2.set_ylim([0, 3.0])
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Differential stress(GPa)')
    ax2.legend()
    ax2.set_title('fig 9')
    fig.savefig(pdffile, bbox_inches='tight')
    assert(os.path.isfile(pdffile))