# Module to keep the commonly-used functions for rheology related calculations
#
# Magali Billen, UCD 2019

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

# Half-space cooling model 
# t is seconds
# z in meters
def temperature_halfspace(z,t,p):
	# Physical constants
	kappa = 1e-6  # thermal diffusivity (m^2/s)
	T_s = 0  # surface temperature (C)
	T_m = 1350 # mantle temperature (C)

	T = T_s + (T_m - T_s)*erf(z/(2*np.sqrt(kappa*t)))
	
	if p == 1:
		fig = plt.figure()
		ax = fig.add_subplot(1,3,1)
		ax.plot(T,z/1000)  # plot in kilometers
		ax.set_ylim(150,0)
		ax.grid(True)
		ax.set_xlabel('Temperature (C)')
		ax.set_ylabel('Depth (km)')
	
	return T
	
# Pressure as a function of depth from the compressibility using Equation 4-321 
# (Turcotte and Schubert, Geodynamics 2nd edition, p. 190). 
# z in meters
def pressure_from_compessibility(z,p):
	# Physical constants
	beta = 4.3e-12  # compressiblity (1/Pa)
	rho_m = 3300    # reference density (kg/m^3)
	g = 9.81        # gravitational acceleration (m/s^2)
	
	P = (-1/beta)*np.log(1 - rho_m*g*beta*z)  # pressure (Pa)
	if p == 1:
		fig = plt.figure()
		ax = fig.add_subplot(1,3,1)
		ax.plot(P/1e9,z/1000,color='blue')
		ax.set_ylim(150,0)
		ax.grid(True)
		ax.set_xlabel('Pressure (GPa)')
		ax.set_ylabel('Depth (km)')
	return P

# Parameters from "Rheology of the Upper Mantle and the Mantle Wedge: A View from the 
# Experimentalists, Greg Hirth and David Kohlstedt, Inside the Subduction Factory, 
# Geophysical Monograph, 138 (2003)	
# T (K), P (Pa), d (microns), sigd (Pa), Coh (ppm-H/Si)

def diffusion_creep(T,P,d,sigd,Coh):
	# Physical Constants
	R = 8.314 # J/mol*K
	mpa = 1e6

	# Flow-law parameters (Hirth and Kohlstedt, wet olivine, 2003)
	n = 1
	p = 3
	r = 1.0
	A = 1.0e6/mpa	# for stress in Pa
	E = 335e3		# J/mol (+/- 75)
	V = 4.0e-6

	edot = A*d**(-p)*sigd**n*Coh**r*np.exp(-(E + P*V)/(R*T))

	return edot
    
# Parameters from "Rheology of the Upper Mantle and the Mantle Wedge: A View from the 
# Experimentalists, Greg Hirth and David Kohlstedt, Inside the Subduction Factory, 
# Geophysical Monograph, 138 (2003)	
# T (K), P (Pa), d (microns), sigd (Pa), Coh (ppm-H/Si)
# This uses version for hydrated samples with constant values of C_OH
# Fro dry conditions, A and E are different. 

def dislocation_creep(T,P,sigd,Coh):
	# Physical Constants
    R = 8.314 # J/mol*K
    mpa = 1e6
    # Flow-law parameters (Hirth and Kohlstedt, wet olivine, 2003)
    n = 3.5    
    r = 1.2
    A = 90/np.power(mpa,n)  # for stress in Pa
    E = 480e3 # J/mol (+/- 40)
    V = 11e-6  # m^3/mol
    
    edot = A*sigd**n*Coh**r*np.exp(-(E + P*V)/(R*T))
    
    return edot

# Values from  "Experimental constraints on the strength of the lithospheric mantle"
# by Mei, Kohlstedt, Dixon and Durham, JGR, v 115, 2010
# Applicable at  673 > T > 1273 K and anhydrous conditions, P=4-9 GPa   
def peierls_creep(T,sigd):
	# Physical Constants
    R = 8.314 # J/mol*K
    mpa = 1e6
    # Flow-law parameters (Hirth and Kohlstedt, wet olivine, 2003)
    n = 2.0
    p = 0.5
    A = 1.4e-7/(mpa**2)  # for stress in Pa  (s^-1 MPa-2)
    E = 320e3 # J/mol (+/- 50e3)
    sigp = 5.9e9  # Pa (+/-0.2 GPa)
        
    edot = A*sigd**n*np.exp(-(E/(R*T))*(1-(sigd/sigp)**p))
    
    return edot
	
