import numpy as np
from math import *
#import sys

#Constants
h = 1.054571596e-34           #Planck's constant over 2*pi [J-s]
hh = h*2*pi
heV = 6.5821192815e-16        #[eV-s]
kb = 1.380648813e-23          # Boltzmann constant [J/K]
kbeV = 8.617332478e-5         # Boltzmann constant [eV/K]
ee = 1.60217656535e-19        # elementary electric charge [C] or [J/ev]
ep0 = 8.854187817e-12         # electric permittivity
c0 = 2.99792458e+8            # Speed of light in vacuum

tau = 100*1E-15                # Yin et al: 55-150fs 1E-15


def epsilon_charge(wl, d, mu = 0.2, T = 300):
    # Fn returns value of epsilon of conductor
    # wl - wavelength in nm, d thickness in m
    #to use: print epsilon_charge(500, 0.92*1E-9)
    # (1.00210490977+1.9835170885j)

    om = 2*pi*c0/wl*1.0E9

    #Drude contribution
    mu_FD = kbeV*T
    if (mu<10*mu_FD):
        sigma_Drude = 1j/(om+1j/tau)*2*(ee**2)*kb*T/(pi*(h**2))*log(2*cosh(mu/(2*kbeV*T)))
    else:
        sigma_Drude = 1j/(om+1j/tau)*(ee**2)/(pi*(h**2))*np.abs(ee*mu)

    # Interband contribution
    xi         = np.zeros((1,100001), dtype = float)
    xi[0]      = np.linspace(0, 20, 100001)
    G          = np.zeros((1,100001), dtype = float)
    Int_term   = np.zeros((1,100001), dtype = float)

    G_hom2 = sinh(heV*om/2/(kbeV*T))/(cosh(mu/(kbeV*T))+cosh(heV*om/2/(kbeV*T)))

    for i in range(len(xi[0])):
        G[0][i] = np.sinh(xi[0][i]/(kbeV*T)) /(np.cosh(mu/(kbeV*T)) + np.cosh(xi[0][i]/(kbeV*T)))
        if np.isnan(G[0][i]):
            G[0][i] = 1
        else:
            G[0][i] = G[0][i]

        Int_term[0][i] = (G[0][i]-G_hom2)/((heV*om)**2-4*xi[0][i]**2)

    #Integration = np.trapz(xi[0],Int_term[0])
    Integration = np.trapz(Int_term[0],xi[0])
    Second_Inter = 4.0*heV*om/pi*Integration
    sigma_Inter = (ee**2)/(4*h)*(G_hom2+1j*Second_Inter)
    sigma = sigma_Drude + sigma_Inter

    epsilon = 1+1j*sigma/(om*ep0*d)

    return epsilon