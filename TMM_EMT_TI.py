# -*- coding: utf-8 -*-
"""
@author: yorzh
"""
import numpy as np
from math import *
import cmath
from TMM_aniso import get_A_B
import pylab as pl

# CONSTANTS
#-------------------------------
c0 = 2.99792458e+8
ep0 = 8.854187817e-12
mu0 = 4*pi*1e-7
#-------------------------------

# NANO-STRUCTURE
#-------------------------------
#~ d_ti = 100.92
#~ d_diel = 100
#~ ff = d_ti/(d_ti+d_diel)
N_layers = 18
N_periods = 8

mu = np.ones((1, N_layers), dtype=float)
d_TMM = np.ones((1, N_layers), dtype=float)
d_TMM[0] = 30.0
d_TMM[0][0] = d_TMM[0][-1] = 10E6

# select wavelengths
wl = np.array([[400, 500, 600]])

# Setup dielectric fn of material for different wavelength
eps_air = np.array([[1.0]])
eps_mO = np.array([4.2, 5.2, 6.2])
eps_dO = np.array([2.2, 1.2, 0.2])
eps_dE = np.array([2.2, 1.2, 0.2])
eps_mE = np.array([2.2, 1.2, 0.2])

# select angles of incidence
#theta_i= np.zeros((1,180), dtype=float)
#theta_i[0] = np.linspace(0,89.9,180) # or:
theta_i= np.array([[0]])

#-------------------------------

# TESTING
#-------------------------------
# Set up auxiliary lists (for easy appending)
AA =[]
BB = []
CC = []

for i in range(len(theta_i[0])):
    for j in range(len(wl[0])):
        #Epsilon of period (metal, dielectric)
        # move it to wavelength loop
        eps_period_O = np.array([[eps_mO[j], eps_dO[j]]])
        eps_period_E = np.array([[eps_mE[j], eps_dE[j]]])
        # Build a long array of total structure:
        # 1. do period x N_periods
        # 2. do (air) + (period x N_periods)
        # 3. do (air + period x N_periods)+ air
        # 4. change shape to match the rest of the code
        eps_period_O = np.tile(eps_period_O, N_periods)
        eps_TMM_O = np.append(eps_air, eps_period_O)
        eps_TMM_O = np.append(eps_TMM_O, eps_air)
        eps_TMM_O = np.array([eps_TMM_O])
        
        eps_period_E = np.tile(eps_period_E, N_periods)
        eps_TMM_E = np.append(eps_air, eps_period_E)
        eps_TMM_E = np.append(eps_TMM_E, eps_air)
        eps_TMM_E = np.array([eps_TMM_E])

        om = 2*pi*c0/wl[0][j]*1E9
        k0 = float(om/c0)        
        kx = k0*sin(theta_i[0][i]*pi/180)
        kz_air = sqrt(pow(k0,2)-pow(kx,2))
        kz_end = sqrt(pow(k0,2)*eps_TMM_O[0][-1]-pow(kx,2)*eps_TMM_O[0][-1]/eps_TMM_E[0][-1])
    
        Ap, Bp = get_A_B(d_TMM, N_layers, wl[0][j], eps_TMM_O, eps_TMM_E, kx, mu)

        Rp_TM = pow(np.abs(Bp[0][0]/Ap[0][0]),2)
        Tr_TM = np.real(kz_end/eps_TMM_O[-1][0])/np.real(kz_air/eps_TMM_O[0][0])*pow(np.abs(Ap[0][-1]/Ap[0][0]), 2)
        Ab_TM = 1 - Rp_TM - Tr_TM
        AA.append(Rp_TM)
        BB.append(Tr_TM)
        CC.append(Ab_TM)

Rp_TM = np.array(AA)
Tr_TM = np.array(BB)
Ab_TM = np.array(CC)

#~ print "Rp:", "%0.3f" % Rp_TM[i]
print "Rp:", Rp_TM, '\n'
print "Tp:", Tr_TM, '\n'
print "Ab:", Ab_TM, '\n'

#add penetration depth
#for l in range (wl.size):
#    n = sqrt(eps_[])
#    kappa = n.imag
#    pd = wl[l]/(4.0*pi*kappa)
#-------------------------------

# POST PROCESSING
#-------------------------------
def plotRp_Tp(ax, ay,theta_i, R, T):
    ax.plot(theta_i[:], R[:], label= 'Reflectance')
    ax.plot(theta_i[:], R[:], color='red', linewidth=0.5, linestyle='-')
    ax.yaxis.tick_left()
    ax.set_ylabel('Reflectance, $R$', color = 'r')
    ax.tick_params('y',colors='r')

    ay = ax.twinx()
    ay.plot(theta_i[:], T[:], label= 'Transmittance')
    ay.plot(theta_i[:], T[:], color='darkblue', linewidth=0.5, linestyle='-')
    ay.set_ylabel('Transmittance, $T$', color = 'b')
    ay.tick_params('y',colors='b')
    ax.minorticks_on()
    ay.minorticks_on()
    ax.set_xlabel('Incidence angle theta, $theta_i$')
    ax.legend(loc=3,fancybox=True)
    ay.legend(loc=2,fancybox=True)
    #ay.legend(loc='best',fancybox=True)

def doFigure(theta_i, Rp_TM, Tr_TM):
    theta = theta_i
    R = Rp_TM
    T = Tr_TM
    fig = pl.figure()
    axR = fig.add_subplot(111)
    ayR = fig.add_subplot(111)
    plotRp_Tp(axR, ayR, theta,R, T)
    fig.tight_layout()
    #pl.subplots_adjust(top=0.91, hspace=0.46)
    #fig.suptitle('Reflectance/Transmittance', fontsize=12, fontweight='bold')
    fig.savefig('plotR_T.pdf')
    return

#doFigure(theta_i[0], Rp_TM, Tr_TM)
#pl.show()
#-------------------------------
