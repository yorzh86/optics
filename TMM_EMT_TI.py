# -*- coding: utf-8 -*-
"""
@author: yorzh
"""
import numpy as np
from math import *
import cmath
from TMM_aniso import get_A_B
import pylab as pl

wl = 12
# wl in eV; eV to nm -> 1eV = 1.2398E+3 nm
c0 = 2.99792458e+8
ep0 = 8.854187817e-12
mu0 = 4*pi*1e-7
mu = np.array([[1,1,1,1]])
N_layers = 4
om = 2*pi*c0/wl*1E6
k0 = float(om/c0)

theta_i= np.array([[10,20,30,40,50,60,70,80]])
        
d_ti = 100.92
# Bi2Se3 conducting 0.92 + 100 nm bulk thickness topological insulator nanometers
d_diel = 100
# thickness dielectric nanometers
ff = d_ti/(d_ti+d_diel)
# filling ratio

eps1 = 1.0
eps3 = 1.0

eps_dO = 2.2
eps_mO = 4.2
eps_dE = 2.2
eps_mE = 4.2

eps_TMM_O = np.array([[eps1, eps_dO, eps_mO, eps3]])
eps_TMM_E = np.array([[eps1, eps_dE, eps_mE, eps3]])
# structure described by _O and _E epsilons
d_TMM = np.array([[1, 1, 1, 1]])
AA =[]
BB = []
CC = []
# auxiliary lists (for easy appending)

#eps_diel_O = np.array([6.05, 6.00, 5.95, 5.95, 5.90, 5.90])
#eps_diel_E = np.array([6.05, 6.00, 5.95, 5.95, 5.90, 5.90])
## Dielectric. Isotropic
#eps_bulk_O = np.array([22.0+5j, 25.0+7j, 27.0+7j, 28.0+11j, 30.0+18j, 20.0+20j])
#eps_bulk_E = np.array([22.0+5j, 25.0+7j, 27.0+7j, 28.0+11j, 30.0+18j, 20.0+20j])
## Bulk isotropic (bulk - interband Yin J)
#eps_cond_O = np.array([24.47+27.0j, -11.0+8.6j, -3.73+2.5j, -0.71+1.8j, 0.11+1.7j, 0.96+1.6j])
#eps_cond_E = np.array([24.47+27.0j, -11.0+8.6j, -3.73+2.5j, -0.71+1.8j, 0.11+1.7j, 0.96+1.6j])
# Conduction surface anisotropic. Needs correction! (black/red lines)
#d_TMM = np.array([1000, 100, 0.92, 100, 0.92, 1000])

for i in range(len(theta_i[0])):
    kx = k0*sin(theta_i[0][i]*pi/180)
    kz_air = sqrt(pow(k0,2)-pow(kx,2))
    kz_end = sqrt(pow(k0,2)*eps_TMM_O[0][-1]-pow(kx,2)*eps_TMM_O[0][-1]/eps_TMM_E[0][-1])

    Ap, Bp = get_A_B(d_TMM, N_layers, wl, eps_TMM_O, eps_TMM_E, kx, mu)

    Rp_TM = pow(np.abs(Bp[0][0]/Ap[0][0]),2)
    Tr_TM = np.real(kz_end/eps_TMM_O[-1][0])/np.real(kz_air/eps_TMM_O[0][0])*pow(np.abs(Ap[0][-1]/Ap[0][0]), 2)
    Ab_TM = 1 - Rp_TM - Tr_TM
    AA.append(Rp_TM)
    BB.append(Tr_TM)
    CC.append(Ab_TM)

Rp_TM = np.array(AA)
Tr_TM = np.array(BB)
Ab_TM = np.array(CC)

#~ print "Rp:", Rp_TM
#~ print "Tp:", Tr_TM
#~ print "Ab:", Ab_TM

#add penetration depth
#for l in range (wl.size):
#    n = sqrt(eps_[])
#    kappa = n.imag
#    pd = wl[l]/(4.0*pi*kappa)

#--------------------------
# here goes plotting magic:
#--------------------------
def plotRp_Tp(ax, ay,theta_i, R, T):
    ax.plot(theta_i[:], R[:],'r.', label= 'Reflectance', markersize = 15)
    ax.plot(theta_i[:], R[:], color='darkred', linewidth=0.8, linestyle='-')
    ax.yaxis.tick_left()
    ax.set_ylabel('Reflectance, $R$', color = 'r')
    ax.tick_params('y',colors='r')
    
    ay = ax.twinx()
    ay.plot(theta_i[:], T[:],'b*', label= 'Transmittance', markersize = 15)
    ay.plot(theta_i[:], T[:], color='darkblue', linewidth=0.8, linestyle='-')
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

doFigure(theta_i[0], Rp_TM, Tr_TM)
pl.show()
