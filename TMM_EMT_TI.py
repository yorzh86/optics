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

# AUXILIARY VARIABLES
#-------------------------------
AA = []
BB = []
CC = []
#-------------------------------

# NANO-STRUCTURE
#-------------------------------
Material_name = "Ag-TiO2"
N_layers = 18
N_periods = 8

# Permeability and thickness of layers
mu = np.ones((1, N_layers), dtype=float)
d_TMM = np.ones((1, N_layers), dtype=float)
d_TMM[0] = 0.03
d_TMM[0][0] = d_TMM[0][-1] = 1E6

# Setup wavelengths
#wl = np.array([[0.4, 0.5, 0.6, 0.7, 0.8]])
wl = np.array([[0.4]])

# Setup dielectric fn for each layer at different wavelengths
eps_air = np.array([[1.0]])
#eps_mO = np.array([[-7.7123+0.0505j, -12.6127+0.0986j,-18.6019+0.1705j, -25.6796+0.2707j,-33.8458+0.4040j]])
#eps_mE = np.array([[-7.7123+0.0505j, -12.6127+0.0986j,-18.6019+0.1705j, -25.6796+0.2707j,-33.8458+0.4040j]])

eps_mO=np.array([[3, 3, 3, 4, 1]])
eps_mE=np.array([[2, 2, 1, 1, 1]])

#eps_dO = np.array([[9.1038+0.0782j, 7.3514,  6.7857, 6.5088, 6.3491]])
#eps_dE = np.array([[11.4552+0.5068j, 9.2018, 8.4019, 8.0158, 7.7948]])

eps_dO = np.array([[1, 2, 3, 4, 5]])
eps_dE = np.array([[2, 2, 3, 3, 4]])

# Select angles of incidence
#theta_i= np.zeros((1,180), dtype=float)
#theta_i[0] = np.linspace(0,89.9,180) # or:
theta_i= np.array([[0, 10, 20, 30]])
#-------------------------------

# TESTING
#-------------------------------
# 1. Loop through all angles of incidence
# 2. Loop through all wavelengths
for i in range(len(theta_i[0])):
    for j in range(len(wl[0])):

        #Epsilon of period (metal, dielectric)
        eps_period_O = np.array([[eps_mO[0][j], eps_dO[0][j]]])
        eps_period_E = np.array([[eps_mE[0][j], eps_dE[0][j]]])

        # Build a long array of nano-structure:
        # 1. do (period x N_periods)
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

        om = 2*pi*c0/wl[0][j]*1E6#*1E6?? 1E9?
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

print
print "Material:", Material_name
print "Tested incidence angles:", theta_i[0]
print "Tested wavelengths,[um]:", wl[0], '\n'
print "Rp:", Rp_TM
print "Tp:", Tr_TM

#print('\x1b[6;30;42m' + 'Success!' + '\x1b[0m')
#add penetration depth
#for l in range (wl.size):
#    n = sqrt(eps_[])
#    kappa = n.imag
#    pd = wl[l]/(4.0*pi*kappa)
#-------------------------------

# POST PROCESSING
#-------------------------------
# 1. Perform plotting magic
def plotRp_Tp(ax, ay, xaxis, R, T):
    ax.plot(xaxis[:], T[:], label= 'Transmittance')
    ax.plot(xaxis[:], T[:], color='r', linewidth=0.4, linestyle='-')
    ax.yaxis.tick_left()
    ax.set_ylabel('Transmittance, $T$', color = 'r')
    ax.tick_params('y',colors='r')

    ay = ax.twinx()
    ay.plot(xaxis[:], R[:], label= 'Reflectance')
    ay.plot(xaxis[:], R[:], color='darkblue', linewidth=0.3, linestyle='--')
    ay.set_ylabel('Reflectance, $R$', color = 'b')
    ay.tick_params('y',colors='b')

    ax.set_xlabel('Incidence angle theta, $theta_i$')
    ax.legend(loc=3,fancybox=True)
    ay.legend(loc=2,fancybox=True)
    #ay.legend(loc='best',fancybox=True)

    ymajor_ticks = np.arange(0, 1.01, 0.2)
    yminor_ticks = np.arange(0, 1.02, 0.02)
    xmajor_ticks = np.arange(300, 800, 100)
    xminor_ticks = np.arange(300, 725, 25)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)
    ay.set_yticks(ymajor_ticks)
    ay.set_yticks(yminor_ticks, minor = True)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)
    #ax(ay).minorticks_on()
    #ax.grid(which='both')
    #ax.grid(which='minor', linestyle = '--')
    #ax.grid(which='major', linestyle = '--')

def doFigure(xaxis, Rp_TM, Tr_TM):
    R = Rp_TM
    T = Tr_TM
    fig = pl.figure()
    axR = fig.add_subplot(111)
    ayR = fig.add_subplot(111)
    plotRp_Tp(axR, ayR, xaxis, R, T)
    fig.tight_layout()
    #pl.subplots_adjust(top=0.91, hspace=0.46)
    #fig.suptitle('Reflectance/Transmittance', fontsize=12, fontweight='bold')
    fig.savefig('plotR_T.pdf')
    return

#~ doFigure(wl[0], Rp_TM, Tr_TM)
#~ pl.show()
#-------------------------------
