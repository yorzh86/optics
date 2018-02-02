#Concerns:
#1. interpolation gives WAAAY better results than Drude model
#2. ZnSe starts only from 490 nm.
#3. High absorbtance???

#todo:
#- plot for different angles
#- re-do Thickness array creation (hard coded 4 periods (d_TMM)

import numpy as np
from math import *
import cmath
from TMM_aniso import get_A_B
from Bi2Se3_drude import drude_O_eps, drude_E_eps
from Bi2Se3_bulk import get_eps_Bi2Se3_bulk
from ZnSe import get_eps_ZnSe
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
Material_name = "Bi2Se3-ZnSe"
N_layers = 18
N_periods = 4

# Permeability and thickness of layers (nm)
mu = np.ones((1, N_layers), dtype=float)
#d_TMM = np.ones((1, N_layers), dtype=float)
#d_TMM[0] = 100*1E-9
#d_TMM[0][0] = d_TMM[0][-1] = 1E12
d_TMM = np.array([
[1E12,   0.92E-9, 100E-9, 100E-9, 0.92E-9,    0.92E-9, 100E-9, 100E-9, 0.92E-9,   \
 0.92E-9, 100E-9, 100E-9, 0.92E-9,   0.92E-9, 100E-9, 100E-9, 0.92E-9,    1E12]])

# Setup wavelengths
wl = np.array([[500,600]])

# Select angles of incidence
#theta_i= np.zeros((1,180), dtype=float)
#theta_i[0] = np.linspace(0,89.9,180) # or:
theta_i= np.array([[0,10,20,30]])

# Setup dielectric fn for each layer at different wavelengths
eps_air = np.array([[1.0]])
eps_dO = np.zeros((1,len(wl[0])), dtype=complex)
eps_dE = np.zeros((1,len(wl[0])), dtype=complex)
eps_condO = np.zeros((1,len(wl[0])), dtype=complex)
eps_condE = np.zeros((1,len(wl[0])), dtype=complex)
eps_bulkO = np.zeros((1,len(wl[0])), dtype=complex)
eps_bulkE = np.zeros((1,len(wl[0])), dtype=complex)

Rp = np.zeros((len(theta_i[0]),len(wl[0])),dtype=float)
Tr = np.zeros((len(theta_i[0]),len(wl[0])),dtype=float)
Ab = np.zeros((len(theta_i[0]),len(wl[0])),dtype=float)

#-------------------------------

# TESTING
#-------------------------------
# 1. Loop through all angles of incidence
# 2. Loop through all wavelengths
for i in range(len(theta_i[0])):
    for j in range(len(wl[0])):
        # Setup dielectric fn for each layer at different wavelengths
        eps_dO[0][j] = complex(get_eps_ZnSe(wl[0][j])[0],get_eps_ZnSe(wl[0][j])[1])
        eps_dE[0][j] = complex(get_eps_ZnSe(wl[0][j])[0],get_eps_ZnSe(wl[0][j])[1])
        eps_condO[0][j] = complex(drude_O_eps(wl[0][j])[0],drude_O_eps(wl[0][j])[1])
        eps_condE[0][j] = complex(drude_E_eps(wl[0][j])[0],drude_E_eps(wl[0][j])[1])      
        eps_bulkO[0][j] = complex(get_eps_Bi2Se3_bulk(wl[0][j])[0],get_eps_Bi2Se3_bulk(wl[0][j])[1])
        eps_bulkE[0][j] = complex(get_eps_Bi2Se3_bulk(wl[0][j])[0],get_eps_Bi2Se3_bulk(wl[0][j])[1])

        #Epsilon of period (dielectric, conduction, bulk, conduction)
        eps_period_O = np.array([
            [eps_dO[0][j], eps_condO[0][j], eps_bulkO[0][j], eps_condO[0][j] ]])
        eps_period_E = np.array([
            [eps_dE[0][j], eps_condE[0][j], eps_bulkE[0][j], eps_condE[0][j] ]])

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

        om = 2*pi*c0/wl[0][j]*1E9
        k0 = float(om/c0)
        kx = k0*sin(theta_i[0][i]*pi/180)
        kz_air = cmath.sqrt(pow(k0,2)-pow(kx,2))
        kz_end = cmath.sqrt(pow(k0,2)*eps_TMM_O[0][-1]-pow(kx,2)*eps_TMM_O[0][-1]/eps_TMM_E[0][-1])

        Ap, Bp = get_A_B(d_TMM, N_layers, wl[0][j], eps_TMM_O, eps_TMM_E, kx, mu)

        Rp_TM = pow(np.abs(Bp[0][0]/Ap[0][0]),2)
        Tr_TM = np.real(kz_end/eps_TMM_O[-1][0])/np.real(kz_air/eps_TMM_O[0][0])*pow(np.abs(Ap[0][-1]/Ap[0][0]), 2)
        Ab_TM = 1 - Rp_TM - Tr_TM

        #AA.append(Rp_TM)
        #BB.append(Tr_TM)
        #CC.append(Ab_TM)
        
        Rp[i][j] = Rp_TM
        Tr[i][j] = Tr_TM
        Ab[i][j] = Ab_TM

#Rp_TM = np.array(AA)
#Tr_TM = np.array(BB)
#Ab_TM = np.array(CC)


#~ print "Rp:", "%0.3f" % Rp_TM[i]

#~ print
#~ print "Material:", Material_name
#~ print "Tested incidence angles:", theta_i[0]
#~ print "Tested wavelengths,[um]:", wl[0], '\n'

#~ print "Rp:", Rp
#~ print
#~ print "Tp:", Tr
#~ print
#~ print "Ab:", Ab

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
    ay.plot(xaxis[:], R[:], color='b', linewidth=0.3, linestyle='-.')
    ay.set_ylabel('Reflectance, $R$', color = 'b')
    ay.tick_params('y',colors='b')

    ax.set_xlabel('Incidence angle theta, 'r'$\theta$')
    ax.legend(loc=3,fancybox=True)
    ay.legend(loc=2,fancybox=True)

    ymajor_ticks = np.arange(0, 1.01, 0.2)
    yminor_ticks = np.arange(0, 1.02, 0.02)

    # 1 angle - many wl:
    #xmajor_ticks = np.arange(300, 800, 100)
    #xminor_ticks = np.arange(300, 725, 25)

    # 1 wl - many angles:
    xmajor_ticks = np.arange(0, 100, 10)
    xminor_ticks = np.arange(0, 91, 1)

    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)
    ay.set_yticks(ymajor_ticks)
    ay.set_yticks(yminor_ticks, minor = True)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)
    pl.text(0.7, 0.85,"$\lambda$ = 400 nm", horizontalalignment = 'center', verticalalignment = 'center',
            transform = ax.transAxes)
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
    fig.savefig('plotR_T.pdf')
    return

# Select theta_i[0] or wl[0]:
#doFigure(theta_i[0], Rp_TM, Tr_TM)
#pl.show()
#-------------------------------
