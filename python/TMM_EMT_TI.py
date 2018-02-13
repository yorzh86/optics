from __future__ import division
import numpy as np
import sys
from math import *
import cmath
from TMM_aniso import get_A_B

from Bi2Se3_drude import drude_E_eps
from Bi2Se3_bulk import  bulk_Wolf
from Bi2Se3_properties import *

#from Bi2Te3_drude import drude_E_eps
#from Bi2Te3_bulk import bulk_Wolf
#from Bi2Se3_properties import *

from ZnSe import eps_ZnSe_Marple
from sigma_epsilon import eps_conductor
import pylab as pl

# CONSTANTS
#-------------------------------
c0 = 2.99792458e+8
ep0 = 8.854187817e-12
mu0 = 4*pi*1e-7
#-------------------------------

# NANO-STRUCTURE
#-------------------------------
N_periods = 10
N_layers = 2+N_periods*4
Plot_resolution = 50

# Permeability
mu = np.ones((1, N_layers), dtype=float)

# Thickness of layers for TMM
d_TMM = np.zeros((1, 4*N_periods+2), dtype = float)
d_air = 0
d_cond = d_conduct()
d_bulk = (50-d_cond*2)*1E-9
d_dielectric = 30*1E-9

aa_ = [d_dielectric, d_cond, d_bulk, d_cond]
d_TMM[0] = d_air
aa_ = np.tile(aa_, N_periods)
for i in range(len(d_TMM[0])-2):
    d_TMM[0][i+1] = aa_[i]

# Thickness of layers for EMT
d_EMT = np.array([ [d_air, (d_cond+d_bulk+d_dielectric+d_cond)*N_periods, d_air]  ])

ff = (d_cond*2+d_bulk)/(d_cond*2+d_bulk+d_dielectric)
ffm = d_bulk/(d_cond*2+d_bulk+d_dielectric)
ffd = d_dielectric/(d_cond*2+d_bulk+d_dielectric)

# Setup wavelengths
wl = np.zeros((1,Plot_resolution), dtype=float)
wl[0] = np.linspace(500, 3500, Plot_resolution)
#wl[0] = np.linspace(500, 20000, Plot_resolution)

# Select angles of incidence
#theta_i= np.zeros((1,180), dtype=float)
#theta_i[0] = np.linspace(0,89.9,180) # or:
theta_i= np.array([[0]])

# Setup dielectric fn for each layer at different wavelengths
eps_air = np.array([[1.0]])

eps_dE = np.zeros((1,len(wl[0])), dtype=complex)
eps_dO = np.zeros((1,len(wl[0])), dtype=complex)

eps_condE = np.zeros((1,len(wl[0])), dtype=complex)
eps_condO = np.zeros((1,len(wl[0])), dtype=complex)

eps_bulkO = np.zeros((1,len(wl[0])), dtype=complex)
eps_bulkE = np.zeros((1,len(wl[0])), dtype=complex)

eps_EMTO_st = np.zeros((1,len(wl[0])), dtype=complex)
eps_EMTE_st = np.zeros((1,len(wl[0])), dtype=complex)

eps_EMTO_i1 = np.zeros((1,len(wl[0])), dtype=complex)
eps_EMTE_i1 = np.zeros((1,len(wl[0])), dtype=complex)

Rp = np.zeros((len(theta_i[0]),len(wl[0])),dtype=float)
Tr = np.zeros((len(theta_i[0]),len(wl[0])),dtype=float)
Ab = np.zeros((len(theta_i[0]),len(wl[0])),dtype=float)

RpEMT_st = np.zeros((len(theta_i[0]),len(wl[0])),dtype=float)
TrEMT_st = np.zeros((len(theta_i[0]),len(wl[0])),dtype=float)
AbEMT_st = np.zeros((len(theta_i[0]),len(wl[0])),dtype=float)

RpEMT_i1 = np.zeros((len(theta_i[0]),len(wl[0])),dtype=float)
TrEMT_i1 = np.zeros((len(theta_i[0]),len(wl[0])),dtype=float)
AbEMT_i1 = np.zeros((len(theta_i[0]),len(wl[0])),dtype=float)

#-------------------------------

# TESTING
#-------------------------------
# 1. Loop through all angles of incidence
# 2. Loop through all wavelengths
for i in range(len(theta_i[0])):
    for j in range(len(wl[0])):
        # Setup dielectric fn for each layer at different wavelengths
        eps_dO[0][j] = complex(eps_ZnSe_Marple(wl[0][j])[0],eps_ZnSe_Marple(wl[0][j])[1])
        eps_dE[0][j] = complex(eps_ZnSe_Marple(wl[0][j])[0],eps_ZnSe_Marple(wl[0][j])[1])


        eps_bulkO[0][j] = bulk_Wolf(wl[0][j])
        eps_bulkE[0][j] = bulk_Wolf(wl[0][j])

        eps_condO[0][j] = eps_conductor(wl[0][j], d_cond, tau(), muf()) + eps_bulkO[0][j]
        eps_condE[0][j] = complex(drude_E_eps(wl[0][j])[0],drude_E_eps(wl[0][j])[1])

        eps_EMTO_st[0][j] = ff*eps_bulkO[0][j]+(1.0-ff)*eps_dO[0][j]
        eps_EMTE_st[0][j] = eps_bulkE[0][j]*eps_dE[0][j]/(ff*eps_dE[0][j]+(1-ff)*eps_bulkE[0][j])

        eps_EMTO_i1[0][j] = ffm*eps_bulkO[0][j]+ ffd*eps_dO[0][j] + (1-ffm-ffd)*eps_condO[0][j]
        eps_EMTE_i1[0][j] = 1.0/(ffm/eps_bulkO[0][j] + ffd/eps_dO[0][j] + (1-ffm-ffd)/eps_condE[0][j])

        eps_EMT_O_st = np.array([ [eps_air, eps_EMTO_st[0][j], eps_air] ])
        eps_EMT_E_st = np.array([ [eps_air, eps_EMTE_st[0][j], eps_air] ])
        eps_EMT_O_i1 = np.array([ [eps_air, eps_EMTO_i1[0][j], eps_air] ])
        eps_EMT_E_i1 = np.array([ [eps_air, eps_EMTE_i1[0][j], eps_air] ])

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

        ApEMT_st, BpEMT_st = get_A_B(d_EMT, 3, wl[0][j], eps_EMT_O_st, eps_EMT_E_st, kx, mu)
        ApEMT_i1, BpEMT_i1 = get_A_B(d_EMT, 3, wl[0][j], eps_EMT_O_i1, eps_EMT_E_i1, kx, mu)

        Rp_TM = pow(np.abs(Bp[0][0]/Ap[0][0]),2)
        Tr_TM = np.real(kz_end/eps_TMM_O[-1][0])/np.real(kz_air/eps_TMM_O[0][0])*pow(np.abs(Ap[0][-1]/Ap[0][0]), 2)
        Ab_TM = 1 - Rp_TM - Tr_TM

        Rp_EMT_st = pow(np.abs(BpEMT_st[0][0]/ApEMT_st[0][0]),2)
        Tr_EMT_st = pow(np.abs(ApEMT_st[0][-1]/ApEMT_st[0][0]),2)
        Ab_EMT_st = 1 - Rp_EMT_st - Tr_EMT_st

        Rp_EMT_i1 = pow(np.abs(BpEMT_i1[0][0]/ApEMT_i1[0][0]),2)
        Tr_EMT_i1 = pow(np.abs(ApEMT_i1[0][-1]/ApEMT_i1[0][0]),2)
        Ab_EMT_i1 = 1 - Rp_EMT_i1 - Tr_EMT_i1

        Rp[i][j] = Rp_TM
        Tr[i][j] = Tr_TM
        Ab[i][j] = Ab_TM

        RpEMT_st[i][j] = Rp_EMT_st
        TrEMT_st[i][j] = Tr_EMT_st
        AbEMT_st[i][j] = Ab_EMT_st

        RpEMT_i1[i][j] = Rp_EMT_i1
        TrEMT_i1[i][j] = Tr_EMT_i1
        AbEMT_i1[i][j] = Ab_EMT_i1


print
print "Material:", material_name()
print "Overal number of layers:", N_layers
print "Period number:", N_periods


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
def plot_Rp_Tp(ax, xaxis, TMM, EMTi, EMTst):
    ax.plot(xaxis[:], TMM[:], label= 'Transmittance TMM', color = 'r', linewidth = 0.4)
    ax.plot(xaxis[:], EMTi[:], label= 'Transmittance EMTi', color = 'b', linewidth = 0.4)
    ax.plot(xaxis[:], EMTst[:], label= 'Transmittance EMTst', color = 'g', linestyle=':')#linewidth = 0.4)
    ax.yaxis.tick_left()
    ax.set_ylabel('Transmittance, $T$')#, color = 'r')

    #ax.set_xlabel('Incidence angle theta, 'r'$\theta$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$')
    ax.legend(loc=2, fancybox=True)

    ymajor_ticks = np.arange(0, 1.01, 0.2)
    yminor_ticks = np.arange(0, 1.02, 0.02)

    # 1 angle - many wl:
    xmajor_ticks = np.arange(500, 4000, 500)
    xminor_ticks = np.arange(500, 3550, 50)
    #xmajor_ticks = np.arange(500, 24000, 4000)
    #xminor_ticks = np.arange(500, 21000, 1000)

    # 1 wl - many angles:
    #xmajor_ticks = np.arange(0, 100, 10)
    #xminor_ticks = np.arange(0, 91, 1)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)
    #pl.text(0.5, 0.85,"$\lambda$ = 400 nm", horizontalalignment = 'center', verticalalignment = 'center',
    #        transform = ax.transAxes)
    pl.text(0.2, 0.7,"$\ theta$ = 0", horizontalalignment = 'center', verticalalignment = 'center',
            transform = ax.transAxes)

def doFigure_RTA(xaxis, T, Ti, Tst):
    fig = pl.figure()
    axR = fig.add_subplot(111)
    ayR = fig.add_subplot(111)
    plot_Rp_Tp(axR,xaxis, T, Ti, Tst)
    #fig.tight_layout()
    fig.suptitle('Transmittance of Bi2Se3 via TMM and EMT(0.5-3.5um)', fontsize=10)
    fig.savefig('../plots/updateFeb11/Bi2Se3/TMM,EMTi,EMTst_T_T_0.5-3.5um,bulk=8.16nm.png', dpi=500)
    return

R =[]
T =[]
A =[]

Rst =[]
Tst =[]
Ast =[]

Ri =[]
Ti =[]
Ai =[]

#EMT: fixed angle - many wl:
for i in range(len(RpEMT_st[0])):
    Rst.append(RpEMT_st[0][i])
    Tst.append(TrEMT_st[0][i])
    Ast.append(AbEMT_st[0][i])

for i in range(len(RpEMT_i1[0])):
    Ri.append(RpEMT_i1[0][i])
    Ti.append(TrEMT_i1[0][i])
    Ai.append(AbEMT_i1[0][i])

#EMT: fixed wl - many angles:
#for i in range(len(RpEMT_i1)):
#    R.append(RpEMT_i1[i][0])
#    T.append(TrEMT_i1[i][0])
#    A.append(AbEMT_i1[i][0])
#------
#TMM: fixed angle - many wl:
for i in range(len(Rp[0])):
    R.append(Rp[0][i])
    T.append(Tr[0][i])
    A.append(Ab[0][i])

#TMM: fixed wl - many angles:
#for i in range(len(Rp)):
#    R.append(Rp[i][0])
#    T.append(Tr[i][0])
#    A.append(Ab[i][0])

# Select theta_i[0] or wl[0]:
#doFigure_RTA(wl[0], T, Ti, Tst)
#pl.show()

# ============== epsilon======================
# Attention - works only for 1 theta and range of wavelengths
def plot_Eps(ax, wl, epsER, epsOR, epsEI, epsOI):
    ax.plot(wl, epsER, label= "Extraordinary, real", color='r')
    ax.plot(wl, epsOR, label= "Ordinary, real", color='b')

    ax.plot(wl, epsEI, color='r', linestyle=':',  label='Extraordinary, imaginary')
    ax.plot(wl, epsOI, color='b', linestyle=':', label='Ordinary, imaginary')
    ax.plot(wl, wl*0, linestyle = '-.', color='grey', linewidth = 0.3)
    ax.set_ylabel('Epsilon, 'r'$\epsilon$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')
    #ax.minorticks_on()

    xmajor_ticks = np.arange(500, 4000, 500)
    xminor_ticks = np.arange(500, 3550, 50)
    #xmajor_ticks = np.arange(500, 24000, 4000)
    #minor_ticks = np.arange(500, 21000, 1000)
    ymajor_ticks = np.arange(-5, 25, 5)
    yminor_ticks = np.arange(-5, 25, 0.5)

    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)

    ax.legend(loc=2, fancybox=True)
    pl.text(0.1, 0.65,"$\ theta$ = 0", horizontalalignment='center', verticalalignment = 'center',
        transform = ax.transAxes)

def doFigure_Eps(wl, epsER, epsOR, epsEI, epsOI):
    fig = pl.figure()
    axEps = fig.add_subplot(111)
    plot_Eps(axEps, wl, epsER, epsOR, epsEI, epsOI)
    fig.suptitle('Dielectric function of Bi2Se3 using improved EMTi.', fontsize=10)
    #fig.tight_layout()
    #fig.savefig('../plots/updateFeb11/Bi2Se3/Eps_EMTi(0.5-3.5um), bulk=48.16.png', dpi=500)
    return

doFigure_Eps(wl[0], eps_EMTE_i1[0].real, eps_EMTO_i1[0].real, eps_EMTE_i1[0].imag, eps_EMTO_i1[0].imag)
pl.show()