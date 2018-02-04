import numpy as np
import sys
from math import *
import cmath
from TMM_aniso import get_A_B
from Bi2Se3_drude import drude_O_eps, drude_E_eps
from Bi2Se3_bulk import get_eps_Bi2Se3_bulk_Yin, get_eps_Bi2Se3_bulk_Wolf
from ZnSe import get_eps_ZnSe_Marple
import pylab as pl

# CONSTANTS
#-------------------------------
c0 = 2.99792458e+8
ep0 = 8.854187817e-12
mu0 = 4*pi*1e-7
#-------------------------------

# NANO-STRUCTURE
#-------------------------------
Material_name = "Bi2Se3-ZnSe"
N_layers = 18
N_periods = 4

# Permeability
mu = np.ones((1, N_layers), dtype=float)

# Thickness of layers for TMM
d_TMM = np.zeros((1, 4*N_periods+2), dtype = float)
d_air = 1E12
d_cond = 0.92*1E-9
d_bulk = 100*1E-9
d_dielectric = 100*1E-9
aa_ = [d_dielectric, d_cond, d_bulk, d_cond]
d_TMM[0] = d_air
aa_ = np.tile(aa_, 4)
for i in range(len(d_TMM[0])-2):
    d_TMM[0][i+1] = aa_[i]

# Thickness of layers for EMT
d_EMT = np.array([ [d_air, (d_cond+d_bulk+d_dielectric)*N_periods, d_air]  ])
ff = (d_cond+d_bulk)/(d_cond+d_bulk+d_dielectric)
ffm = d_bulk/(d_cond+d_bulk+d_dielectric)
ffd = d_dielectric/(d_cond+d_bulk+d_dielectric)

# Setup wavelengths
wl = np.zeros((1,500), dtype=float)
wl[0] = np.linspace(500,1000, 500)
#wl = np.array([[400]])

# Select angles of incidence
#theta_i= np.zeros((1,180), dtype=float)
#theta_i[0] = np.linspace(0,89.9,180) # or:
theta_i= np.array([[0]])

# Setup dielectric fn for each layer at different wavelengths
eps_air = np.array([[1.0]])

eps_dO = np.zeros((1,len(wl[0])), dtype=complex)
eps_dE = np.zeros((1,len(wl[0])), dtype=complex)

eps_condO = np.zeros((1,len(wl[0])), dtype=complex)
eps_condE = np.zeros((1,len(wl[0])), dtype=complex)

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
        eps_dO[0][j] = complex(get_eps_ZnSe_Marple(wl[0][j])[0],get_eps_ZnSe_Marple(wl[0][j])[1])
        eps_dE[0][j] = complex(get_eps_ZnSe_Marple(wl[0][j])[0],get_eps_ZnSe_Marple(wl[0][j])[1])

        eps_condO[0][j] = complex(drude_O_eps(wl[0][j])[0],drude_O_eps(wl[0][j])[1])
        eps_condE[0][j] = complex(drude_E_eps(wl[0][j])[0],drude_E_eps(wl[0][j])[1])

        eps_bulkO[0][j] = get_eps_Bi2Se3_bulk_Wolf(wl[0][j])
        eps_bulkE[0][j] = get_eps_Bi2Se3_bulk_Wolf(wl[0][j])

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
        #Tr_EMT_st = np.real(kz_end/eps_EMT_O_st[-1][0])/np.real(kz_air/eps_EMT_O_st[0][0])*pow(np.abs(ApEMT_st[0][-1]/ApEMT_st[0][0]), 2)
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
print "Material:", Material_name
print "Overal number of layers:", N_layers
print "Period number:", N_periods
#print "TMM Structure:", d_TMM

#print "Rp:", Rp
#print "Tp:", Tr
#print "Ab:", Ab

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
def plot_Rp_Tp(ax, ay, xaxis, R, T):
    ax.plot(xaxis[:], Ti[:], label= 'Transmittance', color = 'r')
    #ax.plot(xaxis[:], Tst[:], label= 'T "Standard"', color = 'r', linestyle=":")
    #ax.plot(xaxis[:], Ti[:], color='r')
    ax.yaxis.tick_left()
    ax.set_ylabel('Transmittance, $T$', color = 'r')
    ax.tick_params('y',colors='r')

    ay = ax.twinx()
    ay.plot(xaxis[:], Ri[:], label= 'Reflectance', color = 'b')
    #ay.plot(xaxis[:], Rst[:], label= 'R "Standard"', color = 'b', linestyle=":")
    #ay.plot(xaxis[:], R[:], color='b')
    ay.set_ylabel('Reflectance, $R$', color = 'b')
    ay.tick_params('y',colors='b')

    #ax.set_xlabel('Incidence angle theta, 'r'$\theta$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$')

    ax.legend(loc=2, fancybox=True)
    ay.legend(loc=1, fancybox=True)
    #ax.legend(bbox_to_anchor=(0, 0.99, 1, 0), loc=2)
    #ay.legend(bbox_to_anchor=(0, 0.9, 1, 0),  loc=2)

    ymajor_ticks = np.arange(0, 1.01, 0.2)
    yminor_ticks = np.arange(0, 1.02, 0.02)

    # 1 angle - many wl:
    xmajor_ticks = np.arange(500, 1100, 100)
    xminor_ticks = np.arange(500, 1010, 10)

    # 1 wl - many angles:
    #xmajor_ticks = np.arange(0, 100, 10)
    #xminor_ticks = np.arange(0, 91, 1)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)
    ay.set_yticks(ymajor_ticks)
    ay.set_yticks(yminor_ticks, minor = True)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)
    #pl.text(0.5, 0.85,"$\lambda$ = 400 nm", horizontalalignment = 'center', verticalalignment = 'center',
    #        transform = ax.transAxes)
    pl.text(0.5, 0.85,"$\ theta$ = 0", horizontalalignment = 'center', verticalalignment = 'center',
            transform = ax.transAxes)

def doFigure_RTA(xaxis, R, x):
    fig = pl.figure()
    axR = fig.add_subplot(111)
    ayR = fig.add_subplot(111)
    plot_Rp_Tp(axR, ayR, xaxis, R, x)
    fig.tight_layout()
    #fig.savefig('../plots/plot_R_T_EMT_i1_vs_st(0 deg).png', dpi=500)
    return

R =[]
T =[]
A =[]

#EMT: fixed angle - many wl:
for i in range(len(RpEMT_i1[0])):
    R.append(RpEMT_i1[0][i])
    T.append(TrEMT_i1[0][i])
    A.append(AbEMT_i1[0][i])

#EMT: fixed wl - many angles:
#for i in range(len(RpEMT_i1)):
#    R.append(RpEMT_i1[i][0])
#    T.append(TrEMT_i1[i][0])
#    A.append(AbEMT_i1[i][0])
#------
#TMM: fixed angle - many wl:
#for i in range(len(Rp[0])):
#    R.append(Rp[0][i])
#    T.append(Tr[0][i])
#    A.append(Ab[0][i])

#TMM: fixed wl - many angles:
#for i in range(len(Rp)):
#    R.append(Rp[i][0])
#    T.append(Tr[i][0])
#    A.append(Ab[i][0])

# Select theta_i[0] or wl[0]:
#doFigure_RTA(wl[0], R, T)
#pl.show()

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
    xmajor_ticks = np.arange(500, 1100, 100)
    xminor_ticks = np.arange(500, 1010, 10)
    ymajor_ticks = np.arange(-1, 5.5, 1)
    yminor_ticks = np.arange(-1, 5.1, 0.1)

    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)
    ax.legend(loc='best',fancybox=True)
    pl.text(0.8, 0.9,"$\ theta$ = 0", horizontalalignment='center', verticalalignment = 'center',
        transform = ax.transAxes)

def doFigure_Eps(wl, epsER, epsOR, epsEI, epsOI):
    fig = pl.figure()
    axEps = fig.add_subplot(111)
    plot_Eps(axEps, wl, epsER, epsOR, epsEI, epsOI)
    fig.suptitle('Dielectric function calculated using improved EMT', fontsize=10)
    #fig.tight_layout()
    fig.savefig('../plots/plot_Eps_EMT_improved1.png', dpi=500)
    return

#doFigure_Eps(wl[0], eps_EMTE_i1[0].real, eps_EMTO_i1[0].real, eps_EMTE_i1[0].imag, eps_EMTO_i1[0].imag)
#pl.show()
