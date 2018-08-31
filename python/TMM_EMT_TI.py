from __future__ import division
import numpy as np
import sys
from math import *
import cmath
from TMM_aniso import get_A_B
import postprocess

from Bi2Se3_lorentz import lorentz_E_eps
from Bi2Se3_bulk import  bulk_Wolf
from Bi2Se3_properties import *

#from Bi2Te3_lorentz import lorentz_E_eps
#from Bi2Te3_bulk import bulk_Wolf
#from Bi2Te3_properties import *

#from ZnSe import eps_ZnSe_Marple
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
Wavelength_resolution = 10
Angle_resolution = 18

# Permeability
mu = np.ones((1, N_layers), dtype=float)

# Thickness of layers for TMM
d_TMM = np.zeros((1, 4*N_periods+2), dtype = float)
d_air = 0
d_cond = d_conduct()
d_bulk = (100*1E-9-d_cond*2)
d_dielectric = 100*1E-9

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
wl = np.zeros((1,Wavelength_resolution), dtype=float)
#wl[0] = np.linspace(500, 3500, Wavelength_resolution)
wl[0] = np.linspace(500, 20000, Wavelength_resolution)
#wl[0] = np.linspace(500, 13360, Wavelength_resolution)

# Select angles of incidence
theta_i= np.zeros((1,Angle_resolution), dtype=float)
theta_i[0] = np.linspace(0,89.9,Angle_resolution) # or:
#theta_i= np.array([[0]])


kx = np.zeros((len(theta_i[0]),len(wl[0])), dtype=complex)
kz_air = np.zeros((len(theta_i[0]),len(wl[0])), dtype=complex)
kz_end = np.zeros((len(theta_i[0]),len(wl[0])), dtype=complex)

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

# Penetration depth array
Pd = np.zeros((2, len(wl[0])), dtype=float)

#-------------------------------

# TESTING
#-------------------------------
# 1. Loop through all angles of incidence
# 2. Loop through all wavelengths
for i in range(len(theta_i[0])):
    for j in range(len(wl[0])):
        # Setup dielectric fn for each layer at different wavelengths
        eps_dO[0][j] = complex(eps_dielectric()[0],eps_dielectric()[1])
        eps_dE[0][j] = complex(eps_dielectric()[0],eps_dielectric()[1])

        eps_bulkO[0][j] = bulk_Wolf(wl[0][j])
        eps_bulkE[0][j] = bulk_Wolf(wl[0][j])

        eps_condO[0][j] = eps_conductor(wl[0][j], d_cond, tau(), muf()) + eps_bulkO[0][j]
        eps_condE[0][j] = complex(lorentz_E_eps(wl[0][j])[0],lorentz_E_eps(wl[0][j])[1])

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

        om = 2*pi*c0/wl[0]*1E9 #wl[0][j]
        k0 = om/c0
        kx[i][j] = k0[j]*sin(theta_i[0][i]*pi/180)
        kz_air[i][j] = cmath.sqrt(pow(k0[j],2)-pow(kx[i][j],2))
        kz_end[i][j] = cmath.sqrt(pow(k0[j],2)*eps_TMM_O[0][-1]-pow(kx[i][j],2)*eps_TMM_O[0][-1]/eps_TMM_E[0][-1])

        Ap, Bp = get_A_B(d_TMM, N_layers, wl[0][j], eps_TMM_O, eps_TMM_E, kx[i][j], mu)
        
        

        ApEMT_st, BpEMT_st = get_A_B(d_EMT, 3, wl[0][j], eps_EMT_O_st, eps_EMT_E_st, kx[i][j], mu)
        ApEMT_i1, BpEMT_i1 = get_A_B(d_EMT, 3, wl[0][j], eps_EMT_O_i1, eps_EMT_E_i1, kx[i][j], mu)

        Rp_TM = pow(np.abs(Bp[0][0]/Ap[0][0]),2)
        Tr_TM = np.real(kz_end[i][j]/eps_TMM_O[-1][0])/np.real(kz_air[i][j]/eps_TMM_O[0][0])*pow(np.abs(Ap[0][-1]/Ap[0][0]), 2)
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
#print "==="
#print
#print "Rp:"
#print Rp
#print "Wavelength:"
#print wl[0]
#print "Theta_i"
#print theta_i[0]
#print
#print "==="

#for i in range(len(eps_EMTE_i1[0])):
#    ni = cmath.sqrt(eps_EMTE_i1[0][i])
#    ki = ni.imag
#    nst = cmath.sqrt(eps_EMTE_st[0][i])
#    kst = nst.imag
#
#    Pd[0][i] = wl[0][i]/(4.0*pi*ki)  # EMT improved
#    Pd[1][i] = wl[0][i]/(4.0*pi*kst) # EMT standard


#-------------------------------

# POST PROCESSING
#-------------------------------
#for 2d plots
R_3 = np.zeros((3, Wavelength_resolution), dtype=float)
T_3 = np.zeros((3, Wavelength_resolution), dtype=float)
A_3 = np.zeros((3, Wavelength_resolution), dtype=float)

for i in range(len(Rp[0])):
    R_3[0][i] = Rp[0][i]          # TMM
    T_3[0][i] = Tr[0][i]
    A_3[0][i] = Ab[0][i]

    R_3[1][i] = RpEMT_st[0][i]    # EMT standard
    T_3[1][i] = TrEMT_st[0][i]
    A_3[1][i] = AbEMT_st[0][i]

    R_3[2][i] = RpEMT_i1[0][i]    # EMT improved
    T_3[2][i] = TrEMT_i1[0][i]
    A_3[2][i] = AbEMT_i1[0][i]

emtE_ir = eps_EMTE_i1[0].real
emtO_ir = eps_EMTO_i1[0].real
emtE_ii = eps_EMTE_i1[0].imag
emtO_ii = eps_EMTO_i1[0].imag

emtE_str = eps_EMTE_st[0].real
emtO_str = eps_EMTO_st[0].real
emtE_sti = eps_EMTE_st[0].imag
emtO_sti = eps_EMTO_st[0].imag

condEr = eps_condE[0].real
condOr = eps_condO[0].real
condEi = eps_condE[0].imag
condOi = eps_condO[0].imag

#fixed wl - many angles: !!! REDO!!!!!!!!
#for i in range(len(RpEMT_st[0])):
#    R_3[0][i] = RpEMT_st[i][0] #EMT_st
#    T_3[0][i] = TrEMT_st[i][0]
#
#    R_3[1][i] = RpEMT_i1[i][0] #EMT_i1
#    T_3[1][i] = TrEMT_i1[i][0]
#
#    R_3[2][i] = Rp[i][0] #TMM
#    T_3[2][i] = Tr[i][0]

directory = '../plots/September/2/'
prop1 = "Transmittance"
prop2 = "Reflectance_Contour"
prop3 = "Absorbtance"

EMTi = '_EMTi_eps.png'
EMTst = '_EMTst_eps.png'
l1 = "Extraordinary, real"
l2 = "Ordinary, real"
l3 = "Extraordinary, imaginary"
l4 = "Ordinary, imaginary"

fnEMTi =  directory +"diel="  +  str(d_dielectric*1E9)+ "_bulk=" +str(d_bulk*1E9) + EMTi
fnEMTst =  directory +"diel=" + str(d_dielectric*1E9)+ "_bulk=" +str(d_bulk*1E9) + EMTst

fn1 = directory +"diel=" + str(d_dielectric*1E9)+ "_bulk=" +str(d_bulk*1E9) + prop1+ ".png"
fn2 = directory +"diel=" + str(d_dielectric*1E9)+ "_bulk=" +str(d_bulk*1E9) + prop2+ ".png"
fn21= directory +"diel=" + str(d_dielectric*1E9)+ "_bulk=" +str(d_bulk*1E9) + prop3+ ".png"

fn3 = directory + "diel=" + str(d_dielectric*1E9)+"bulk_eps_test1111.png"
fn4 = directory + "diel=" + str(d_dielectric*1E9)+"CONDUCTION_eps.png"
fn5 = directory +"diel="+str(d_dielectric*1E9)+"_bulk=" +str(d_bulk*1E9)+"_Pd_"

argsEpsI4  = [emtE_ir, emtO_ir, emtE_ii, emtO_ii, fnEMTi, l1, l2, l3, l4, 3]
argsEpsST4 = [emtE_str, emtO_str, emtE_sti, emtO_sti, fnEMTst, l1, l2, l3, l4, 3]

argsEpsCOND = [condEr, condOr, condEi, condOi, fn4, l1, l2, l3, l4, 3]
argsEpsBulk =[eps_bulkO[0].real, eps_bulkO[0].imag,"Epsilon real","Epsilon imaginary", fn3, 3]

postprocess.basic_info(material_name(), N_layers, N_periods)

#postprocess.doFigure_Eps4(wl[0], argsEpsI4)
#postprocess.doFigure_Eps4(wl[0], argsEpsST4)

#postprocess.doFigure_RTA(wl[0], T_3[0], T_3[2], T_3[1], fn1, prop1,1)
#postprocess.doFigure_RTA(wl[0], R_3[0], R_3[2], R_3[1], fn2, prop2,1)
#postprocess.doFigure_RTA(wl[0], A_3[0], A_3[2], A_3[1], fn21, prop3,1)

#postprocess.doFigure_Eps4(wl[0], argsEpsCOND)
#postprocess.doFigure_Eps2(wl[0], argsEpsBulk)

# ======== Write to a file:============
args = [wl[0], emtE_ir, emtO_ir, emtE_ii, emtO_ii]
argsA = [wl[0],A_3[0], A_3[2], A_3[1]]
argsR = [wl[0],R_3[0], R_3[2], R_3[1]]

argsB = [wl[0], eps_bulkO[0].real, eps_bulkO[0].imag]
argsC = [wl[0], condEr, condOr, condEi, condOi]

titleAR = 'Wavelength,[nm]'+ '\t'+'TMM'+ '\t'+ "EMTi"+ '\t'+"EMTst"
titleEPSi = 'Wavelength,[nm]'+ '\t'+'ExtraO_real'+ '\t'+ "Ordinary_real"+ '\t'+"ExtraO_imag" \
 + '\t'+"Ordinary_imag"
titleBULK = 'Wavelength,[nm]'+ '\t'+'EPS_real'+ '\t'+ "EPS_imag"
titleCOND = 'Wavelength,[nm]'+ '\t'+'ExtraO_real'+ '\t'+ "Ordinary_real"+ '\t'+"ExtraO_imag" \
 + '\t'+"Ordinary_imag"
#postprocess.writeToFile(fn3[:-4] +".xls", titleBULK, argsB) #FIG2 - BULK
#postprocess.writeToFile(fn4[:-4] +".xls", titleCOND, argsC) #FIG2 - COND

#postprocess.writeToFile(fnEMTi[:-4] +".xls", titleEPSi, args)    #FIG3 -6
#postprocess.writeToFile(fn2[:-4] +".xls", titleAR, argsR)        #FIG4 -7
#postprocess.writeToFile(fn21[:-4] +".xls", titleAR, argsA)       #FIG5 -8


postprocess.doContourPlot(wl[0], theta_i[0], Rp, 'Contour_Rp_100&100nm_10x18_semilog')
postprocess.doContourPlot(wl[0], theta_i[0], Tr, 'Contour_Tr_100&100nm_10x18_semilog')
#postprocess.doContourPlot(wl[0], theta_i[0], Ab, 'Contour_Ab_1001')

