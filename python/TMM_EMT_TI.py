from __future__ import division
import numpy as np
import sys
from math import *
import cmath
from TMM_aniso import get_A_B
from postprocess import *
#from Poynting import calcPoynting
import datetime
import os

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
def calculateRpTrAb(substrate, ti, total, wl_r=10, angle_r=10):

    N_periods = int(round(total/(substrate+ti)))
    N_layers = int(3+N_periods*4)
    Wavelength_resolution = wl_r
    Angle_resolution = angle_r

    # Permeability
    mu = np.ones((1, N_layers), dtype=float)

    # Thickness of layers for TMM
    d_TMM = np.zeros((1, N_layers), dtype = float)
    d_air = 0
    d_dielectric = substrate*1E-9
    d_cond = d_conduct()
    d_bulk = (ti*1E-9-d_cond*2) #should be 43 layers instead of 42 


    aa_ = [d_dielectric, d_cond, d_bulk, d_cond]
    d_TMM[0] = d_air
    aa_ = np.tile(aa_, N_periods)
    for i in range(len(d_TMM[0])-3):
        d_TMM[0][i+1] = aa_[i]  # + dielectric before air
    d_TMM[0][-2] = d_dielectric
        
    
    diff_norm = np.cumsum(d_TMM[0])/total

    # Thickness of layers for EMT
    d_EMT = np.array([ [d_air, (d_cond+d_bulk+d_dielectric+d_cond)*N_periods, d_air]  ])

    ff = (d_cond*2+d_bulk)/(d_cond*2+d_bulk+d_dielectric)
    ffm = d_bulk/(d_cond*2+d_bulk+d_dielectric)
    ffd = d_dielectric/(d_cond*2+d_bulk+d_dielectric)

    # Setup wavelengths
    wl = np.zeros((1,Wavelength_resolution), dtype=float)
    wl[0] = np.logspace(np.log10(500), np.log10(20000), Wavelength_resolution)
    #wl[0] = np.linspace(500, 2000, Wavelength_resolution)
    #wl= np.array([[500, 1000]])

    # Select angles of incidence
    theta_i= np.zeros((1,Angle_resolution), dtype=float)
    theta_i[0] = np.linspace(0,89.9,Angle_resolution)
    #theta_i= np.array([[0, 60]])

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
    Pd = np.zeros((len(theta_i[0]), len(wl[0])), dtype=float)

    # Auxiliary arrays to calculate A,B coefficients
    kx = np.zeros((len(theta_i[0]),len(wl[0])), dtype=complex)
    kz_air = np.zeros((len(theta_i[0]),len(wl[0])), dtype=complex)
    kz_end = np.zeros((len(theta_i[0]),len(wl[0])), dtype=complex)
    kz_EMT = np.zeros((len(theta_i[0]),len(wl[0])), dtype=complex)
    
    # what is this?
#    step = 100
#    z_min = -0.5
#    z_max = 2
#    slab = total
#    z_range = np.zeros((len(d_TMM), step*2.5), dtype = float)
#    z_range[0] = np.linspace(z_min,-1.0/step,step*2.5)
#    
#    ii = 1
#    for ii in range(len(d_TMM)):
#        z_range[ii] = np.linspace()
#   for zz = 2:length(diff)
#       z_range{zz} = linspace(sum(diff(1:(zz-1)))/slab,sum(diff(1:zz))/slab-(diff(zz)/slab)/step,step)
#
#   z_range{length(diff)+1} = linspace(1,z_max,step)
#   z_norm = [z_range{1:end}]

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
            eps_TMM_O = np.append(eps_TMM_O, eps_dO[0][j])
            eps_TMM_O = np.append(eps_TMM_O, eps_air)
            eps_TMM_O = np.array([eps_TMM_O])

            eps_period_E = np.tile(eps_period_E, N_periods)
            eps_TMM_E = np.append(eps_air, eps_period_E)
            eps_TMM_E = np.append(eps_TMM_E, eps_dE[0][j])
            eps_TMM_E = np.append(eps_TMM_E, eps_air)
            eps_TMM_E = np.array([eps_TMM_E])
            
            #print eps_EMTO_i1[0]
            #sys.exit()
            
            om = 2*pi*c0/wl[0]*1E9 #wl[0][j]
            k0 = om/c0
            kx[i][j] = k0[j]*sin(theta_i[0][i]*pi/180)
            kz_air[i][j] = cmath.sqrt(pow(k0[j],2)-pow(kx[i][j],2))
            kz_end[i][j] = cmath.sqrt(pow(k0[j],2)*eps_TMM_O[0][-1]-pow(kx[i][j],2)*eps_TMM_O[0][-1]/eps_TMM_E[0][-1])
            kz_EMT[i][j] = cmath.sqrt(pow(k0[j],2)*eps_EMTO_i1[0][j]-pow(kx[i][j],2)*eps_EMTO_i1[0][j]/eps_EMTE_i1[0][j])         
            
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
            
            Pd[i][j] = 1/(2*np.imag(kz_EMT[i][j]))*1e9
            
            #a[][],b[][],c[][],d[][],e[][],f[][] = calcPoynting(z_norm, d_norm, total_d, kx[i][j], eps_TMM_O, eps_TMM_E, wl[0][j], A_TE=0, B_TE=0, A_TM =0, B_TM=0)

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


    today = datetime.date.today()
    month_day = today.strftime('%b')+'/'+ today.strftime('%d')+'/'
    
    #create folder
    path = '../plots/'+month_day
    if (os.path.isdir(path) == False):
        os.makedirs(path)
    
    foldername = '../plots/'+month_day


    directory = foldername
    prop1 = "Transmittance"
    prop2 = "Reflectance_Contour"
    prop3 = "Absorbtance"

    EMTi = '_EMTi_eps.png'
    EMTst = '_EMTst_eps.png'
    l1 = "Extraordinary, real"
    l2 = "Ordinary, real"
    l3 = "Extraordinary, imaginary"
    l4 = "Ordinary, imaginary"

    fnEMTi =  directory +str(material_name()[:6])+"_diel="  +  str(d_dielectric*1E9)+ "_bulk=" +str(d_bulk*1E9) + EMTi
    fnEMTst =  directory +"diel=" + str(d_dielectric*1E9)+ "_bulk=" +str(d_bulk*1E9) + EMTst

    fn1 = directory +"diel=" + str(d_dielectric*1E9)+ "_bulk=" +str(d_bulk*1E9) + prop1+ ".png"
    fn2 = directory +"diel=" + str(d_dielectric*1E9)+ "_bulk=" +str(d_bulk*1E9) + prop2+ ".png"
    fn21= directory +"diel=" + str(d_dielectric*1E9)+ "_bulk=" +str(d_bulk*1E9) + prop3+ ".png"

    fn3 = directory + "diel=" + str(d_dielectric*1E9)+"bulk_eps_test1111.png"
    fn4 = directory + "diel=" + str(d_dielectric*1E9)+"CONDUCTION_eps.png"
    fn5 = directory +"diel="+str(d_dielectric*1E9)+"_bulk=" +str(d_bulk*1E9)+"_Pd_"
    fn6 = directory +"diel="+str(d_dielectric*1E9)+"_bulk=" +str(d_bulk*1E9)+"_Rp_"
    fn7 = directory +"diel="+str(d_dielectric*1E9)+"_bulk=" +str(d_bulk*1E9)+"_Tr_"

    argsEpsI4  = [emtE_ir, emtO_ir, emtE_ii, emtO_ii, fnEMTi, l1, l2, l3, l4, 3]
    argsEpsST4 = [emtE_str, emtO_str, emtE_sti, emtO_sti, fnEMTst, l1, l2, l3, l4, 3]

    argsEpsCOND = [condEr, condOr, condEi, condOi, fn4, l1, l2, l3, l4, 3]
    argsEpsBulk =[eps_bulkO[0].real, eps_bulkO[0].imag,"Epsilon real","Epsilon imaginary", fn3, 3]

    basic_info(material_name(), N_layers, N_periods)

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
    argsRP =[wl[0], theta_i[0], Rp] 
    argsTr =[wl[0], theta_i[0], Tr] 

    titleAR = 'Wavelength,[nm]'+ '\t'+'TMM'+ '\t'+ "EMTi"+ '\t'+"EMTst"
    titleEPSi = 'Wavelength,[nm]'+ '\t'+'ExtraO_real'+ '\t'+ "Ordinary_real"+ '\t'+"ExtraO_imag" \
     + '\t'+"Ordinary_imag"
    titleBULK = 'Wavelength,[nm]'+ '\t'+'EPS_real'+ '\t'+ "EPS_imag"
    titleRP = " Structure of the file: 1. all theta_i, 2. all wl, 3. each row has values of Rp for single angle for each wl " 
    titleTr = " Structure of the file: 1. all theta_i, 2. all wl, 3. each row has values of Tr for single angle for each wl " 
    
    titleCOND = 'Wavelength,[nm]'+ '\t'+'ExtraO_real'+ '\t'+ "Ordinary_real"+ '\t'+"ExtraO_imag" \
     + '\t'+"Ordinary_imag"
    #postprocess.writeToFile(fn3[:-4] +".xls", titleBULK, argsB) #FIG2 - BULK
    #postprocess.writeToFile(fn4[:-4] +".xls", titleCOND, argsC) #FIG2 - COND

    #writeToFile(fnEMTi[:-4] +".xls", titleEPSi, args)    #FIG3 -6
    #postprocess.writeToFile(fn2[:-4] +".xls", titleAR, argsR)        #FIG4 -7
    #postprocess.writeToFile(fn21[:-4] +".xls", titleAR, argsA)       #FIG5 -8

    cfg = 'cfg.'+str(d_dielectric*1E9)[:3]+'x'+str(d_bulk*1E9)[:5]+'nm__'
    rsl = 'res.'+str(Wavelength_resolution)+'x'+str(Angle_resolution)


    styles = ['viridis', 'plasma', 'inferno', 'magma',
              'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
                'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
                'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                'hot', 'afmhot', 'gist_heat', 'copper',
                'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic',
                'Pastel1', 'Pastel2', 'Paired', 'Accent',
                'Dark2', 'Set1', 'Set2', 'Set3',
                'tab10', 'tab20', 'tab20b', 'tab20c',
                'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
                'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
                'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']


#    for element in styles:
#        doContourPlot(wl[0], theta_i[0], Rp, element+'_'+rsl+'.png', element)
    
    doContourPlot(wl[0], theta_i[0], Rp, foldername, str(material_name()[:6])+'_Rp_'+
                              cfg +rsl, 1, style = 'jet')
    doContourPlot(wl[0], theta_i[0], Tr, foldername, str(material_name()[:6])+'_Tr_'+
                              cfg +rsl, 2, style = 'jet')
#    doContourPlot(wl[0], theta_i[0], Pd, foldername, str(material_name()[:6])+'_Pd_'+
#                             cfg +rsl, 3)
    
    writeToFile_Contour(fn6[:-1]+ ".xls", titleRP, argsRP)
    writeToFile_Contour(fn7[:-1]+ ".xls", titleTr, argsTr)
    return

#testing
calculateRpTrAb(100, 100, 2000, 200 , 180)