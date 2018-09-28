from __future__ import division
import sys
import cmath
import numpy as np
import cmath
f

def Poynting(z_norm, d, d_norm, kx, eps_O, eps_E, wl, theta_i, A_TM =0, B_TM=0, A_TE=0, B_TE=0):
    #fn returns S_xTE, S_zTE, S_xTM, S_zTM, Tr_TM, E_TMmag
    c0 = 2.99792458E8
    eps0 =8.854187817e-12
    mu0 = 4.0*np.pi/1E7  # vac permeability
    
    om = 2.0*np.pi*c0/wl/1E9   #wl in nanometers
    k0 = om/c0
    
    z = z_norm*d
    z_0 = np.zeros((0, d*d_norm), dtype=float)
    
    #some loop for indexing that I dont undestand
    # is it to find the number of layers?or number of last layer?
#    for ii = 1:length(diff_norm)-1
#    if (z_norm < diff_norm(1))
#        ind = 1;
#    elseif (z_norm < diff_norm(ii+1))&&(z_norm >= diff_norm(ii))
#        ind = ii+1;
#    elseif (z_norm >= diff_norm(end))
#        ind = length(diff_norm)+1;
#        end
#    end
#    
    
    kz_air = cmath.sqrt(k0**2 - kx**2)
    kzTM[i] = cmath.sqrt(k0**2*eps_TMM_O[0][-1]-(kx[i][j]**2)*eps_TMM_O[0][-1]/eps_TMM_E[0][-1])
    if (kz_p[0][l].imag < 0):
        kz_end[0][l] = np.conj(kz_p[0][l])
   
    kzTM(ind)=sqrt(k0^2*ep_o(ind)-ep_o(ind)/ep_e(ind)*k_x^2)
   if (imag(kzTM(ind))<0):
       kzTM(ind)=-sqrt(k0^2*ep_o(ind)-ep_o(ind)/ep_e(ind)*k_x^2)

        
    kzTE(ind) = sqrt(ep_o(ind)*k0^2-k_x^2)
    if (imag(kzTE[i])<0):
        kzTE[i] = -sqrt(eps_O[i]*k0**2-kx**2)

        
    #for TM
    if (A_TE==0 and B_TE ==0):        
        Hy = (A_TM[i]*np.exp(1j*kz_TM[i]*(z-z_0[i]))+B_TM[i]*np.exp(-1j*kz_TM[i] \
         *(z-z_0[i])))*np.exp(1j*kx)
        Hy_conj = np.conj(Hy)
        
        dHy_dz = (kz_TM[i])*A_TM[i]*np.exp(1j*kz_TM[i])*(z-z_0[i])-kzTM[i]*B_TM[i] \ 
                        *np.exp(-1j*kzTM[i]*(z-z_0[i])))*np.exp(1j*kx)
        dHy_dx = kx*Hy
        
        Hy0 = (A_TM[0]+B_TM[0])*np.exp(1j*kx)
        dHy_dz0 = (kzTM[0]*A_TM[0]-kzTM[0]*B_TM[0])*np.exp(1j*kx)
        dHy_dx0 = k_x*Hy0
        
        Ex0 = 1.0/(om*eps0*eps_O[0]*eps_E[0])*(eps_E[0]*dHy_dz0)
        Ez0 = -1.0/(om*eps0*eps_O[0]*eps_E[0])*(eps_O[0]*dHy_dx0)
        E_TM0 = np.sqrt(Ex0*np.conj(Ex0)+Ez0*np.conj(Ez0))
        
        Tr_TM = np.real(kz_TM[i]/eps_O[i])/np.real(kz_air/eps_O[0] \
                       *np.abs(A_TM[i]/A_TM[0])**2    # return this
        E_TMmag = Hy**2/E_TM0**2      # return this
        
        Ex = 1.0/(om*eps0*eps_O[i]*eps_E[i])*(eps_E[i]*dHy_dz)
        Ez = -1.0/(om*eps0*eps_O[i]*eps_E[i])*(eps_O[i]*dHy_dx)
        S_xTM = -0.5*real(Ez*Hy_conj) # return this
        S_zTM = 0.5*real(Ex*Hy_conj)  # return this
    
    #for TE
    else if (A_TM==0 and B_TM==0):
        Ey = (A_TE[i]*np.exp(1j*kzTE[i]*(z-z_0[i]))+B_TE[i]*np.exp(-1j*kzTE[i]*(z-z_0[i])))*np.exp(1j*kx)
        dEy_dz = (1j*kzTE[i]*A_TE[i]*np.exp(1j*kzTE[i]*(z-z_0[i]))-1j*kzTE[i]*B_TE[i]*np.exp(-1j*kzTE[i])*(z-z_0[i])))*np.exp(1j*kx)
        dEy_dx = 1j*kx*(A_TE[i]*np.exp(1j*kzTE[i]*(z-z_0[i]))+B_TE[i]*np.exp(-1j*kzTE[i]*(z-z_0[i])))*np.exp(1j*kx)
        Hx = -1.0/(1j*om*mu0)*dEy_dz
        Hz = 1.0/(1j*om*mu0)*dEy_dx
        Hx_conj = conj(Hx)
        Hz_conj = conj(Hz)
        
        S_xTE = 0.5*np.real(Ey*Hz_conj)   # return this
        S_zTE = -0.5*np.real(Ey*Hx_conj)  # return this
        
        Ey0 = (A_TE(1)+B_TE(1))*exp(1i*k_x)
        E_TEmag = Ey^2/Ey0^2           # return this
        
    return A_TM, B_TM, A_TE, B_TE, 

    