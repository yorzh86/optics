from __future__ import division
import numpy as np
from math import *
import cmath

def get_A_B(d, N_layers, wl, eps_O, eps_E, kx, mu):

    c0 = 2.99792458E8
    eps0 =8.854187817e-12

    # axiliary array;includes PD^-1D elements
    kz_p= np.zeros((1, N_layers), dtype=complex)
    kz_n= np.zeros((1, N_layers), dtype=complex)
    # wave component
    z = np.zeros((1, N_layers), dtype=complex)
    # impedance
    P  = np.zeros((N_layers-1, 2, 2), dtype = complex)
    # propagation matrix
    DD = np.zeros((N_layers-1, 2, 2), dtype = complex)
    # dynamical matrix for s polarization
    M = np.ones((N_layers-1, 2, 2), dtype = complex)
    # auxiliary M matrix
    Mp = np.ones((N_layers, 2, 2), dtype = complex)
    # the p-polarization M matrix (matlab M_prod)
    Ap = np.zeros((1, N_layers), dtype=complex)
    Bp = np.zeros((1, N_layers), dtype=complex)
    # matrixes for A and B coefficients

    for l in range(N_layers):
        om = 2*pi*c0/wl*1.0E9
        k0 = om/c0
        kz_p[0][l] = cmath.sqrt(eps_O[0][l]*mu[0][l]*pow(k0,2) - pow(kx,2)*eps_O[0][l]/eps_E[0][l])
        kz_n[0][l] = -kz_p[0][l]
        if (kz_p[0][l].imag < 0):
            kz_p[0][l] = np.conj(kz_p[0][l])

        z[0][l] = 1.0/(om*eps0)*kz_p[0][l]/eps_O[0][l]

    Mp[N_layers-1] = np.identity(2)
    for l in range(N_layers-2, -1, -1):
        r_ij = (z[0][l]-z[0][l+1])/(z[0][l]+z[0][l+1])
        r_ji = (z[0][l+1]-z[0][l])/(z[0][l]+z[0][l+1])
        t_ij = (2*z[0][l])/(z[0][l]+z[0][l+1])
        t_ji = (2*z[0][l+1])/(z[0][l]+z[0][l+1])

        DD[l] = np.array ([
            [1.0, -r_ji],
            [r_ij, t_ij*t_ji-r_ij*r_ji],
            ])*1.0/t_ij

        #~ print "------"
        #~ print "l:", l
        #~ print "kz_p[0][l]:", kz_p[0][l]
        #~ print "d[0][l]:", d[0][l]
        #~ print "product:", -1j*kz_p[0][l]*d[0][l]
        #~ print "exp(product):", cmath.exp(-1j*kz_p[0][l]*d[0][l]) Error here!!!
        #~ print "------"

        P[l] = np.array([
            [cmath.exp(-1j*kz_p[0][l]*d[0][l]), 0.0],
            [0.0, cmath.exp(-1j*kz_n[0][l]*d[0][l])],
            ])
        M[l] = P[l].dot(DD[l])
        Mp[l] = np.dot(M[l],Mp[l+1])

    Bp[0][N_layers-1] = 0.0
    Ap[0][N_layers-1] = 1.0/Mp[0][0][0]
    Bp[0][0] = Mp[0][1][0]*Ap[0][N_layers-1]
    Ap[0][0]= 1.0

    for l in range (N_layers-1, -1, -1):
        Ap[0][l] = Mp[l][0][0]*Ap[0][N_layers-1]
        Bp[0][l] = Mp[l][1][0]*Ap[0][N_layers-1]

    return Ap, Bp
