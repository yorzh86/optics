from __future__ import division
import numpy as np
import cmath
import math
import pylab as pl


def bulk_Wolf(x1):
    # fn calculates epsilon of bulk part of Bi2Se3 for 300K
    # using 3x-Lorentz-Drude model
    # given energy in [nm] (units conversion performed below)
    # source: M.S.Wolf"INFRARED AND OPTICAL STUDIES OF TOPOLOGICAL INSULATORS"
    # to use: print bulk_Wolf(500)
    ##(-8.61086993518+5.93331696103j)

    w = 1.2398/x1*1E3  #[nm to eV]
    x = w*8065.54429   #[eV to cm^-1]

    eps_inf = 1.0
    wpD     = 5651.5
    w0D     = 0.0
    gammaD  = 111.86
    wp = np.array([ [66024] ])
    w0 = np.array([ [8386.6] ])
    gamma = np.array([ [10260] ])

    eps1 = eps_inf - math.pow(wpD,2)/(math.pow(x,2)+1j*gammaD*x)
    eps2 = 0.0+0.0j

    # Attention - I do not divide by 4pi to get correct results
    for i in range(len(wp[0])):
         eps2 += math.pow(wp[0][i],2)/(math.pow(w0[0][i],2)-math.pow(x,2)-1j*gamma[0][i]*x)
    return eps1+eps2

#========
# Uncomment to test. Compare with figure 5.12 300K drude (page 62)
# Attention, change fn to work with cm^-1(delete all unit conversion(w,x) and change
# x1 to x)
#========
#~ R = np.zeros((1, 1000), dtype=complex)
#~ freq = np.linspace(1,25000, 1000)
#~ for i in range(len(freq)):
    #~ a = get_eps_Bi2Te3_bulk_Wolf(freq[i])
    #~ R[0][i] = math.pow(np.abs((cmath.sqrt(a)-1)/(cmath.sqrt(a)+1)),2)
#~ #
#~ b = get_eps_Bi2Te3_bulk_Wolf(100)
#~ c = math.pow(np.abs((cmath.sqrt(b)-1)/(cmath.sqrt(b)+1)),2)
#~ print 'R, at 100 cm^-1, (should be ~0.95+):', c
#~ #
#~ b1 = get_eps_Bi2Te3_bulk_Wolf(300)
#~ c1 = math.pow(np.abs((cmath.sqrt(b1)-1)/(cmath.sqrt(b1)+1)),2)
#~ print 'R, at 300 cm^-1, (should be ~0.95+):', c1
#~ #
#~ b3 = get_eps_Bi2Te3_bulk_Wolf(500)
#~ c3 = math.pow(np.abs((cmath.sqrt(b3)-1)/(cmath.sqrt(b3)+1)),2)
#~ print 'R, at 500 cm^-1, (should be ~0.9):', c3
#~ #
#~ b2 = get_eps_Bi2Te3_bulk_Wolf(1000)
#~ c2 = math.pow(np.abs((cmath.sqrt(b2)-1)/(cmath.sqrt(b2)+1)),2)
#~ print 'R, at 1000 cm^-1, (should be ~0.5):', c2
#~ pl.plot(freq, R[0])
