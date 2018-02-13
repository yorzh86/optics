from __future__ import division
import numpy as np
import cmath
import math
import pylab as pl

wl_eps_real = np.array([
    [0.024621213, 16.02564],
    [0.077651516, 16.179487],
    [0.1344697, 16.461538],
    [0.17992425, 16.717949],
    [0.2935606, 16.76923],
    [0.3219697, 16.384615],
    [0.35416666, 15.897436],
    [0.37689394, 15.076923],
    [0.40151516, 14.128205],
    [0.42992425, 12.692307],
    [0.45454547, 11.589744],
    [0.47537878, 10.230769],
    [0.49810606, 8.897436],
    [0.5113636, 7.794872],
    [0.5359849, 6.6666665],
    [0.56439394, 5.3589745],
    [0.6003788, 3.974359],
    [0.625, 2.974359],
    [0.6685606, 1.8974359],
    [0.7026515, 0.974359],
    [0.7462121, 0.33333334],
    [0.8125, -0.35897437],
    [0.8863636, -0.35897437],
    [0.9943182, 0.102564104],
    [1.0890151, 0.5897436],
    [1.1723485, 0.84615386],
    [1.2386364, 1.0769231],
    [1.3674242, 1.1025641],
    [1.4469697, 0.51282054],
    [1.501894, -0.15384616],
    [1.5397727, -0.5641026],
    [1.5719697, -0.82051283],
    [1.810606, -0.7692308],
    [1.8579545, -1.1025641],
    [1.9166666, -1.3076923],
    [1.9640151, -1.6153846],
    [2.0113637, -1.8461539],
    [2.1742425, -1.8717948],
    [2.2026515, -1.6410257],
    [2.3333333, -1.6666666],
    [2.3655303, -1.4358975],
    [2.5132575, -1.4102564],
    [2.6496212, -1.1025641],
    [2.7632575, -0.84615386],
    [2.8806818, -0.6666667],
    [3.0359848, -0.41025642],
    [3.2405303, -0.35897437],
    [3.4734848, -0.17948718],
    [3.6401515, -0.025641026],
    [3.905303, 0.23076923],
    [4.2708335, 0.53846157],
    [4.998106, 0.61538464],
])

wl_eps_imag = np.array([
    [0.018939395, 1.1025641],
    [0.071969695, 1.9230769],
    [0.115530305, 2.6923077],
    [0.1534091, 3.6923077],
    [0.19318181, 4.769231],
    [0.24431819, 6.1025643],
    [0.2840909, 7.74359],
    [0.31060606, 8.897436],
    [0.34280303, 10.307693],
    [0.375, 11.897436],
    [0.41098484, 13.307693],
    [0.44507575, 14.641026],
    [0.48674244, 15.461538],
    [0.5265151, 15.461538],
    [0.5833333, 14.48718],
    [0.62310606, 13.51282],
    [0.67045456, 12.538462],
    [0.71022725, 11.358974],
    [0.7537879, 10.358974],
    [0.7859849, 9.461538],
    [0.82575756, 8.48718],
    [0.8863636, 7.3589745],
    [0.9375, 6.5897436],
    [0.9886364, 6.1025643],
    [1.032197, 5.74359],
    [1.094697, 5.3846154],
    [1.1590909, 5.1025643],
    [1.2121212, 4.8717947],
    [1.2689394, 5.1282053],
    [1.3238636, 5.282051],
    [1.375, 5.6153846],
    [1.4280303, 6.0],
    [1.4810606, 5.974359],
    [1.5587121, 5.4358974],
    [1.6155303, 4.8717947],
    [1.6761364, 4.5897436],
    [1.7291666, 4.3076925],
    [1.7992424, 4.3076925],
    [1.8920455, 4.25641],
    [2.003788, 3.5384614],
    [2.0984848, 2.9487178],
    [2.2007575, 2.3846154],
    [2.3238637, 1.8717948],
    [2.4242425, 1.6153846],
    [2.5795455, 1.1794872],
    [2.6666667, 0.9230769],
    [2.7708333, 0.7692308],
    [2.8863637, 0.7692308],
    [2.9375, 0.51282054],
    [3.2443182, 0.53846157],
    [3.3238637, 0.33333334],
    [3.873106, 0.3846154],
    [4.380682, 0.3846154],
    [4.998106, 0.30769232],
])

def bulk_Wolf(x1):
    # fn calculates epsilon of bulk part of Bi2Se3 for 300K
    # using 3x-Lorentz-Drude model
    # given energy in [nm] (units conversion performed below)
    # source: M.S.Wolf"INFRARED AND OPTICAL STUDIES OF TOPOLOGICAL INSULATORS"
    # to use: print bulk_Wolf(500)
    #(0.689172927554+0.0609276491179j)

    w = 1.2398/x1*1E3  #[nm to eV]
    x = w*8065.54429   #[eV to cm^-1]
    eps_inf = 1.0
    wpD     = 908.66
    w0D     = 0.0
    gammaD  = 7.43
    wp = np.array([ [675.9, 100, 11249] ])
    w0 = np.array([ [63.03, 126.94, 2029.5] ])
    gamma = np.array([ [17.5, 10, 3920.5] ])

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
#R = np.zeros((1, 1000), dtype=complex)
#freq = np.linspace(1,25000, 1000)
#for i in range(len(freq)):
#    a = get_eps_Bi2Se3_bulk_Wolf(freq[i])
#    R[0][i] = math.pow(np.abs((cmath.sqrt(a)-1)/(cmath.sqrt(a)+1)),2)
#
#b = get_eps_Bi2Se3_bulk_Wolf(100)
#c = math.pow(np.abs((cmath.sqrt(b)-1)/(cmath.sqrt(b)+1)),2)
#print 'R, at 100 cm^-1, (should be ~0.95):', c
#
#b1 = get_eps_Bi2Se3_bulk_Wolf(300)
#c1 = math.pow(np.abs((cmath.sqrt(b1)-1)/(cmath.sqrt(b1)+1)),2)
#print 'R, at 300 cm^-1, (should be ~0.38):', c1
#
#b3 = get_eps_Bi2Se3_bulk_Wolf(500)
#c3 = math.pow(np.abs((cmath.sqrt(b3)-1)/(cmath.sqrt(b3)+1)),2)
#print 'R, at 500 cm^-1, (should be ~0.43):', c3
#
#b2 = get_eps_Bi2Se3_bulk_Wolf(1000)
#c2 = math.pow(np.abs((cmath.sqrt(b2)-1)/(cmath.sqrt(b2)+1)),2)
#print 'R, at 1000 cm^-1, (should be ~0.5):', c2
#pl.plot(freq, R[0])


def bulk_Yin(x1):
    # fn interates through known values of eps(eV) and interpolates
    # source: Yin et al. Figure 5:"Plasmonics of Topological Insulators at Optical Frequencies"
    # to simulate bulk - we look at interband contributions (solid lines)

    # How to use:
    #x = 490 #nm -> 2.53 eV
    #a = bulk_Yin(x)
    #print "Bi2Se3 bulk epsilon(real, imag) for",x, "[nm] is:", a

    #convert nm to eV
    x = 1.2398/x1*1E3 #WIKI correct!!
    for i in range(len(wl_eps_imag)):
        if (wl_eps_imag[i][0]>x):
            low_wl_i =  wl_eps_imag[i-1][0]
            up_wl_i  =  wl_eps_imag[i][0]
            low_eps_i = wl_eps_imag[i-1][1]
            up_eps_i  =  wl_eps_imag[i][1]
            break
    eps_imag = low_eps_i+(x-low_wl_i)*(up_eps_i-low_eps_i)/(up_wl_i-low_wl_i)

    for j in range(len(wl_eps_real)):
        if (wl_eps_real[j][0]>x):
            low_wl_r =  wl_eps_real[j-1][0]
            up_wl_r  =   wl_eps_real[j][0]
            low_eps_r = wl_eps_real[j-1][1]
            up_eps_r  =  wl_eps_real[j][1]
            break
    eps_real = low_eps_r+(x-low_wl_r)*(up_eps_r-low_eps_r)/(up_wl_r-low_wl_r)

    return [eps_real, eps_imag]
