from __future__ import division

import numpy as np
import scipy.optimize as opt
import math

#for red line (through plane)
def drude_E_eps(w1):
    #Calculates epsilon, given energy[nm]. Other parameters w0,wp,etc can also be changed
    #a = drude_E_eps(500)
    #print a

    w0 = 0.9310534
    wp = 1.21041352
    eps_inf = 1.58780666
    gamma= 0.1

    #convert [nm] to [eV] #WIKI correct!!
    w = 1.2398/w1*1E3
    eps_r = (math.pow(wp,2)*(math.pow(w0,2) - math.pow(w,2))) / \
        ((math.pow(w0,2) - math.pow(w,2))**2 + math.pow(w*gamma, 2))+eps_inf
    eps_i = math.pow(wp,2)*gamma*w/((math.pow(w0,2) - math.pow(w,2))**2 + \
        math.pow(w*gamma, 2))
    return [eps_r, eps_i]

#================OPTIMIZATION PART=========================
# doing only for red lines (through plane)
# Array of energies real and imaginary
w_ri = np.array([
    [0.60225785,0.60225624],
    [0.6812816,0.6812236],
    [0.74196064,0.7418592],
    [0.77723914,0.7728821],
    [0.81675106,0.81659615],
    [0.86190754,0.8617204],
    [0.88307464,0.88428247],
    [0.90141946,0.90684456],
    [0.9098863,0.9153054],
    [0.92258656,0.9279965],
    [0.9310534,0.93081677],
    [0.93952024,0.9392776],
    [0.9550427,0.95478904],
    [0.97338754,0.97312075],
    [1.0213664,1.0252956],
    [1.0538225,1.0549084],
    [1.0961567,1.0915718],
    [1.1540134,1.1536175],
    [1.2231593,1.2255342],
    [1.3021832,1.3030914],
    [1.3924961,1.3919297],
    [1.5053872,1.5005099],
    [1.6366233,1.6330621],
    [1.7523367,1.750103],
    [1.8680501,1.8657339],
    [3.1804101,3.1799762],
])

eps_ri = np.array([
    [6.3534193,0.92022264],
    [7.098172,1.5139147],
    [8.181449,2.4341373],
    [9.197021,3.4730983],
    [10.1787405,5.432282],
    [10.212593,10.211503],
    [8.316858,13.239332],
    [6.116452,15.584415],
    [4.0176034,15.2282],
    [1.8849019,14.1595545],
    [-0.3493568,13.387755],
    [-2.041977,12.734694],
    [-3.396073,11.2801485],
    [-4.5132027,9.0834875],
    [-4.242383,4.6011133],
    [-3.3622208,3.5324676],
    [-2.380501,2.4935064],
    [-1.5680434,1.7217069],
    [-0.7217333,1.3061224],
    [-0.14624238,1.0983303],
    [0.19228165,1.0092764],
    [0.46310088,0.83116883],
    [0.6323629,0.7717996],
    [0.7339201,0.742115],
    [0.7677725,0.7124304],
    [1.140149,0.3265306],
])

new_eps_ri = np.ones((26,2), dtype=float)
# x: array([ 0.        ,  2.02942628,  0.4       ])

#~ FF = 0.0
#~ w0 = 0.9310534
#~ #x: array([ 0.9440551 ,  2.22047041,  0.11243039])
#~ wp = 0.9440551
#~ eps_inf = 2.22047041
#~ gamma = 0.11243039

#~ for i in range(len(w_ri)):
    #~ eps_r = (math.pow(wp,2)*(math.pow(w0,2) - math.pow(w_ri[i][0],2))) / \
        #~ ((math.pow(w0,2) - math.pow(w_ri[i][0],2))**2 + math.pow(w_ri[i][0]*gamma, 2))+eps_inf

    #~ eps_i = math.pow(wp,2)*gamma*w_ri[i][1]/((math.pow(w0,2) - math.pow(w_ri[i][1],2))**2 + \
        #~ math.pow(w_ri[i][1]*gamma, 2))

    #~ new_eps_ri[i][0] = eps_r
    #~ new_eps_ri[i][1] = eps_i

#~ for i in range(len(w_ri)):
    #~ FF += math.pow((eps_ri[i][0] - new_eps_ri[i][0]),2) + math.pow((eps_ri[i][1] - new_eps_ri[i][1]),2)

#~ print "Fitness function value: (through plane, init_guess: wp=",wp,"epsinf=","gamma=", gamma,"):", FF

# comment the previous loops
def objective(x):
    w0 = 0.9310534  #through plane, red line
    FF = 0.0

    wp = x[0]
    eps_inf = x[1]
    gamma = x[2]
    for i in range(len(w_ri)):
        eps_r = float((math.pow(wp,2)*(math.pow(w0,2) - math.pow(w_ri[i][0],2)))) / \
            float(((math.pow(w0,2) - math.pow(w_ri[i][0],2))**2 + math.pow(w_ri[i][0]*gamma, 2)))+eps_inf

        eps_i = float(math.pow(wp,2)*gamma*w_ri[i][1])/float(((math.pow(w0,2) - math.pow(w_ri[i][1],2))**2 + \
            math.pow(w_ri[i][1]*gamma, 2)))

        new_eps_ri[i][0] = eps_r
        new_eps_ri[i][1] = eps_i

    for i in range(len(w_ri)):
        FF += math.pow((eps_ri[i][0] - new_eps_ri[i][0]),2) + math.pow((eps_ri[i][1] - new_eps_ri[i][1]),2)

    return FF

#Initial guess obtained manually through Excel
x0 = [0, 0, 0.3]
# Boundaries for wp, eps_inf and gamma
b1 = (0, 3)
b2 = (0, 0.1)
bnds = (b1,b1,b2)
#sol = opt.minimize(objective, x0, bounds=bnds)
#print sol
#x0 = 0, 0, 0.3
#fun 150
#: array([ 1.21041352,  1.58780666,  0.1       ])
