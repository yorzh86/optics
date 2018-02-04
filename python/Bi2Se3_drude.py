#todo:
# check bulk eps
# do new structures + periods.

#concerns:
# interpolation will give WAY better results than drude model

import numpy as np
import scipy.optimize as opt
import math

#Black - parallel (extraordinary)
#Red   - through plane (ordinary)

#for red line (ordinary, through plane) works for [0.6...4.2 eV]
def drude_O_eps(w1, w0=1.0585527, wp=2.63263661, eps_inf=2.55029383, gamma=0.13258352):
    #Calculates epsilon, given energy[nm]. Other parameters w0,wp,etc can also be changed
    #a = drude_O_eps(500, w0=1.0585527, eps_inf=1.0)
    #print a
    
    #convert [nm] to [eV] #WIKI correct!!
    w = 1.2398/w1*1E3
    eps_r = (math.pow(wp,2)*(math.pow(w0,2) - math.pow(w,2))) / \
        ((math.pow(w0,2) - math.pow(w,2))**2 + math.pow(w*gamma, 2))+eps_inf
    eps_i = math.pow(wp,2)*gamma*w/((math.pow(w0,2) - math.pow(w,2))**2 + \
        math.pow(w*gamma, 2))
    return [eps_r, eps_i] 

#for black line (extra-ordinary, in plane) works for [0.6...4.2 eV]
def drude_E_eps(w1, w0=1.7211499, wp=8.46042266, eps_inf=0.18658023, gamma=0.98665155):
    #Calculates epsilon, given energy[nm]. Other parameters w0,wp,etc can also be changed
    #a = drude_E_eps(500, eps_inf=1.0, gamma = 0.05)
    #print a
    
    #convert [nm] to [eV] #WIKI correct!!
    w = 1.2398/w1*1E3
    eps_r = (math.pow(wp,2)*(math.pow(w0,2) - math.pow(w,2))) / \
        ((math.pow(w0,2) - math.pow(w,2))**2 + math.pow(w*gamma, 2))+eps_inf
    eps_i = math.pow(wp,2)*gamma*w/((math.pow(w0,2) - math.pow(w,2))**2 + \
        math.pow(w*gamma, 2))
    return [eps_r, eps_i]

# ========Optimization part:===========
#Black - in plane (extraordinary)
#Red   - through plane (ordinary)

# Array grabbed from plot for red line Bi2Se3 (real and imaginary)
O_w_ri = np.array([
    [0.62624466, 0.6247448],
    [0.7771068, 0.77266574],
    [0.86592436, 0.86809313],
    [0.932381, 0.93081313],
    [0.9519844, 0.9573789],
    [0.9767079, 0.977389],
    [0.9892141, 0.9875487],
    [1.0069029, 1.0043309],
    [1.0277587, 1.021432],
    [1.0539126, 1.0585527],
    [1.0722132, 1.0748636],
    [1.0908935, 1.0908935],
    [1.1022494, 1.1035907],
    [1.1144665, 1.116231],
    [1.1342001, 1.1353745],
    [1.1514072, 1.1546656],
    [1.1809858, 1.1806554],
    [1.2080476, 1.206707],
    [1.2475146, 1.2426565],
    [1.3018408, 1.3049433],
    [1.3982708, 1.3968095],
    [1.5245011, 1.5249014],
    [1.6458377, 1.646428],
    [1.7746327, 1.771249],
    [2.1314154, 2.1358812],
    [2.3246994, 2.3264074],
    [2.5105395, 2.5103636],
    [2.706297, 2.7041655],
    [2.8946295, 2.891416],
    [3.023487, 3.026097],
    [3.191995, 3.1969132],
    [3.3629854, 3.3644445],
    [3.5141416, 3.5122654],
    [3.6578727, 3.6502314],
    [3.7966428, 3.7947674],
    [3.939313, 3.939313],
    [4.123746, 4.126554],
    ])

# Array grabbed from plot for black line Bi2Se3 (real and imaginary)
E_w_ri = np.array([
    [0.6037759, 0.60985625],
    [0.7021417, 0.7010267],
    [0.78033704, 0.7823409],
    [0.86105645, 0.861191],
    [0.93672764, 0.93511295],
    [1.0123988, 1.0114989],
    [1.0805163, 1.0829569],
    [1.1713098, 1.1765914],
    [1.2519844, 1.2505133],
    [1.3376362, 1.3342916],
    [1.3904895, 1.3958932],
    [1.4332349, 1.4328542],
    [1.4709061, 1.4698153],
    [1.5110716, 1.5141684],
    [1.5461738, 1.5437372],
    [1.5862794, 1.5856262],
    [1.6388748, 1.6324435],
    [1.681456, 1.684189],
    [1.7265613, 1.7211499],
    [1.769191, 1.7605749],
    [1.8193818, 1.8197125],
    [1.8796655, 1.8788501],
    [1.9400238, 1.9453799],
    [2.0180733, 2.0119097],
    [2.1087883, 2.1087883],
    [2.1995332, 2.1967146],
    [2.2928433, 2.292813],
    [2.4063277, 2.4086242],
    [2.5122435, 2.517043],
    [2.6282556, 2.620534],
    [2.7417326, 2.7412732],
    [2.8602426, 2.8620124],
    [2.9913511, 2.9901438],
    [3.1149058, 3.110883],
    [3.2636492, 3.263655],
    [3.3846612, 3.3893223],
    [3.5233195, 3.524846],
    [3.659446, 3.6579056],
    [3.793056, 3.7958932],
    [3.916581, 3.9166324],
    [4.0350575, 4.034908],
    [4.1484933, 4.1482544],
    ])

# a loop through energy array[eV] to get epsE and epsO
#~ for i in range(len(E_w_ri)):
    #~ print drude_E_eps(E_w_ri[i][0])
#~ for i in range(len(O_w_ri)):
    #~ print drude_O_eps(O_w_ri[i][0])
    
O_old_eps_ri = np.array([
    [12.526751, 1.318346],
    [16.254206,2.6698534],
    [21.100832, 4.972083],
    [26.655085, 9.309189],
    [29.372524, 13.376189],
    [30.08233, 17.646776],
    [28.665588, 21.985453],
    [24.4735, 27.069695],
    [11.836687, 36.6963],
    [-0.315187454, 50.72908],
    [-16.684837, 49.101475],
    [-19.34486, 43.47388],
    [-20.34486, 37.16842],
    [-18.218311, 30.049398],
    [-17.095352, 21.981031],
    [-15.381959, 16.014355],
    [-13.431773, 11.877986],
    [-11.009263, 8.622971],
    [-8.763343, 5.9778285],
    [-6.3986254, 4.1454577],
    [-3.7366822, 2.5155919],
    [-1.8410742, 2.1727777],
    [-0.71340847, 1.8979565],
    [0.11932759, 1.7586299],
    [0.9626235, 1.6799271],
    [1.089686, 1.6064318],
    [1.3345139, 1.533133],
    [1.520746, 1.3239466],
    [1.5885242, 1.3861425],
    [1.6535476, 1.3143176],
    [1.7204076, 1.2414118],
    [1.7283274, 1.1686044],
    [1.8534387, 1.0963864],
    [1.860096, 1.0244632],
    [1.9255786, 0.9523435],
    [2.0588393, 1.0158167],
    [2.0588393, 0.9424197],
])

E_old_eps_ri = np.array([
    [20.430677, 17.809725],
    [21.57374, 16.490486],
    [22.596828, 16.186047],
    [23.680172, 16.186047],
    [24.643005, 16.693445],
    [25.605839, 17.657505],
    [26.68963, 19.128963],
    [27.651926, 21.767443],
    [28.011131, 24.81184],
    [27.344297, 29.073996],
    [26.13552, 33.08245],
    [24.50469, 35.41649],
    [22.330935, 37.80127],
    [19.734676, 40.490486],
    [16.776527, 42.012684],
    [13.21475, 43.890064],
    [7.842182, 44.9556],
    [3.5561783, 45.412262],
    [-0.66957, 44.60042],
    [-4.171091, 43.078224],
    [-7.67288, 40.33827],
    [-10.993993, 36.53277],
    [-13.10821, 32.372093],
    [-14.438569, 28.515856],
    [-14.743514, 23.7812433],
    [-14.5657, 19.484144],
    [-13.663837, 16.389006],
    [-12.581656, 13.547568],
    [-11.619897, 11.51797],
    [-10.417116, 9.894292],
    [-9.455625, 8.524313],
    [-8.615003, 7.408034],
    [-7.8351717, 6.5961943],
    [-6.934383, 5.835095],
    [-6.2758684, 5.1754756],
    [-5.7370596, 4.566596],
    [-5.1385317, 4.1606765],
    [-4.720949, 3.653277],
    [-4.2429323, 3.4503171],
    [-3.8249023, 3.1966174],
    [-3.527383, 3.0443974],
    [-3.2296848, 2.9429176],
])

w0O = 1.0585527   #ordinary, red line
w0E = 1.7211499   #extraordinary, black line
FFE = 0.0
FFO = 0.0
E_new_eps_ri = np.zeros((42,2), dtype=float)
O_new_eps_ri = np.zeros((37,2), dtype=float)

#initial values to test ff_initial or for optimization 1st guess
wp =  2.25
eps_inf = 2.72
gamma = 0.082

# IMPORTANT!!! OPTIMIZATION PROCEDURE IS IN EV!!!
    
for i in range(len(O_w_ri)):
    eps_r = (math.pow(wp,2)*(math.pow(w0O,2) - math.pow(O_w_ri[i][0],2))) / \
        ((math.pow(w0O,2) - math.pow(O_w_ri[i][0],2))**2 + math.pow(O_w_ri[i][0]*gamma, 2))+eps_inf
    eps_i = math.pow(wp,2)*gamma*O_w_ri[i][1]/((math.pow(w0O,2) - math.pow(O_w_ri[i][1],2))**2 + \
        math.pow(O_w_ri[i][1]*gamma, 2))
    O_new_eps_ri[i][0] = eps_r
    O_new_eps_ri[i][1] = eps_i

for i in range(len(E_w_ri)):
    eps_r = (math.pow(wp,2)*(math.pow(w0E,2) - math.pow(E_w_ri[i][0],2))) / \
        ((math.pow(w0E,2) - math.pow(E_w_ri[i][0],2))**2 + math.pow(E_w_ri[i][0]*gamma, 2))+eps_inf
    eps_i = math.pow(wp,2)*gamma*E_w_ri[i][1]/((math.pow(w0E,2) - math.pow(E_w_ri[i][1],2))**2 + \
        math.pow(E_w_ri[i][1]*gamma, 2))
    E_new_eps_ri[i][0] = eps_r
    E_new_eps_ri[i][1] = eps_i

for i in range(len(O_w_ri)):
    FFO += math.pow((O_old_eps_ri[i][0] - O_new_eps_ri[i][0]),2) + math.pow((O_old_eps_ri[i][1] - O_new_eps_ri[i][1]),2)

for i in range(len(E_w_ri)):
    FFE += math.pow((E_old_eps_ri[i][0] - E_new_eps_ri[i][0]),2) + math.pow((E_old_eps_ri[i][1] - E_new_eps_ri[i][1]),2)

#print "Fitness function value: (ordinary, init_guess: wp=2.25, eps_inf=2.72, gamma=0.082):", FFO #1570.7874974
#print "Fitness function value: (extra-ordinary, init_guess: wp = 2.25, eps_inf=2.72, gamma=0.082):", FFE #28652.518502

def objective(x):
    w0 = 1.0585527  #ordinary, red line
    #w0 = 1.7211499  #extraordinary, black line 
    FF = 0.0
    
    wp = x[0]
    eps_inf = x[1]
    gamma = x[2]
        
    for i in range(len(O_w_ri)):
        eps_r = (math.pow(wp,2)*(math.pow(w0,2) - math.pow(O_w_ri[i][0],2))) / \
            ((math.pow(w0,2) - math.pow(O_w_ri[i][0],2))**2 + math.pow(O_w_ri[i][0]*gamma, 2))+eps_inf
        eps_i = math.pow(wp,2)*gamma*O_w_ri[i][1]/((math.pow(w0,2) - math.pow(O_w_ri[i][1],2))**2 + \
            math.pow(O_w_ri[i][1]*gamma, 2))
        O_new_eps_ri[i][0] = eps_r
        O_new_eps_ri[i][1] = eps_i
    
    for i in range(len(O_w_ri)):
        FF += math.pow((O_old_eps_ri[i][0] - O_new_eps_ri[i][0]),2) + math.pow((O_old_eps_ri[i][1] - O_new_eps_ri[i][1]),2)
    
    return FF

#Initial guess obtained manually through Excel
x0 = [2.25, 1.0, 0.082]
# Boundaries for wp, eps_inf and gamma
b1 = (0.0, 5)
b2 = (0.0, 0.3)
bnds = (b1,b1,b2)
sol = opt.minimize(objective, x0, bounds=bnds)

# Red (ordinary, through plane)
# initial,                                 fn = 1570, sol = [2.25, 2.72, 0.082]
# best, bounds(0:5, 0:5, 0:0.3)            fn = 468,  sol = [2.63263661, 2.55029383, 0.13258352]
# compromise?, bounds(0:5, 0:5, 0.03:0.09) fn = 1190, sol = [2.32592269, 2.68246867, 0.09]

# Black (extraordinary, in plane)
# initial,                                 fn = 28652, sol = [2.25, 2.72, 0.082]
# best, bounds(0:10, 0:10, 0:1.0)          fn = 613,  sol = [8.46042266,  0.18658023,  0.98665155]]
