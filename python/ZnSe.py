from __future__ import division
import numpy as np
import cmath

#Array of wavelength[nm] and epsilon.real
# [486...1098]nm
wl_eps_real = np.array([
    [486.78458, 0.99541223],
    [487.88348, 1.7356131],
    [489.5089, 2.2385597],
    [490.06662, 2.8269281],
    [490.62323, 3.3868268],
    [491.71103, 3.8328402],
    [492.5311, 4.2788568],
    [493.7464, 4.5564227],
    [494.5582, 4.784171],
    [497.3759, 4.954956],
    [499.39362, 5.21116],
    [501.00757, 5.4104295],
    [502.22153, 5.652408],
    [503.03735, 5.986918],
    [503.8508, 6.257371],
    [504.66663, 6.5918813],
    [505.47604, 6.755573],
    [508.69586, 6.940588],
    [514.7119, 6.7269936],
    [519.94324, 7.0045114],
    [524.7494, 6.6556997],
    [529.9837, 7.01151],
    [536.00806, 7.0185556],
    [538.0102, 6.8619485],
    [538.80664, 6.684003],
    [540.0034, 6.4704657],
    [541.2007, 6.271163],
    [543.1961, 5.9366198],
    [546.4192, 6.207044],
    [546.8283, 6.4063277],
    [548.0404, 6.598484],
    [548.8487, 6.733706],
    [552.8747, 6.9970036],
    [557.28564, 6.819015],
    [558.0846, 6.705127],
    [559.2824, 6.520059],
    [560.4789, 6.299404],
    [560.87244, 6.085876],
    [561.66675, 5.8509912],
    [564.8651, 5.4666114],
    [566.8817, 5.6943455],
    [567.6941, 5.936329],
    [568.5054, 6.1498427],
    [570.11957, 6.3562293],
    [570.9306, 6.5626254],
    [572.1416, 6.726312],
    [575.7641, 6.939792],
    [579.77374, 6.768926],
    [581.7702, 6.4628525],
    [582.56537, 6.24932],
    [584.1645, 6.05713],
    [584.95886, 5.822245],
    [585.35077, 5.5660124],
    [586.14703, 5.3809495],
    [589.3526, 5.1887407],
    [591.366, 5.3310657],
    [592.1773, 5.544579],
    [592.5867, 5.7509804],
    [594.2017, 5.978719],
    [595.4151, 6.206463],
    [596.6275, 6.405737],
    [597.03796, 6.6406074],
    [601.86584, 6.8683085],
    [606.6746, 6.590671],
    [607.873, 6.4198384],
    [609.0708, 6.234771],
    [609.8681, 6.0781775],
    [610.6622, 5.836175],
    [612.26056, 5.622633],
    [613.0562, 5.423335],
    [613.8516, 5.21692],
    [616.2511, 4.946429],
    [619.0701, 5.152801],
    [620.2825, 5.352075],
    [621.09375, 5.5655885],
    [621.9053, 5.7862196],
    [623.1184, 6.006846],
    [626.3501, 6.5050282],
    [630.37744, 6.803913],
    [635.59216, 6.64015],
    [637.5886, 6.3340764],
    [639.98535, 5.992411],
    [642.37836, 5.551101],
    [644.3729, 5.1952057],
    [646.7691, 4.839305],
    [651.19727, 5.1168327],
    [652.4088, 5.2947545],
    [652.4158, 5.479808],
    [654.841, 5.8925905],
    [656.46063, 6.241326],
    [658.88074, 6.5188775],
    [662.90674, 6.782175],
    [668.12335, 6.668234],
    [671.3276, 6.440438],
    [674.123, 6.020476],
    [675.7135, 5.600528],
    [677.70703, 5.216162],
    [680.10406, 4.881614],
    [682.90564, 4.6253533],
    [686.93353, 4.9384727],
    [688.9539, 5.265851],
    [690.57513, 5.657291],
    [693.403, 6.098539],
    [695.4207, 6.354743],
    [701.4569, 6.674956],
    [707.07404, 6.532541],
    [711.4726, 6.02715],
    [715.06744, 5.5075345],
    [717.4577, 4.99505],
    [719.45447, 4.696094],
    [722.65765, 4.439828],
    [727.0828, 4.639064],
    [729.90967, 5.051842],
    [732.73944, 5.542912],
    [735.5685, 6.0126295],
    [737.98615, 6.226124],
    [742.41644, 6.560591],
    [746.836, 6.6103606],
    [752.0496, 6.418128],
    [756.04846, 5.9625645],
    [759.2455, 5.5425973],
    [760.84033, 5.2365284],
    [762.43384, 4.8948727],
    [764.02954, 4.610156],
    [768.4346, 4.2755837],
    [774.46515, 4.446331],
    [777.6973, 4.9587483],
    [780.52606, 5.4213486],
    [782.5448, 5.7060223],
    [786.5773, 6.140138],
    [788.99603, 6.3821025],
    [795.0274, 6.574202],
    [800.2437, 6.4531436],
    [804.25226, 6.2538075],
    [809.053, 5.762647],
    [811.8484, 5.342685],
    [814.24054, 4.880023],
    [817.8402, 4.488521],
    [820.2383, 4.1824427],
    [825.05164, 4.025802],
    [831.49695, 4.5452986],
    [833.91736, 4.8299675],
    [837.5493, 5.292558],
    [839.96967, 5.577227],
    [843.196, 5.9330606],
    [848.83325, 6.3244534],
    [854.05817, 6.431153],
    [861.2793, 6.2246614],
    [867.6921, 5.8829484],
    [872.08875, 5.3277354],
    [875.282, 4.8081245],
    [880.88135, 4.195958],
    [882.0787, 3.9966557],
    [888.0952, 3.797296],
    [892.11664, 3.9395971],
    [897.3566, 4.4448733],
    [901.39044, 4.9145765],
    [904.6178, 5.29888],
    [910.2623, 5.882443],
    [915.49475, 6.1884313],
    [922.73065, 6.373399],
    [930.75714, 6.2238374],
    [936.7697, 5.917716],
    [942.7747, 5.4123063],
    [947.573, 4.857089],
    [951.9732, 4.3944035],
    [955.97394, 3.9886615],
    [963.18726, 3.5757644],
    [971.63336, 3.9030666],
    [977.67737, 4.4296856],
    [984.53125, 5.134231],
    [988.1594, 5.497177],
    [994.1994, 5.9170346],
    [1000.2329, 6.1660733],
    [1007.46606, 6.2798667],
    [1015.4923, 6.123188],
    [1022.71045, 5.8384047],
    [1027.9152, 5.4112964],
    [1033.5194, 4.927244],
    [1037.5175, 4.450328],
    [1044.726, 3.909317],
    [1047.13, 3.7598221],
    [1053.9419, 3.3540473],
    [1062.783, 3.5034087],
    [1070.4325, 4.008656],
    [1078.0862, 4.6206656],
    [1084.9309, 5.0832176],
    [1091.7749, 5.524418],
    [1098.2162, 5.9371533],
])

#Array of wavelength[nm] and epsilon.imag
# [487...1097]nm
wl_eps_imag = np.array([
    [487.5591, 0.23620936],
    [490.37466, 0.35005504],
    [495.99744, 0.35710576],
    [502.42096, 0.30009004],
    [508.04105, 0.2359664],
    [514.47, 0.32129943],
    [523.2997, 0.17172843],
    [532.13873, 0.27126774],
    [542.57544, 0.13591257],
    [552.2185, 0.25679466],
    [563.0563, 0.107199855],
    [575.5125, 0.27787068],
    [587.55396, 0.09979181],
    [600.8126, 0.24910079],
    [614.4613, 0.092355184],
    [630.1307, 0.27010533],
    [646.1885, 0.09909627],
    [663.464, 0.26970991],
    [681.93146, 0.09867227],
    [701.6169, 0.2763748],
    [722.8942, 0.06971659],
    [744.9906, 0.2829777],
    [769.8824, 0.07627664],
    [795.5938, 0.30372974],
    [823.2955, 0.061408147],
    [853.4259, 0.32439604],
    [887.95386, 0.053523704],
    [922.5022, 0.32357663],
    [963.0542, 0.05263283],
    [1004.8325, 0.35106975],
    [1053.4153, 0.04444349],
    [1097.1993, 0.27879965],
])

def eps_ZnSe_Marple(x1):
    # Calculates epsilon of ZnSe when input energy in [nm]
    # source: https://refractiveindex.info/?shelf=main&book=ZnSe&page=Marple
    # Use:
    # print eps_ZnSe_Marple(500)
    # [7.46715328467, 0] - 0 stands for imaginary part
    x = x1*1E-3
    eps = 4.0+(1.90*x**2)/(x**2 - 0.113)
    return [eps, 0]

def eps_ZnSe_Desai(x):
    # fn interates through known values of eps(eV) and interpolates
    # source Desai:"Optical and Dispersion Analysis of ZnSe thing film"
    # ZnSe thickness = 7300 A (730 nm)

    # How to use:
    #x = 492.5 #nm
    #a, b = eps_ZnSe(x)
    #print "ZnSe epsilon(real, imag) for",x, "[nm] is:", a,b

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
