from __future__ import division
import numpy as np
import sys
from math import *
import cmath
from TMM_EMT_TI import calculateRpTrAb
#import postprocess


#1. Read configurations ('Bi2Se3.configs' or 'Bi2Te3.configs')
cfg = np.array(pl.loadtxt('Bi2Se3.configs',skiprows=3))
#=============================
#N| d_Bulk| d_Substr| d_Total|
#=============================
#1     100     100     2000
#2     30     50     2000
#3     100     10     2000
#4     12     10     2000

#2. Setup which configurations to run and what resolution to use
for i in range(len(cfg[0])):
    calculateRpTrAb(cfg[i][1], cfg[i][2], cfg[i][3], 10, 18)
    # (substrate, ti, total, wl_reso, angle_reso - optional)