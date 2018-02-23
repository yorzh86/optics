from __future__ import division
import numpy as np
import sys
import math

def d_conduct():
    a = 0.92*1E-9
    return a

def material_name():
    a = "Bi2Se3-ZnSe"
    return a

def muf():
    muf = 0.189 #Yin table 1
    return muf

def tau():
    tau = 100*1E-15 #Yin table 1 55-155
    return tau

def eps_dielectric():
    eps = math.pow(2.5, 2)
    #eps = math.sqrt(2.5)
    return [eps,0]