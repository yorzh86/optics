from __future__ import division
import numpy as np
import sys
import math

def d_conduct():
    a = 1.9*1E-9
    return a

def material_name():
    a = "Bi2Te3-ZnSe"
    return a

def muf():
    muf = 0.194 #Yin table 1
    return muf

def tau():
    tau = 54*1E-15 #Yin table 1
    return tau