#Calculating the dissociation constants of water and morpholine  and hydrolysis constants of iron species at the working temperature
import numpy as np
import matplotlib.pyplot as plt
import math
                      
                

# Fitting a polynomial on the data of Nesmer and Hitch (0-300C)(J. Sol. Chemistry, vol. 6, No. 4, 1977)
# to calculate the hydrolysis constant of morpholine at the working temperature 
# C4H9ON+H2O=C4H9ONH+ +OH- ,K=[C4H9ONH+]*[OH-]/[C4H9ON]
# deg = Order of the fitted polynomial
# T = The temperature at which the K is needed to be calculated (degree C)
def logKFitT(deg,T):    
 temp = [0.0, 25., 50., 75., 100., 125., 150., 175., 200., 225., 250., 275., 300.] #(Degree C)
 logK = [-5.837, -5.505, -5.309, -5.208, -5.178, -5.204, -5.278, -5.394, -5.549, -5.740, -5.965, -6.222, -6.512]
 fitCoeff = np.polyfit(temp, logK, deg)
 fitPoly = np.poly1d(fitCoeff)
 Kmorph = math.pow(10,fitPoly(T))
 return Kmorph

# Calculating Kw at Temperature T
def KwatT(T):
    Tkelvin = T + 273.15 
    logKwConst = -55.86334915 + 0.3276529697411*Tkelvin + (-0.0009593113681)*math.pow(Tkelvin,2) + 0.0000013266059*math.pow(Tkelvin,3) + (-0.0000000007226)*math.pow(Tkelvin,4)
    Kw = math.pow(10,logKwConst)
    return Kw

# Calculating hydrolysis constant of Fe(OH)3 at Temperature T
def KFeOH3atT(T):
    Tkelvin = T + 273.15 
    logKFeOH3Const = -97.4709 + 0.413186*Tkelvin + (-0.000935775)*math.pow(Tkelvin,2) + 0.0000010496*math.pow(Tkelvin,3) + (-0.0000000004667)*math.pow(Tkelvin,4) 
    KFeOH3 = math.pow(10,logKFeOH3Const)
    return KFeOH3

# Calculating hydrolysis constant of Fe(OH)2 at Temperature T 
def KFeOH2atT(T):
    Tkelvin = T + 273.15 
    logKFeOH2Const = -67.8248 + 0.311499*Tkelvin + (-0.00073731)*math.pow(Tkelvin,2) + 0.000000863467*math.pow(Tkelvin,3) + (-0.0000000004)*math.pow(Tkelvin,4)
    KFeOH2 = math.pow(10,logKFeOH2Const)
    return KFeOH2

# Calculating hydrolysis constant of FeOH at Temperature T
def KFeOHatT(T):
    Tkelvin = T + 273.15 
    logKFeOHConst = -41.5545 + 0.230763*Tkelvin + (-0.00062963)*math.pow(Tkelvin,2) + 0.00000081013*math.pow(Tkelvin,3) + (-0.0000000004)*math.pow(Tkelvin,4)
    KFeOH = math.pow(10,logKFeOHConst)
    return KFeOH
 
 
 
def DebyeHuckelConst(T):
    DHConst = 0.5027+ (-0.0009028)*T + 0.00003315*math.pow(T,2) + (-0.0000001709)*math.pow(T,3) + 0.000000000329*math.pow(T,4) 
    return DHConst
