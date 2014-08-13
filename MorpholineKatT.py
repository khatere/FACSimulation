import numpy as np
import matplotlib.pyplot as plt
import KeqatT
import math
from scipy.optimize import fsolve
deg = 5                      # Order of the fitted polynomial
T= 25                      # The temperature at which the K is needed to be calculated
Kmorpholine = KeqatT.logKFitT(deg,T)
Kw = KeqatT.KwatT(T)
KFeOH3 = KeqatT.KFeOH3atT(T)
KFeOH2 = KeqatT.KFeOH2atT(T)
KFeOH = KeqatT.KFeOHatT(T)
DebyeHuckConst = KeqatT.DebyeHuckelConst(T)
print  Kmorpholine, Kw, DebyeHuckConst, KFeOH3, KFeOH2, KFeOH

#30-60 mg/kg of added Morpholine (morpholine molar mass: 87.1 g/mol) equals to  0.000344-0.000688 mol/L take an average value of 0.000516 mol/L
ConcC4H9ONTotal = 0.000516 #mol/L

#solve the following equations for ConcH, ConcOH, ConcC4H9ONH, IStrength, gamma1
# Kw = ConcH*gamma1*ConcOH*gamma1;   
# ConH + ConcC4H9ONH = ConOH;  charge neutrality  
# Kmorpholine = ConcC4H9ONH*gamma1*ConOH*gamma1/ConcC4H9ON;
# ConcC4H9ONTotal = ConcC4H9ONH + ConcC4H9ON = ConcC4H9ONH + gamma1^2*ConcC4H9ONH*ConOH/Kmorpholine
# IStrength = ((z^2)*ConcH + (z^2)*ConcOH +(z^2)*ConcC4H9ONH)/2; z=1;
# log(gamma1) =  -DebyeHuckconst*z^2*[(sqrt(I)/(1+sqrt(I)))-Beta*I];   (Davis equation)         
#Since the concentrations are very small, I multiply them by 1e5 and introduced new concentrations and IStrength.
def equations(p):
    out = [(p[0]/100000) - (p[1]/100000) + (p[2]/100000)]
    out.append(((p[0]/100000) + (p[1]/100000) + (p[2]/100000))/2 - (p[3]/100000))
    out.append((p[0]/100000)*(p[1]/100000)*p[4]*p[4] - Kw)
    out.append((p[2]/100000) + (p[2]/100000)*(p[1]/100000)*p[4]*p[4]/Kmorpholine - ConcC4H9ONTotal)
    out.append(math.pow(10,(-1*DebyeHuckConst*(math.sqrt(p[3]/100000)/(1 + math.sqrt(p[3]/100000)) - 0.2*p[3]/100000))) - p[4])
    return out
P1 =  fsolve(equations, [10, 10, 10, 10, 10]) # fsolve is sensitive to initial guess
[ConcH, ConcOH, ConcC4H9ONH, IStrength, gamma1] = [ P1[0]*1e-5, P1[1]*1e-5, P1[2]*1e-5, P1[3], P1[4]] 
print 'ConcH=', ConcH, 'ConcOH=', ConcOH, 'ConcC4H9ONH=', ConcC4H9ONH, 'IStrength=', IStrength, 'gamma1=', gamma1
pH = -1*math.log(ConcH,10)
print 'pH=', pH



