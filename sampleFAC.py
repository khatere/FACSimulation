import os
import csv
import random as rand
import matplotlib.pyplot as plt
import numpy as np
import math
import ConcentrationCal
import constantssample
import ActivationEnergyCal

# # Assigning the constants and initial values
# values = constantssample.Const()
# Temperature = values.Temperature()
# TKelvin = values.Temperature() + 273.15
# ConcC4H9ONTotal = values.ConcC4H9ONTotal()
# CSat_MP = values.CSat_MP()
# CSat_OB = values.CSat_OB()
# CSat_PO = values.CSat_PO()
# CB = values.CB()
# kd = values.kd()
# t = values.t()
# ksp = values.ksp()
# FF = values.FF()
# CrContent = values.CrContent()
# U = values.U()
# Ce_OB = values.Ce_OB() 
# Ce_MP = 5.0e-8
# CeH2_MP = 8.0e-7 # gr/cm^3; if given in cm^3/kg: 10.0 cm^3/kgH2O *P/(RT)*Molar mass H2*ro water => 10*101325*2.016*0.926/(8.1314*298.15*1e6*1e3)
# CeH2_OB = 8.0e-7
# CSat_Cr = 1.0e-10
# Delta = values.DOXs()
# Thickness = list()
# corrosionnrate = list()
# FeCr2O4 = 0
# DummyTime = t
# Time = list()
# j = 0 
# spall = list()
# mu, sigma = 0, 10 # mean and standard deviation
# # ActivationEnergyValues = ActivationEnergyCal.ActivationEnergy(Ce_MP,CeH2_MP, Ce_OB, CeH2_OB, Temperature, ConcC4H9ONTotal)
# # ActivationEFe_Fe2_MP = ActivationEnergyValues.getActivationEnergyFe_Fe2_MP()
# # ActivationEH2_Hp_MP = ActivationEnergyValues.getActivationEnergyH2_Hp_MP()
# # ActivationEH2_Hp_OB = ActivationEnergyValues.getActivationEnergyH2_Hp_OB()
# # ActivationEFe2p_Fe3O4_OB = ActivationEnergyValues.getActivationEnergyFe2p_Fe3O4_OB()
# ActivationEFe_Fe2_MP = 264385.58
# ActivationEH2_Hp_MP = 100502.6341
# ActivationEH2_Hp_OB = 377712.8108
# ActivationEFe2p_Fe3O4_OB = 291642.3046
# 
# print  ActivationEFe_Fe2_MP,ActivationEH2_Hp_MP, ActivationEH2_Hp_OB, ActivationEFe2p_Fe3O4_OB 
# 
# 
# #corrosion rate equation 
# def DDeltaFunction(Delta, CSat_MP, CSat_PO, CSat_OB, CB, Ce_OB, kd):
#     A = values.S() * (CSat_MP + CSat_PO)
#     B = ((kd * values.f() * values.S() * CSat_OB) + values.km() * CB) / (kd * values.f() + values.km())
#     C = 0.0#(values.DeltaP() * values.ChiP())/(values.PhiP() * values.DP()* values.RhoP() * (1 - values.PhiP()))
#     D = 0.0#(values.MFe() * values.DeltaP() * values.ChiP()) / (values.MP() * values.RhoFe() * values.PhiP() * values.DP()) # Why DFe, should not be DP????????????????????
#     E = (0.476 * (1.101 + values.PhiOX()) * values.ChiOX() * Delta) /(values.PhiOX() * values.DFe() * values.RhoOX() * (1 - values.PhiOX()))
#     F = (0.476 * (1.101 + values.PhiOX())) / (kd * values.f() + values.km())
#     DM = (A - B) / (C - D + E + F)
#     DDelta = (0.476 * DM * (1 - values.PhiOX()) - (kd * values.f() * (values.S() * CSat_OB - Ce_OB))) / 0.723 # variation in oxide thickness
#     return  DDelta
#  
#  
# def CorrosionRateFunction(Delta, CSat_MP, CSat_PO, CSat_OB, CB):
#      AA = values.S() * (CSat_MP + CSat_PO)
#      BB = ((kd * values.f() * values.S() * CSat_OB) + values.km() * CB) / (kd * values.f() + values.km())
#      CC = 0.0#(values.DeltaP() * values.ChiP())/(values.PhiP() * values.DP() * values.RhoP() * (1 - values.PhiP()))
#      DD = 0.0#(values.MFe() * values.DeltaP() * values.ChiP()) / (values.MP() * values.RhoFe() * values.PhiP() * values.DP()) # Why DFe, should not be DP????????????????????
#      EE = (0.476 * (1.101 + values.PhiOX()) * values.ChiOX() * Delta) /(values.PhiOX() * values.DFe() * values.RhoOX() * (1 - values.PhiOX()))
#      FF = (0.476 * (1.101 + values.PhiOX())) / (kd * values.f() + values.km())
#      DM = (AA - BB) / (CC - DD + EE + FF)
#      print AA,BB,CC,DD,EE,FF
#     # DDelta = ((0.476 * DM * (1 - values.PhiOX())) - (Kd * F * (values.S() * CSat_OB - CB))) / 0.723 # variation in oxide thickness
#      return  DM
# for i in range (0,70000):
#     x = rand.normalvariate(mu, sigma)
#     spall.append(x)
#     
# for i in range(0,70000):
#     
#     DDeltaOX = (DDeltaFunction(Delta, CSat_MP, CSat_PO, CSat_OB, CB, Ce_OB, kd))
#     DMOX = CorrosionRateFunction(Delta, CSat_MP, CSat_PO, CSat_OB, CB)
#          
#     # Runge Kutta method for calculating the Oxide thickness
# 
#     K1 = values.t() * float(DDeltaOX)
#     K2 = values.t() * DDeltaFunction(Delta + (K1 / 2), CSat_MP, CSat_PO, CSat_OB, CB,Ce_OB, kd)
#     K3 = values.t() * DDeltaFunction(Delta + (K2 / 2), CSat_MP, CSat_PO, CSat_OB, CB,Ce_OB, kd)
#     K4 = values.t() * DDeltaFunction(Delta + K3, CSat_MP, CSat_PO, CSat_OB, CB,Ce_OB, kd)
#     DeltaOX = (Delta + (K1 + (2 * K2) + (2 * K3) + K4) / 6)
#     
#                 
#     # Build up of Cr in oxide layer
#     #FeCrBuildUp =  ( 1 / 0.464) * DMOX * (CrContent / 100) * values.t()
#     #FeCr2O4 = FeCr2O4 + FeCrBuildUp 
#     #Fe3O4 = DeltaOX - FeCr2O4
#     #gFe_Fe3O4 = 0.77 * Fe3O4
#     #gFe_FeCr2O4 = 0.2495 * FeCr2O4
#     #gCr_FeCr2O4 = 0.464 * FeCr2O4
#     #PCr_Fe = (gCr_FeCr2O4 / (gFe_Fe3O4 + gFe_FeCr2O4)) 
#     #CSat_PO =  PCr_Fe  * CSat_Cr + (1 - PCr_Fe ) * CSat_PO
#     #print  PCr_Fe, CSat_PO
#     #Erosion effect
#     ParticleDiameter = float(spall[j])
#     SpallingTime = ksp * math.pow(ParticleDiameter * 1e-7,2) / (FF * math.pow(U,2) * values.PhiOX() * kd * values.f() * (values.S() * CSat_OB - Ce_OB))
#     
# #     SpallingTime = ksp * math.pow(ParticleDiameter * 1e-7,2) / (FF * math.pow(U,2) * values.PhiOX() * .2 * values.f() * (values.S() * 1.16e-7 ))
#     
#     if (DummyTime > SpallingTime):# and ParticleDiameter * 1e-7 * values.RhoOX() < DeltaOX / 2):
#         DeltaOX = DeltaOX - ParticleDiameter * 1e-7 * values.RhoOX() 
#         DummyTime = t
#         j = j + 1
#         #FeCr2O4removed = (FeCr2O4 / DeltaOX) * (ParticleDiameter * 1e-7 * values.RhoOX())
#         #FeCr2O4 = FeCr2O4 - FeCr2O4removed
#         #Fe3O4 = DeltaOX - FeCr2O4
#         #gFe_Fe3O4 = 0.77 * Fe3O4
#         #gFe_FeCr2O4 = 0.2495 * FeCr2O4
#         #gCr_FeCr2O4 = 0.464 * FeCr2O4
#         #PCr_Fe = (gCr_FeCr2O4 / (gFe_Fe3O4 + gFe_FeCr2O4)) 
#         #CSat_PO =  PCr_Fe * CSat_Cr + (1 - PCr_Fe ) * CSat_PO
#         #CrConc = 0.464 * FeCr2O4 / (DeltaOX / 5.2)
# #         y = 0.0
# #         while y == 0.0:
# #             x = rand.normalvariate(mu, sigma)
# #             if (x > 10 and x * 1e-7 * values.RhoOX() <  (DeltaOX / 2)):
# #                 spall = x 
# #                 y = 1
#     else:
#         DummyTime = DummyTime + t    
# 
#     
#     if DeltaOX <= 0.0:
#          DeltaOX = 0.0
# #     X = ConcentrationCal.Concentrations(Ce_MP, CeH2_MP,Ce_OB, CeH2_OB, DeltaOX, Temperature, ConcC4H9ONTotal, ActivationEFe_Fe2_MP, ActivationEH2_Hp_MP, ActivationEH2_Hp_OB, ActivationEFe2p_Fe3O4_OB)
# #     Y = X.ConcentrationsCalculator()
#     Thickness.append((DeltaOX *10000.0) / values.RhoOX())
#     corrosionnrate.append((DMOX *10.0 * 3600.0 * 24.0 *365.0)/ values.RhoOX()) 
# #     CSat_MP = Y[0] 
# #     CSat_OB = Y[1]
# #     CeH2_MP = Y[2]
# #     CeH2_OB = Y[3]
# #     Ce_OB = Y[4]
# #     Ce_MP = Y[5]
# #     kd = Y[6]
#     Delta = DeltaOX
#     Time.append(t*i/(3600*24))
# 
# print   DDeltaOX, DMOX*10 * 3600 * 24 *365/ values.RhoOX(), DeltaOX, Thickness, CSat_PO, CSat_MP,Ce_MP, CSat_OB, Ce_OB, CeH2_MP,CeH2_OB,kd
# # 
# plt.figure()
# font ={'family':'Computer Modern Roman', 'size':16}
# plt.rc('font', **font)
# plt.rc('text', usetex= True)
# plt.plot(Time,Thickness, '-')
# plt.title('Thickness vs time' )
# plt.xlabel('Day')
# plt.ylabel('thickness (um)')
# plt.show()
#  
# plt.figure()
# font ={'family':'Computer Modern Roman', 'size':16}
# plt.rc('font', **font)
# plt.rc('text', usetex= True)
# plt.plot(Time,corrosionnrate, '-')
# plt.title('corrosion rate vs time' )
# plt.xlabel('Day')
# plt.ylabel('corrosion rate (mm/year)')
# plt.show()

import os
import csv
import random as rand
import matplotlib.pyplot as plt
import numpy as np
import math
import ConcentrationCal
import Constants
import ActivationEnergyCal
 
# Assigning the constants and initial values
values = Constants.Const()
Temperature = values.Temperature()
TKelvin = values.Temperature() + 273.15
ConcC4H9ONTotal = values.ConcC4H9ONTotal()
CSat_MP = values.CSat_MP()
CSat_OB = values.CSat_OB()
CSat_PO = values.CSat_PO()
CB = values.CB()
kd = values.kd()
t = values.t()
ksp = values.ksp()
FF = values.FF()
CrContent = values.CrContent()
U = values.U()
Ce_OB = values.Ce_OB() 
Ce_MP = 1.0e-8
CeH2_MP = 8.0e-7 # gr/cm^3; if given in cm^3/kg: 10.0 cm^3/kgH2O *P/(RT)*Molar mass H2*ro water => 10*101325*2.016*0.926/(8.1314*298.15*1e6*1e3)
CeH2_OB = 8.0e-7
CSat_Cr = 1.0e-10
Delta = values.DOXs()
Thickness = list()
corrosionnrate = list()
FeCr2O4 = 0
DummyTime = t
j = 0 
spall = list()
mu, sigma = 0.0, 50.0 # mean and standard deviation
# ActivationEnergyValues = ActivationEnergyCal.ActivationEnergy(Ce_MP,CeH2_MP, Ce_OB, CeH2_OB, Temperature, ConcC4H9ONTotal)
# ActivationEFe_Fe2_MP = ActivationEnergyValues.getActivationEnergyFe_Fe2_MP()
# ActivationEH2_Hp_MP = ActivationEnergyValues.getActivationEnergyH2_Hp_MP()
# ActivationEH2_Hp_OB = ActivationEnergyValues.getActivationEnergyH2_Hp_OB()
# ActivationEFe2p_Fe3O4_OB = ActivationEnergyValues.getActivationEnergyFe2p_Fe3O4_OB()
ActivationEFe_Fe2_MP = 464385.58
ActivationEH2_Hp_MP = 100502.6341
ActivationEH2_Hp_OB = 377712.0#377712.0#132314.8545
ActivationEFe2p_Fe3O4_OB = 291642.3046
 
# print  ActivationEFe_Fe2_MP,ActivationEH2_Hp_MP, ActivationEH2_Hp 
Time = list()


#corrosion rate equation 
def DDeltaFunction(Delta, CSat_MP, CSat_PO, CSat_OB, CB, Ce_OB, kd):
    if Delta <= 0:
        Delta = 0.0
    A = values.S() * (CSat_MP + CSat_PO)
    B = ((kd * values.f() * values.S() * CSat_OB) + values.km() * CB) / ((kd * values.f()) + values.km())
    C = (values.DeltaP() * values.ChiP())/(values.PhiP() * values.DP()* values.RhoP() * (1 - values.PhiP()))
    D = (values.MFe() * values.DeltaP() * values.ChiP()) / (values.MP() * values.RhoFe() * values.PhiP() * values.DP()) # Why DFe, should not be DP????????????????????
    E = (0.476 * (1.101 + values.PhiOX()) * values.ChiOX() * Delta) /(values.PhiOX() * values.DFe() * values.RhoOX() * (1 - values.PhiOX()))
    F = (0.476 * (1.101 + values.PhiOX())) / ((kd * values.f()) + values.km())
    DM = (A - B) / (C - D + E + F)
    DDelta = (0.476 * DM * (1.0 - values.PhiOX()) - (kd * values.f() * (values.S() * CSat_OB - Ce_OB))) / 0.723 # variation in oxide thickness
    if Delta == 0:
        DDelta = 0.476 * DM * (1.0 - values.PhiOX())
    if DDelta < 0: 
        DDelta = 0.476 * DM * (1.0 - values.PhiOX())   
    return  DDelta
  
  
def CorrosionRateFunction(Delta, CSat_MP, CSat_PO, CSat_OB, CB):
    if Delta <= 0:
        Delta = 0.0 
    AA = values.S() * (CSat_MP + CSat_PO)
    BB = ((kd * values.f() * values.S() * CSat_OB) + values.km() * CB) / ((kd * values.f()) + values.km())
    CC = (values.DeltaP() * values.ChiP())/(values.PhiP() * values.DP() * values.RhoP() * (1 - values.PhiP()))
    DD = (values.MFe() * values.DeltaP() * values.ChiP()) / (values.MP() * values.RhoFe() * values.PhiP() * values.DP()) # Why DFe, should not be DP????????????????????
    EE = (0.476 * (1.101 + values.PhiOX()) * values.ChiOX() * Delta) /(values.PhiOX() * values.DFe() * values.RhoOX() * (1 - values.PhiOX()))
    FF = (0.476 * (1.101 + values.PhiOX())) / ((kd * values.f()) + values.km())
    DM = (AA - BB) / (CC - DD + EE + FF)
    print AA,BB,CC,DD,EE,FF
    # DDelta = ((0.476 * DM * (1 - values.PhiOX())) - (Kd * F * (values.S() * CSat_OB - CB))) / 0.723 # variation in oxide thickness
    return  DM
 
i=0
while i < 20000:
    x = rand.normalvariate(mu, sigma)
    if x > 10.0:
        spall.append(x)
        i = i+1
     
for i in range(0,20000):
     
    DDeltaOX = (DDeltaFunction(Delta, CSat_MP, CSat_PO, CSat_OB, CB, Ce_OB, kd))
    DMOX = CorrosionRateFunction(Delta, CSat_MP, CSat_PO, CSat_OB, CB)
          
    # Runge Kutta method for calculating the Oxide thickness
    # Diff = 0.476 * DMOX * (1.101 + values.PhiOX())
    # Diss = kd * values.f() * (values.S() * Ce_OB - CB)
 
    K1 = values.t() * float(DDeltaOX)
    K2 = values.t() * DDeltaFunction(Delta + (K1 / 2.0), CSat_MP, CSat_PO, CSat_OB, CB,Ce_OB, kd)
    K3 = values.t() * DDeltaFunction(Delta + (K2 / 2.0), CSat_MP, CSat_PO, CSat_OB, CB,Ce_OB, kd)
    K4 = values.t() * DDeltaFunction(Delta + K3, CSat_MP, CSat_PO, CSat_OB, CB,Ce_OB, kd)
    DeltaOX = (Delta + (K1 + (2.0 * K2) + (2.0 * K3) + K4) / 6.0)
#     print K1, K2, K3, K4, DeltaOX
    if DeltaOX <= 0.0:
        DeltaOX = 1.0e-5 
                 
    # Build up of Cr in oxide layer
#     FeCrBuildUp =  ( 1 / 0.464) * DMOX * (CrContent / 100) * values.t()
#     FeCr2O4 = FeCr2O4 + FeCrBuildUp 
#     Fe3O4 = DeltaOX - FeCr2O4
#     gFe_Fe3O4 = 0.77 * Fe3O4
#     gFe_FeCr2O4 = 0.2495 * FeCr2O4
#     gCr_FeCr2O4 = 0.464 * FeCr2O4
#     PCr_Fe = (gCr_FeCr2O4 / (gFe_Fe3O4 + gFe_FeCr2O4)) 
#     CSat_PO =  PCr_Fe  * CSat_Cr + (1 - PCr_Fe ) * CSat_PO
    #print  PCr_Fe, CSat_PO
    #Erosion effect
    ParticleDiameter = float(spall[j])
    SpallingTime = 1*ksp * math.pow(ParticleDiameter * 1e-7,1) / (FF * math.pow(U,2) * values.PhiOX() * kd * values.f() * (values.S() * CSat_OB - Ce_OB))
     
#     SpallingTime = ksp * math.pow(ParticleDiameter * 1e-7,2) / (FF * math.pow(U,2) * values.PhiOX() * .2 * values.f() * (values.S() * 1.16e-7 ))
     
    if (DummyTime > SpallingTime and ParticleDiameter * 1e-7 * values.RhoOX() < (DeltaOX / 4) ):
        DeltaOX = DeltaOX - ParticleDiameter * 1e-7 * values.RhoOX() 
        DummyTime = values.t()
        j = j + 1
#         FeCr2O4removed = (FeCr2O4 / DeltaOX) * (ParticleDiameter * 1e-7 * values.RhoOX())
#         FeCr2O4 = FeCr2O4 - FeCr2O4removed
#         Fe3O4 = DeltaOX - FeCr2O4
#         gFe_Fe3O4 = 0.77 * Fe3O4
#         gFe_FeCr2O4 = 0.2495 * FeCr2O4
#         gCr_FeCr2O4 = 0.464 * FeCr2O4
#         PCr_Fe = (gCr_FeCr2O4 / (gFe_Fe3O4 + gFe_FeCr2O4)) 
#         CSat_PO =  PCr_Fe * CSat_Cr + (1 - PCr_Fe ) * CSat_PO
#         CrConc = 0.464 * FeCr2O4 / (DeltaOX / 5.2)

    else:
        DummyTime = DummyTime + t    
 
     

    X = ConcentrationCal.Concentrations(Ce_MP, CeH2_MP,Ce_OB, CeH2_OB, DeltaOX, Temperature, ConcC4H9ONTotal, ActivationEFe_Fe2_MP, ActivationEH2_Hp_MP, ActivationEH2_Hp_OB, ActivationEFe2p_Fe3O4_OB)
    Y = X.ConcentrationsCalculator()
    Thickness.append(DeltaOX *10000 / values.RhoOX())
    corrosionnrate.append(DMOX *10 * 3600 * 24 *365/ values.RhoOX()) 
    CSat_MP = float(Y[0]) 
    CSat_OB = float(Y[1])
    CeH2_MP = float(Y[2])
    CeH2_OB = float(Y[3])
    Ce_OB = float(Y[4])
    Ce_MP = float(Y[5])
    kd = float(Y[6])
    Emix =  float(Y[7])
    Ecor = float(Y[8])
    Delta = DeltaOX
    Time.append(t*i/(3600*24))
   
#     print CSat_PO, CSat_MP,Ce_MP, CSat_OB, Ce_OB, CeH2_MP,CeH2_OB,kd, Emix, Ecor
print   DDeltaOX, DMOX*10 * 3600 * 24 *365/ values.RhoOX(), DeltaOX, Thickness, CSat_PO, CSat_MP,Ce_MP, CSat_OB, Ce_OB, CeH2_MP,CeH2_OB,kd
# 
plt.figure()
font ={'family':'Computer Modern Roman', 'size':16}
plt.rc('font', **font)
plt.rc('text', usetex= True)
plt.plot(Time,Thickness, '-')
plt.title('Thickness vs time' )
plt.xlabel('Day')
plt.ylabel('thickness (um)')
plt.show()
  
plt.figure()
font ={'family':'Computer Modern Roman', 'size':16}
plt.rc('font', **font)
plt.rc('text', usetex= True)
plt.plot(Time,corrosionnrate, '-')
plt.title('corrosion rate vs time' )
plt.xlabel('Day')
plt.ylabel('corrosion rate (mm/year)')
plt.show()

# mu, sigma = 0.0, 200.0
# spall =list()
# i=0
# while i < 20000:
#     x = rand.normalvariate(mu, sigma)
#     if x > 10.0:
#         spall.append(x)
#         i = i+1
# count, bins, ignored = plt.hist(spall, 20)
# # plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2) ), linewidth=2, color='r')
# plt.show()
        