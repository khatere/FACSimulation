 
import os
import csv
import random as rand
import matplotlib.pyplot as plt
import numpy as np
import math
import MyFACECPCal
import Constants
import ActivationEnergyCal
import PHCalculator
  
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
Ce_MP = values.Ce_MP()
Ce_OP = values.Ce_OP()
# Ce_MP = 1.0e-8
# CeH2_MP = 8.0e-7 # gr/cm^3; if given in cm^3/kg: 10.0 cm^3/kgH2O *P/(RT)*Molar mass H2*ro water => 10*101325*2.016*0.926/(8.1314*298.15*1e6*1e3)
# CeH2_OB = 8.0e-7
CSat_Cr = 1.0e-10
Delta = values.DOXs()
Thickness = list()
corrosionnrate = list()
FeCr2O4 = 0
DummyTime = t
j = 0 
spall = list()
mu, sigma = 0.0, 100.0 # mean and standard deviation
# ActivationEnergyValues = ActivationEnergyCal.ActivationEnergy(Ce_MP,CeH2_MP, Ce_OB, CeH2_OB, Temperature, ConcC4H9ONTotal)
# ActivationEFe_Fe2_MP = ActivationEnergyValues.getActivationEnergyFe_Fe2_MP()
# ActivationEH2_Hp_MP = ActivationEnergyValues.getActivationEnergyH2_Hp_MP()
# ActivationEH2_Hp_OB = ActivationEnergyValues.getActivationEnergyH2_Hp_OB()
# ActivationEFe2p_Fe3O4_OB = ActivationEnergyValues.getActivationEnergyFe2p_Fe3O4_OB()
ActivationEFe_Fe2_MP = 190600.0#264385.58
ActivationEH2_Hp_MP = 150000.0#100502.6341
ActivationEH2_Hp_OB = 210000.0#377712.0#377712.0#132314.8545
ActivationEFe2p_Fe3O4_OB = 190000.0#291642.3046
ActivationEFe2p_Fe3O4_OP = 178000.0
ActivationEH2_Hp_OP = 185000.0
# PHvalue = PHCalculator.equilibriumConstants(Temperature, ConcC4H9ONTotal)
# pH = PHvalue.PHCalculation() 
pH = 5.88  
print  ActivationEFe_Fe2_MP,ActivationEH2_Hp_MP, ActivationEH2_Hp_OB, ActivationEFe2p_Fe3O4_OB  
Time = list()

 
#corrosion rate equation 
def DDeltaFunction(Delta, Ce_MP, CSat_OB, CB, Ce_OB, kd):
    if Delta <= 0:
        Delta = 0.0
    A = values.S() * (Ce_MP)
    B = ((kd * values.f() * values.S() * CSat_OB) + values.km() * CB) / ((kd * values.f()) + values.km())
#     C = 20*(values.DeltaP() * values.ChiP())/(values.PhiP() * values.DP()* values.RhoP() * (1 - values.PhiP()))
#     D = 20*(values.MFe() * values.DeltaP() * values.ChiP()) / (values.MP() * values.RhoFe() * values.PhiP() * values.DP()) # Why DFe, should not be DP????????????????????
    D = 1.0*(values.ChiP() * values.DeltaP()) /(values.PhiP() * values.DP() * values.RhoP() * (1 - values.PhiP()))
    E = (0.476 * (1.101 + values.PhiOX()) * values.ChiOX() * Delta) /(values.PhiOX() * values.DFe() * values.RhoOX() * (1 - values.PhiOX()))
    F = (0.476 * (1.101 + values.PhiOX())) / ((kd * values.f()) + values.km())
    DM = (A - B) / ( D + E + F)
    DDelta = (0.476 * DM * (1.0 - values.PhiOX()) - (kd * values.f() * (values.S() * CSat_OB - Ce_OB))) / 0.723 # variation in oxide thickness
    if Delta == 0:
        DDelta = 0.476 * DM * (1.0 - values.PhiOX())
    if DDelta < 0: 
        DDelta = 0.476 * DM * (1.0 - values.PhiOX())   
    return  DDelta
   
   
def CorrosionRateFunction(Delta, Ce_MP, CSat_OB, CB):
    if Delta <= 0:
        Delta = 0.0 
    AA = values.S() * (Ce_MP)
    BB = ((kd * values.f() * values.S() * CSat_OB) + values.km() * CB) / ((kd * values.f()) + values.km())
#     CC = 20*(values.DeltaP() * values.ChiP())/(values.PhiP() * values.DP() * values.RhoP() * (1 - values.PhiP()))
#     DD = 20*(values.MFe() * values.DeltaP() * values.ChiP()) / (values.MP() * values.RhoFe() * values.PhiP() * values.DP()) # Why DFe, should not be DP????????????????????
    DD = 1.0*(values.ChiP() * values.DeltaP()) /(values.PhiP() * values.DP() * values.RhoP() * (1 - values.PhiP()))
    EE = (0.476 * (1.101 + values.PhiOX()) * values.ChiOX() * Delta) /(values.PhiOX() * values.DFe() * values.RhoOX() * (1 - values.PhiOX()))
    FF = (0.476 * (1.101 + values.PhiOX())) / ((kd * values.f()) + values.km())
    DM = (AA - BB) / ( DD + EE + FF)
#     print AA,BB,DD,EE,FF
    # DDelta = ((0.476 * DM * (1 - values.PhiOX())) - (Kd * F * (values.S() * CSat_OB - CB))) / 0.723 # variation in oxide thickness
    return  DM
  
i=0
while i < 7000:
    x = rand.normalvariate(mu, sigma)
    if x > 10.0 and x < 100.0:
        spall.append(x)
        i = i+1
DeltaOX = Delta      
for i in range(0,7000):
    ParticleDiameter = float(spall[j])        
    SpallingTime = 1*ksp * math.pow(ParticleDiameter * 1e-7,1) / (FF * math.pow(U,2) * values.PhiOX() * kd * values.f() * (values.S() * CSat_OB - Ce_OB))
#     print DeltaOX, SpallingTime, ParticleDiameter * 1e-7 * values.RhoOX(), DeltaOX/2  
#     SpallingTime = ksp * math.pow(ParticleDiameter * 1e-7,2) / (FF * math.pow(U,2) * values.PhiOX() * .2 * values.f() * (values.S() * 1.16e-7 ))
      
    if (DummyTime > SpallingTime and ParticleDiameter * 1e-7 * values.RhoOX() < (DeltaOX/2 )):
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
        CrConc = 0.464 * FeCr2O4 / (DeltaOX / 5.2)
 
    else:
        DummyTime = DummyTime + t   
#     print DeltaOX  
    DDeltaOX = (DDeltaFunction(DeltaOX, Ce_MP,  CSat_OB, CB, Ce_OB, kd))
    DMOX = CorrosionRateFunction(DeltaOX, Ce_MP, CSat_OB, CB)
    if DMOX <0.0:
        DMOX = 1e-8       
    # Runge Kutta method for calculating the Oxide thickness
    # Diff = 0.476 * DMOX * (1.101 + values.PhiOX())
    # Diss = kd * values.f() * (values.S() * Ce_OB - CB)
  
    K1 = values.t() * float(DDeltaOX)
    K2 = values.t() * DDeltaFunction(DeltaOX + (K1 / 2.0), Ce_MP, CSat_OB, CB,Ce_OB, kd)
    K3 = values.t() * DDeltaFunction(DeltaOX + (K2 / 2.0), Ce_MP,  CSat_OB, CB,Ce_OB, kd)
    K4 = values.t() * DDeltaFunction(DeltaOX + K3, Ce_MP, CSat_OB, CB,Ce_OB, kd)
    DeltaOX = (DeltaOX + (K1 + (2.0 * K2) + (2.0 * K3) + K4) / 6.0)
#     print K1, K2, K3, K4, DeltaOX
    if DeltaOX <= 0.0:
        DeltaOX = 1.0e-5 
#     print DeltaOX              
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
     
  
      
    
    
#     print 0.476 * (1.101 + values.PhiOX()) * DMOX, kd * values.f() * values.S() * CSat_OB + values.km() * CB,values.km() + kd * values.f()
#     Ce_MP = ((0.476 * (1.101 + values.PhiOX())) / ((kd * values.f()) + values.km()) + (0.476 * (1.101 + values.PhiOX()) * values.ChiOX() * Delta) /(values.PhiOX() * values.DFe() * values.RhoOX() * (1 - values.PhiOX()))) * DMOX + ((kd * values.f() * values.S() * CSat_OB) + values.km() * CB) / ((kd * values.f()) + values.km())
    
    X = MyFACECPCal.Concentrations( CSat_OB, Ce_MP, Ce_OP, Ce_OB,  DeltaOX, DMOX,  Temperature, ConcC4H9ONTotal, ActivationEFe_Fe2_MP, ActivationEH2_Hp_MP, ActivationEH2_Hp_OB, ActivationEFe2p_Fe3O4_OB,ActivationEH2_Hp_OP,ActivationEFe2p_Fe3O4_OP, pH)
    Y = X.ConcentrationsCalculator()
    CSat_MP = float(Y[0])
    CSat_OP = float(Y[1]) 
    CSat_OB = float(Y[2])  
    kd = float(Y[3])
#     kdOP = float(Y[4])
    Emixed =  float(Y[5])
    EmixedOP = float(Y[5]) 
    Ecorr = float(Y[7])
    CeH2_MP = float(Y[8])
    CeH2_OP = float(Y[9])
    CeH2_OB = float(Y[10])
    corrosion = float(Y[11])
    Ce_MP = float(Y[12])
    Ce_OP = float(Y[13])
    Delta = DeltaOX
    Ce_OB = (0.476 * (1.101 + values.PhiOX()) * DMOX + (kd * values.f() * values.S() * CSat_OB) + values.km() * CB) / (values.km() + kd * values.f() )

#     Ce_MP = CSat_MP + 0.476 * (1 - values.PhiOX())*DMOX / kd / 90
    
    Thickness.append(DeltaOX *10000 / values.RhoOX())
    corrosionnrate.append(DMOX *10 * 3600 * 24 *365/ values.RhoOX()) 
    
    Time.append(t*i/(3600*24))
    
    
#     print CSat_PO, Ce_OP,Ce_MP, CSat_OB, Ce_OB, CeH2_MP,CeH2_OB,kd, Emixed, Ecorr,corrosion,i
print   DDeltaOX, DMOX*10 * 3600 * 24 *365/ values.RhoOX(), DeltaOX, Thickness, CSat_PO, Ce_OP, CSat_MP,Ce_MP, CSat_OB, Ce_OB, CeH2_MP,CeH2_OB,kd,Emixed, Ecorr
# 
plt.figure()
font ={'family':'Computer Modern Roman', 'size':16}
plt.rc('font', **font)
plt.rc('text', usetex= True)
plt.plot(Time,Thickness, '-')
plt.title('Thickness vs time')
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

# mu, sigma = 0.0, 70.0
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
# dirName = '\Users\Khatere\Documents\sampython'
# fileName = 'randomdata.dat'
# 
# if not os.path.exists(dirName):
#     os.makedirs(dirName)
#     
# dataList = list()
# [dataList.append([Thickness[i], Time[i]]) for i in range(len(Thickness))]
# print dataList
# with open(os.path.join(dirName, fileName),'w') as csvfile:
#     writer = csv.writer(csvfile, delimiter = ',')
#     writer.writerow(['Thickness', 'Time'])
#     writer.writerows(dataList)





