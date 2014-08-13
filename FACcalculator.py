import os
import csv
import random as rand
import matplotlib.pyplot as plt
import numpy as np
import math
import ConcentrationCal

dirName = '\Users\Khatere\Documents\sampython'
fileName = 'InputFileFAC.csv'
 
with open(os.path.join(dirName,fileName),'r') as csvfile:
    dataReader = csv.reader(csvfile, delimiter = ',')
    InputData = list()
    for row in dataReader:
        InputData.append(row)
    InitialValues = [float(InputData[i][1]) for i in range(len(InputData))]
    Dataname = [InputData[i][0] for i in range(len(InputData))]
    print 'InitialValues =', InitialValues
    #print 'Dataname = ', Dataname[0]

CSat_MP = InitialValues[0]
CSat_PO = InitialValues[1]
kd = InitialValues[2]
km = InitialValues[3]
f = InitialValues[4]
S = InitialValues[5]
CB = InitialValues[6]
DeltaP = InitialValues[7]
ChiP = InitialValues[8]
PhiP = InitialValues[9]
DFe = InitialValues[10]
RhoP = InitialValues[11]
MFe = InitialValues[12]
MP = InitialValues[13]
RhoFe = InitialValues[14] 
PhiOX = InitialValues[15]
ChiOX = InitialValues[16]
RhoOX = InitialValues[17]
CPO = InitialValues[18]
DCORs = InitialValues[19]
DOXs = InitialValues[20]
t = InitialValues[21]
Tstep = InitialValues[22]
DP = InitialValues[23]
CSat_OB = InitialValues[24]
Ce_OB = InitialValues[25]
CrContent = InitialValues[26]
Paverage = InitialValues[27]
SDP = InitialValues[28]
TS = InitialValues[29]
TFs = InitialValues[30]
U = InitialValues[31]
Tphi = InitialValues[32]
TKd = InitialValues[33]
TFSA = InitialValues[34]
TFeSAT = InitialValues[35]
Tfev = InitialValues[36]
RhoH20 = InitialValues[37]
Ksp = InitialValues[38]
Cell = InitialValues[39]

#DeltaOX = 2
DeltaOX = [0 for i in range(1,10)]
DDeltaOX = [0 for i in range(1,10)]
DMOX = [0 for i in range(1,10)]
# corrosion rate equation 
def DDeltaFunction(Delta, CSat_MP, CSat_PO, CSat_OB, CB):
    AA = S * (CSat_MP - CSat_PO)
    BB = ((kd * f * S * CSat_OB) + km * CB) / ((kd * f) + km)
    CC = (DeltaP + ChiP)/(PhiP * DP * RhoP * (1 - PhiP))
    DD = (MFe * DeltaP * ChiP * (1 - PhiP)) / (MP * RhoFe * PhiP * DFe * (1 - PhiP)) # Why DFe, should not be DP????????????????????
    EE = (0.476 * (1.101 + PhiOX) * ChiOX * Delta) /(PhiOX * DFe * RhoOX * (1 - PhiOX))
    FF = (0.476 * (1.101 + PhiOX)) / ((kd * f) + km)
    DM = (AA- BB) / (CC - DD + EE + FF)
    DDelta = ((0.476 * DM * (1 - PhiOX)) - (kd * f * (S * CSat_OB - CB))) / 0.723 # variation in oxide thickness
    return DDelta
 
def CorrosionRateFunction(Delta, CSat_MP, CSat_PO, CSat_OB, CB):
    AA = S * (CSat_MP - CSat_PO)
    BB = ((kd * f * S * CSat_OB) + km * CB) / ((kd * f) + km)
    CC = (DeltaP + ChiP)/(PhiP * DP * RhoP * (1 - PhiP))
    DD = (MFe * DeltaP * ChiP * (1 - PhiP)) / (MP * RhoFe * PhiP * DFe * (1 - PhiP)) # Why DFe, should not be DP????????????????????
    EE = (0.476 * (1.101 + PhiOX) * ChiOX * Delta) /(PhiOX * DFe * RhoOX * (1 - PhiOX))
    FF = (0.476 * (1.101 + PhiOX)) / ((kd * f) + km)
    DM = (AA - BB) / (CC - DD + EE + FF)
   # DDelta = ((0.476 * DM * (1 - PhiOX)) - (Kd * F * (S * CSat_OB - CB))) / 0.723 # variation in oxide thickness
    return  DM
# random number generator 
for i in range(0, 8):
    DDeltaOX[i] = DDeltaFunction(DeltaOX[i])
    DMOX[i] = CorrosionRateFunction(DeltaOX[i])
# Runge Kutta method for calculating the Oxide thickness
    Diff = 0.476 * DMOX[i] * (1.101 + PhiOX)
    Diss = kd * f * (S * Ce_OB - CB)
    K1 = t * DDeltaOX[i]
    K2 = t * DDeltaFunction(DeltaOX[i] + (K1 / 2))
    K3 = t * DDeltaFunction(DeltaOX[i] + (K2 / 2))
    K4 = t * DDeltaFunction(DeltaOX[i] + K3)
    DeltaOX[i+1] = DeltaOX[i] + (K1 + (2 * K2) + (2 * K3) + K4) / 6
    X = ConcentrationCal.Concentrations(Ce_MP, CeH2_MP,Ce_OB, CeH2_OB)
    CSat_MP = X[0] 
    CSat_OB = X[1]
    CeH2_MP = X[2]
    CeH2_OB = X[3]
    Ce_OB = X[4]
    kd = X[5]
#DeltaOX = DeltaOXe
print DeltaOX, DDeltaOX, DMOX
print DP


