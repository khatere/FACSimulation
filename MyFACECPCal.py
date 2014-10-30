import os
import csv
import random as rand
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import fsolve
import Constants
import GibbsEnergyClass
import PHCalculator
class Concentrations(object):
    
    def __init__(self, CSat_OB, Ce_MP, Ce_OP, Ce_OB,  DeltaOX, DMOX, Temperature, ConcC4H9ONTotal, ActivationEFe_Fe2_MP, ActivationEH2_Hp_MP, ActivationEH2_Hp_OB, ActivationEFe2p_Fe3O4_OB, ActivationEH2_Hp_OP, ActivationEFe2p_Fe3O4_OP,pH):
         
        self.Ce_MP = Ce_MP
        self.Ce_OB = Ce_OB
        self.Ce_OP = Ce_OP
        self.CSat_OB = CSat_OB
        self.DeltaOX = DeltaOX
        self.DMOX = DMOX
        self.Temperature = Temperature
        self.ConcC4H9ONTotal = ConcC4H9ONTotal
        self.ActivationEFe_Fe2_MP = ActivationEFe_Fe2_MP
        self.ActivationEH2_Hp_MP = ActivationEH2_Hp_MP
        self.ActivationEH2_Hp_OB = ActivationEH2_Hp_OB
        self.ActivationEFe2p_Fe3O4_OB = ActivationEFe2p_Fe3O4_OB
        self.ActivationEH2_Hp_OP = ActivationEH2_Hp_OP
        self.ActivationEFe2p_Fe3O4_OP = ActivationEFe2p_Fe3O4_OP
        self.pH = pH 

    def ConcentrationsCalculator(self,  DeltaOX = None, DMOX= None, Ce_MP = None, Ce_OP = None, Ce_OB = None, pH = None):
        gibbsValues = GibbsEnergyClass.gibbs(self.Temperature)
        DeltaGFe_Fe2_MP = gibbsValues.getDeltaGFe_Fe2()
        DeltaGFe2p_Fe3O4_OP = 60.0#gibbsValues.getDeltaGFe2p_Fe304()
        DeltaGFe2p_Fe3O4_OB = 60.0#gibbsValues.getDeltaGFe2p_Fe304()
        DeltaGH2_Hp_OP = gibbsValues.getDeltaGfH2_Hp()
        DeltaGH2_Hp_OB = gibbsValues.getDeltaGfH2_Hp()
        DeltaGH2_Hp_MP = gibbsValues.getDeltaGfH2_Hp()
        kHenry = gibbsValues.getH2HenryConstant()
       
#         pH = 6.88 # for now not using the PHCalculator function
       
        
        #pH = 9.58
        value = Constants.Const()
        n = value.n()
        F = value.F()                 
        R = value.R()                  
        h = value.h()           
        kboltzman = value.kboltzman()   
        Beta = value.Beta()           
        PhiOX =  value.PhiOX()
        ChiOX =  value.ChiOX()
        RhoOX =  value.RhoOX()
        DFe =  value.DFe()
        PhiP =  value.PhiP()
        ChiP =  value.ChiP()
        RhoP =  value.RhoP()
        DeltaP = value.DeltaP()
        DP =  value.DP()
        kd_OB =  value.kd()
        kd_OP =  value.kd()
        hH2 =  0.3#value.km()  # for now equals to km 
        Temperature = value.Temperature()
        ConcC4H9ONTotal = value.ConcC4H9ONTotal() 
        FeMolarMass = value.FeMolarMass() 
        H2MolarMass = value.H2MolarMass() 
        RhoH2O = value.RhoH2O()
        CH2coolant = value.CH2coolant()          
        TKelvin = self.Temperature + 273.15;
        
        if Ce_MP == None:
                Ce_MP = self.Ce_MP
        if Ce_OP == None:
                Ce_OP = self.Ce_OP      
#         if CeH2_MP == None:
#                 CeH2_MP = self.CeH2_MP
        if Ce_OB == None:
                Ce_OB = self.Ce_OB   
#         if CeH2_OB == None:
#                 CeH2_OB = self.CeH2_OB
        if DeltaOX == None:
                DeltaOX = self.DeltaOX           
        if DMOX == None:
                DMOX = self.DMOX 
#         if CSat_MP == None:
#                 CSat_MP = self.CSat_MP
#         if CSat_OB == None:
#                 CSat_OB = self.CSat_OB  
        if pH == None:
                pH = self.pH                
        CeHp = math.pow (10, -pH)        
        #Diffusivity of H2
        DH2 = 0.0000222 * TKelvin * math.exp(-12400 / (R * TKelvin))
        
        
        Corrosionrate = DMOX                      
        AP = (RhoP * (1 - PhiP) * DH2 * PhiP)/(DeltaP * ChiP)
        AOX= (RhoOX * (1 - PhiOX) * DH2 * PhiOX)/(DeltaOX * ChiOX)
        
        # H2 concentration at M_P interface assuming high mass transfer to the bulk
        CH2coolant = CH2coolant * 1e-6 * H2MolarMass * 101325 * RhoH2O / (8.314 * 298.15 * 1000)  # cm^3/kg water *1e-6 * H2MolarMass *P/(RT) * rowater/1000 = g H2/cm^3 H2O
        
        CeH2_MP = CH2coolant + (0.005659 * Corrosionrate * (7.32 - PhiOX) *(1 / AOX + 1 / hH2)) + 3.58e-2 *Corrosionrate / AP
        CeH2_OP = CH2coolant + (0.005659 * Corrosionrate * (7.32 - PhiOX) *(1 / AOX + 1 / hH2))
        # H2 concentration at O_B interface
        
        CeH2_OB = CH2coolant + (0.005659 *  Corrosionrate * (7.32 - PhiOX) / hH2)
#         print DeltaOX,CeH2_MP,CeH2_OP,CeH2_OB
        #Dissolution of magnetite at O_B interface Fe3O4 + 2H2O + 2H+ + 2e- = 3 Fe(OH)2 
        
        EeFe2p_Fe3O4_OB = -1 * DeltaGFe2p_Fe3O4_OB * 1000/ (n * F) - math.log(10) * 2 * R * TKelvin / (n * F) * pH - 3 * R * TKelvin / (n * F) * math.log(Ce_OB * 1000 / 90.0) # Equilibrium potential for dissolution reaction at O_B interface
        i0Fe2p_Fe3O4_OB = F * kboltzman * TKelvin / h * (math.exp(-1 * self.ActivationEFe2p_Fe3O4_OB  / (R * TKelvin )) * math.pow(Ce_OB * 1000 / 90.0, 1.0)*  math.exp(-1 * Beta * n * F * EeFe2p_Fe3O4_OB / (R * TKelvin)));
#         print DeltaGFe2p_Fe3O4_OB, Ce_OB
        
#Hydrogen consumption at O_B interface 2H+ + 2e- =  H2 
#Standard DeltaGH2_Hp is zero at any temperature
        
        EeH2_Hp_OB = -1 * DeltaGH2_Hp_OB * 1000 / (n * F ) - math.log(10) * 2 * R * TKelvin * pH / (n * F)  - R * TKelvin / (n * F) * math.log((CeH2_OB )) # At O_B interface
        
        i0H2_Hp_OB = F * kboltzman * TKelvin / h * math.exp(-1 * self.ActivationEH2_Hp_OB / (R * TKelvin)) * math.pow(CeH2_OB * 1000 / H2MolarMass  , 1.0) * math.exp(-1 * Beta * n * F * EeH2_Hp_OB/(R * TKelvin));
        
# Mixed Potential at the O_B interface
        i0_H1= i0H2_Hp_OB * math.exp(Beta * n * F * EeH2_Hp_OB / (R * TKelvin))
        i0_H2 = i0H2_Hp_OB * math.exp(-1 * (1 - Beta) * n * F * EeH2_Hp_OB / (R * TKelvin))
        i0_Fe1 = i0Fe2p_Fe3O4_OB * math.exp(Beta * n * F * EeFe2p_Fe3O4_OB / (R * TKelvin))
        i0_Fe2 = i0Fe2p_Fe3O4_OB * math.exp(-1 * (1 - Beta) * n * F * EeFe2p_Fe3O4_OB / (R * TKelvin))                                                         
        #just for Beta = 0.5
        Emixed = R * TKelvin  / ( n * F * Beta * 2) * math.log((i0_H1 + i0_Fe1) / (i0_H2 + i0_Fe2))
        print EeFe2p_Fe3O4_OB,i0Fe2p_Fe3O4_OB,EeH2_Hp_OB,i0H2_Hp_OB, Emixed,Corrosionrate
#Adjusting the concentration of iron species at O_B interface (in mol/l should be calculated back to g/cm^3)
        CSat_OBNew = (math.exp((Emixed + DeltaGFe2p_Fe3O4_OB * 1000/ (n * F) + 2 * R * TKelvin * math.log(10) * pH / (n * F)) * (-1) * n * F / (3 * R * TKelvin)))* 90 / 1000
#         print CSat_OBNew
#Adjusting the dissolution rate constant of magnetite
        ked_OB = kd_OB * math.exp(-1 * (1 - Beta) * n * F * (Emixed - EeFe2p_Fe3O4_OB) / (R * TKelvin))  
        

#         # Fe concentration at M_P interface
# 
#         CMP_Diffused = ((0.476 * (1.101 + PhiOX) * ChiOX * DeltaOX) * Corrosionrate / (PhiOX * DFe * RhoOX * (1 - PhiOX))) * FeMolarMass / 1000
#         
#         Ce_MP = CSat_MP - CMP_Diffused 

        # At O_P interface
        #Precipitation of magnetite at O_P interface   Fe3O4 + 2H2O + 2H+ + 2e- = 3Fe(OH)2 
        # R * 6 is due to the fact that the n=1/3 in these equation instead of 2
        EeFe2p_Fe3O4_OP = -1 * DeltaGFe2p_Fe3O4_OP * 1000 / (n * F ) - math.log(10) * 2 * R * TKelvin / (n * F) * pH - 3 * R * TKelvin / (n * F) * math.log(Ce_OP * 1000 / FeMolarMass)
        i0Fe2p_Fe3O4_OP = F * kboltzman * TKelvin / h * (math.exp(-1 * self.ActivationEFe2p_Fe3O4_OP  / (R * TKelvin )) * math.pow(Ce_OP * 1000 / 90.0, 1.0)*  math.exp(-1 * Beta * n * F * EeFe2p_Fe3O4_OP / (6 * R * TKelvin)));
        
        #Hydrogen production 2H+ + 2e- = H2 at O_P; 
               
        EeH2_Hp_OP = -1 * DeltaGH2_Hp_OP * 1000 / (n * F ) - math.log(10) * 2 * R * TKelvin * pH / (n * F)  - R * TKelvin / (n * F) * math.log((CeH2_OP )) # At M_P interface, CH2_MP converted to mol/l, g/cm^3 * 1000/Molar mass H2 and then to atm: mol/l/KHenry(mol/l.atm) => atm
        
        i0H2_Hp_OP = F * kboltzman * TKelvin / h * (math.exp(-1 * self.ActivationEH2_Hp_OP / (R * TKelvin)) *  math.pow(CeH2_OP * 1000  / H2MolarMass , 1.0) * math.exp(-1 * Beta * n * F * EeH2_Hp_OP / (6 * R * TKelvin)))
        
       # Mixed Potential at the O_P interface
        i0_H1OP= i0H2_Hp_OP * math.exp(Beta * n * F * EeH2_Hp_OP / (6 * R * TKelvin))
        i0_H2OP = i0H2_Hp_OP * math.exp(-1 * (1 - Beta) * n * F * EeH2_Hp_OP / (6 * R * TKelvin))
        i0_Fe1OP = i0Fe2p_Fe3O4_OP * math.exp(Beta * n * F * EeFe2p_Fe3O4_OP / (6 * R * TKelvin))
        i0_Fe2OP = i0Fe2p_Fe3O4_OP * math.exp(-1 * (1 - Beta) * n * F * EeFe2p_Fe3O4_OP / (6 * R * TKelvin))                                                         
        #just for Beta = 0.5
        EmixedOP = 6 * R * TKelvin  / ( n * F * Beta * 2) * math.log((i0_H1OP + i0_Fe1OP) / (i0_H2OP + i0_Fe2OP))
        #Adjusting the dissolution rate constant of magnetite
        ked_OP = kd_OP * math.exp(-1 * (1 - Beta) * n * F * (EmixedOP - EeFe2p_Fe3O4_OP) / (6 * R * TKelvin))
        #Adjusting the concentration of iron species at O_P interface (in mol/l should be calculated back to g/cm^3)
        CSat_OPNew = (math.exp((EmixedOP + DeltaGFe2p_Fe3O4_OP * 1000/ (n * F) + 2 * R * TKelvin * math.log(10) * pH / (n * F)) * (-1) * n * F / (3 * R * TKelvin)))* 90 / 1000
        Ce_OPNew = CSat_OPNew + (0.476* (1-PhiOX)*Corrosionrate / ked_OP) * 90 / 1000
        Ce_MPNew = (Corrosionrate*(value.ChiP() * value.DeltaP()) /(value.PhiP() * value.DP() * value.RhoP() * (1 - value.PhiP())) + Ce_OPNew *1000/90)*56/1000
        print EeFe2p_Fe3O4_OP,i0Fe2p_Fe3O4_OP,EeH2_Hp_OP,i0H2_Hp_OP, EmixedOP,Ce_MPNew,Ce_OP
#         # At M-P interface
# 
#         EeFe_Fe2p_MP = -1 * DeltaGFe_Fe2_MP * 1000 / (n * F) - math.log(10) * 2 * R * TKelvin * pH / (n * F)  + R * TKelvin / (n * F) * math.log(Ce_MP * 1000 / FeMolarMass) # Equilibrium potential for the oxidation reaction at M_P interface, concentration unit changed from g/cm^3 to mol/L 
#         #in i0 equationsthe unit of concentration should be mol/cm^3 according to Dr. Cook thesis and Olga program ( the power 2/3 is not clear and the units does not match the A/cm^2 )  but according to Lisa Lang program Concentration in mol/lit should be multiplied by diameter of pipe/4 (cm). This would make sense in terms of units. I tried all the variations and they make insignificant changes in the results
#         i0Fe_Fe2p_MP = (F * kboltzman * TKelvin / h) * (math.exp(-1 * self.ActivationEFe_Fe2_MP  / (R * TKelvin)) * math.pow(Ce_MP * 1000  / FeMolarMass,1.0) * math.exp(-1 * Beta * n * F * EeFe_Fe2p_MP / (R * TKelvin)));
        
#Precipitation of magnetite at M_P interface   Fe3O4 + 2H2O + 2H+ + 2e- = 3Fe(OH)2 
        
#       EeFe2p_Fe3O4_MP = -1 * DeltaGFe2p_Fe3O4_MP * 1000 / (n * F ) - math.log(10) * 2 * R * TKelvin / (n * F) * pH - 3 * R * TKelvin / (n * F) * math.log(Ce_MP * 1000 / FeMolarMass);
        
        EeFe_Fe2p_MP = -1 * DeltaGFe_Fe2_MP* 1000 / (n * F) - math.log(10) * 2 * R * TKelvin / (n * F) * pH + R * TKelvin / (n * F) * math.log(Ce_MPNew * 1000 / FeMolarMass) # Equilibrium potential for the oxidation reaction at M_P interface
        #in i0 equationsthe unit of concentration should be mol/cm^3 according to Dr. Cookthesis and Olga program ( the power 2/3 is not clear and the units does not match the A/cm^2 )  but according to Lisa Lang program Concentration in mol/lit should be multiplied by diameter of pipe/4 (cm). This would make sense in terms of units. I tried all the variations and they make insignificant changes in the results
        i0Fe_Fe2p_MP = F * kboltzman * TKelvin / h * (math.exp(-1 * self.ActivationEFe_Fe2_MP  / (R * TKelvin)) * math.pow(Ce_MPNew * 1000 / FeMolarMass,1.0) * math.exp(-1 * Beta * n * F * EeFe_Fe2p_MP / (R * TKelvin)));
        
        #Hydrogen production 2H+ + 2e- = H2 at M_P; 
               
        EeH2_Hp_MP = -1 * DeltaGH2_Hp_MP * 1000 / (n * F ) - math.log(10) * 2 * R * TKelvin * pH / (n * F)  - R * TKelvin / (n * F) * math.log((CeH2_MP )) # At M_P interface, CH2_MP converted to mol/l, g/cm^3 * 1000/Molar mass H2 and then to atm: mol/l/KHenry(mol/l.atm) => atm
        
        i0H2_Hp_MP = F * kboltzman * TKelvin / h * (math.exp(-1 * self.ActivationEH2_Hp_MP / (R * TKelvin)) *  math.pow(CeH2_MP * 1000  / H2MolarMass , 1.0) * math.exp(-1 * Beta * n * F * EeH2_Hp_MP / (R * TKelvin)))
        
        

# having the corrosion rate, the Ecorr can be calculated by solving the following equation:
# Corrosionrate = iFe_Fe2p_MP * FeMolarMass/ (n * F)
# iFe_Fe2p_MP = i0Fe_Fe2p_MP * (math.exp(Beta * n * F / (R * TKelvin) * (Ecorr - EeFe_Fe2p_MP)) - math.exp(-1 * (1 - Beta) * n * F / (R * TKelvin) * (Ecorr - EeFe_Fe2p_MP)))
# i0Fe_Fe2p_MP * (math.exp(Beta * n * F / (R * TKelvin) * (Ecorr - EeFe_Fe2p_MP)) - math.exp(-1 * (1 - Beta) * n * F / (R * TKelvin) * (Ecorr - EeFe_Fe2p_MP)))* FeMolarMass/ (n * F) - Corrosionrate =0.0
        A1 = Corrosionrate * (n * F) / (FeMolarMass * i0Fe_Fe2p_MP)
        C1 = Beta * n * F / (R * TKelvin)
        C2 = math.exp(C1 * (EeFe_Fe2p_MP))
        X1 = (A1 + math.pow(math.pow(A1,2) + 4, 0.5)) / 2
        X2 = (A1 - math.pow(math.pow(A1,2) + 4, 0.5)) / 2
        
#         C1 = Beta * n * F / (R * TKelvin)
#         C2 = (1-Beta) * n * F / (R * TKelvin)
#         C3 = math.exp(C1 * (EeFe_Fe2p_MP))
#         C4 = math.exp(C2 * (EeFe_Fe2p_MP))
#         print C1+C2, A1*C3, C2, C3*C4

        
        if X1 > 0: 
            Ecorr = math.log(X1) / C1 + EeFe_Fe2p_MP
        elif X2 > 0:
            Ecorr = math.log(X2) / C1 + EeFe_Fe2p_MP
        else:
            print ('error')
        if (X1 > 0 and X2 > 0):
            Ecorr = math.log(X2) / C1 + EeFe_Fe2p_MP 
            
          
#         print A1,C1, C2, X1, X2, Ecorr, Emixed,
        print EeFe_Fe2p_MP,i0Fe_Fe2p_MP,EeH2_Hp_MP,i0H2_Hp_MP,Ecorr
        

# Current calculations at M_P interface
# iH2_Hp_MP = i0H2_Hp_MP * (math.exp(Beta * n * F / (R * TKelvin) * (EeH2_Hp_MP - Ecorr)) - math.exp(-1 * (1 - Beta) * n * F / (R * TKelvin) * (EeH2_Hp_MP - Ecorr)))
        

        
#Adjusting Ce_MP for the Ecorr (in mol/l should be calculated back to g/cm^3)
        CSat_MPNew = (math.exp((Ecorr-EeFe_Fe2p_MP + DeltaGFe_Fe2_MP * 1000 / (n * F ) + math.log(10) * 2 * R * TKelvin * pH / (n * F)) * (1) * n * F / (R * TKelvin)))* FeMolarMass / 1000
        
#         CSat_MPNew = Ce_MPNew + CMP_Diffused         
        #Solubility at M-P: CMPsat assuming that the solubility at the M_P interface is summation of what is diffused at this interface and the equilibrium concentration considering the developed Ecorr
       
# Current calculations at O_B interface
# iFe3O4_Fe2p_OB = i0Fe2p_Fe3O4_OB * (math.exp(Beta * n * F / (R * TKelvin) * (EeFe2p_Fe3O4_OB - Emixed)) - math.exp(-1 * (1 - Beta) * n * F / (R * TKelvin) * (EeFe2p_Fe3O4_OB - Emixed)));
# iH2_Hp_OB = i0H2_Hp_OB * (math.exp(Beta * n * F / (R * TKelvin) * (Emixed - EeH2_Hp_OB)) - math.exp(-1 * (1 - Beta) * n * F / (R * TKelvin) * (Emixed - EeH2_Hp_OB)));
        
        

         
    


        
        
        
        return CSat_MPNew, CSat_OPNew, CSat_OBNew, ked_OB, ked_OP, Emixed, EmixedOP, Ecorr, CeH2_MP, CeH2_OP, CeH2_OB, Corrosionrate,Ce_MPNew,Ce_OPNew    