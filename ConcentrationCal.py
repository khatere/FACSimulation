import numpy as np
import math
import GibbsEnergyClass
import PHCalculator
import Constants
import ActivationEnergyCal

class Concentrations(object):
    
    def __init__(self, Ce_MP, CeH2_MP,Ce_OB, CeH2_OB, DeltaOX, Temperature, ConcC4H9ONTotal, ActivationEFe_Fe2_MP, ActivationEH2_Hp_MP, ActivationEH2_Hp_OB, ActivationEFe2p_Fe3O4_OB):
         
        self.Ce_MP = Ce_MP
        self.CeH2_MP = CeH2_MP
        self.Ce_OB = Ce_OB
        self.CeH2_OB = CeH2_OB
        self.DeltaOX = DeltaOX
        self.Temperature = Temperature
        self.ConcC4H9ONTotal = ConcC4H9ONTotal
        self.ActivationEFe_Fe2_MP = ActivationEFe_Fe2_MP
        self.ActivationEH2_Hp_MP = ActivationEH2_Hp_MP
        self.ActivationEH2_Hp_OB = ActivationEH2_Hp_OB
        self.ActivationEFe2p_Fe3O4_OB = ActivationEFe2p_Fe3O4_OB
         

    def ConcentrationsCalculator(self, Ce_MP = None, CeH2_MP = None, Ce_OB = None, CeH2_OB = None, DeltaOX = None):
        gibbsValues = GibbsEnergyClass.gibbs(self.Temperature)
        DeltaGFe_Fe2_MP = gibbsValues.getDeltaGFe_Fe2()
        DeltaGFe2p_Fe3O4_MP = gibbsValues.getDeltaGFe2p_Fe304()
        DeltaGFe2p_Fe3O4_OB = gibbsValues.getDeltaGFe2p_Fe304()
        DeltaGH2_Hp_MP = gibbsValues.getDeltaGfH2_Hp()
        DeltaGH2_Hp_OB = gibbsValues.getDeltaGfH2_Hp()
        kHenry = gibbsValues.getH2HenryConstant()
#         PHvalue = PHCalculator.equilibriumConstants(self.Temperature, self.ConcC4H9ONTotal)
#         pH = PHvalue.PHCalculation() 
        pH = 6.88 # for now not using the PHCalculator function
        CeHp = math.pow (10, -pH)
        
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
        kd_OB =  value.kd()
        hH2 =  value.km()  # for now equals to km 
        Temperature = value.Temperature()
        ConcC4H9ONTotal = value.ConcC4H9ONTotal() 
        FeMolarMass = value.FeMolarMass() 
        H2MolarMass = value.H2MolarMass() 
        RhoH2O = value.RhoH2O()
        CH2coolant = value.CH2coolant()          
        TKelvin = self.Temperature + 273.15;
        
        if Ce_MP == None:
                Ce_MP = self.Ce_MP
        if CeH2_MP == None:
                CeH2_MP = self.CeH2_MP
        if Ce_OB == None:
                Ce_OB = self.Ce_OB   
        if CeH2_OB == None:
                CeH2_OB = self.CeH2_OB
        if DeltaOX == None:
                DeltaOX = self.DeltaOX           
                
        
        # At the oxide-Bulk interface O_B
        # Corrosion reaction at the M_P interface; All the reaction are considered in the following direction:  O + ne- -> R
        #  2H+ +2e- + Fe(OH)2 =  Fe + 2H2O 
        # Fe concentration should be converted to mol/L from g/cm^3 in potential equations g/cm^3 * 1000 / Molar mass Fe = mol/L
        EeFe_Fe2p_MP = -1 * DeltaGFe_Fe2_MP* 1000 / (n * F) - math.log(10) * 2 * R * TKelvin / (n * F) * pH + R * TKelvin / (n * F) * math.log(Ce_MP * 1000 / FeMolarMass) # Equilibrium potential for the oxidation reaction at M_P interface
        #in i0 equationsthe unit of concentration should be mol/cm^3 according to Dr. Cookthesis and Olga program ( the power 2/3 is not clear and the units does not match the A/cm^2 )  but according to Lisa Lang program Concentration in mol/lit should be multiplied by diameter of pipe/4 (cm). This would make sense in terms of units. I tried all the variations and they make insignificant changes in the results
        i0Fe_Fe2p_MP = F * kboltzman * TKelvin / h * (math.exp(-1 * self.ActivationEFe_Fe2_MP  / (R * TKelvin)) * math.pow(Ce_MP  / FeMolarMass,0.67) * math.exp(-1 * Beta * n * F * EeFe_Fe2p_MP / (R * TKelvin)));
        
        #Precipitation of magnetite at M_P interface   Fe3O4 + 2H2O + 2H+ + 2e- = 3Fe(OH)2 
        
        EeFe2p_Fe3O4_MP = -1 * DeltaGFe2p_Fe3O4_MP * 1000 / (n * F ) - math.log(10) * 2 * R * TKelvin / (n * F) * pH - 3 * R * TKelvin / (n * F) * math.log(Ce_MP * 1000 / FeMolarMass);
        
        #Hydrogen production 2H+ + 2e- = H2 at M_P; 
        
        
        EeH2_Hp_MP = -1 * DeltaGH2_Hp_MP * 1000 / (n * F ) - math.log(10) * 2 * R * TKelvin / (n * F) * pH - R * TKelvin / (n * F) * math.log((CeH2_MP * 1000  / H2MolarMass) / kHenry) # At M_P interface, CH2_MP converted to mol/l, g/cm^3 * 1000/Molar mass H2 and then to atm: mol/l/KHenry(mol/l.atm) => atm
        
        i0H2_Hp_MP = F * kboltzman * TKelvin / h * (math.exp(-1 * self.ActivationEH2_Hp_MP / (R * TKelvin)) *  math.pow(CeHp/1000, 0.67) * math.exp(-1 * Beta * n * F * EeH2_Hp_MP / (R * TKelvin)))
        
        
        # Mixed Potential at the M_P interface
        Ecorr = R * TKelvin /( n * F) * math.log((i0H2_Hp_MP * math.exp(Beta * n * F * EeH2_Hp_MP / (R * TKelvin)) + i0Fe_Fe2p_MP * math.exp(Beta * n * F * EeFe_Fe2p_MP / (R * TKelvin)))/(i0H2_Hp_MP * math.exp(-1 * (1 - Beta) * n * F * EeH2_Hp_MP / (R * TKelvin)) + i0Fe_Fe2p_MP * math.exp(-1 * (1 - Beta) * n * F * EeFe_Fe2p_MP / (R * TKelvin))));
        
        # Current calculations at M_P interface
        iFe_Fe2p_MP = i0Fe_Fe2p_MP * (math.exp(Beta * n * F / (R * TKelvin) * (Ecorr - EeFe_Fe2p_MP)) - math.exp(-1 * (1 - Beta) * n * F / (R * TKelvin) * (Ecorr - EeFe_Fe2p_MP)));
        iH2_Hp_MP = i0H2_Hp_MP * (math.exp(Beta * n * F / (R * TKelvin) * (EeH2_Hp_MP - Ecorr)) - math.exp(-1 * (1 - Beta) * n * F / (R * TKelvin) * (EeH2_Hp_MP - Ecorr)));
        
        Corrosionrate = iFe_Fe2p_MP * FeMolarMass/ (n * F)
        
        #Adjusting Ce_MP for the Ecorr (in mol/l should be calculated back to g/cm^3)
        Ce_MPNew = (math.exp((Ecorr + DeltaGFe_Fe2_MP * 1000 / (n * F ) + math.log(10) * 2 * R * TKelvin / (n * F) * pH) * (n * F) / (R * TKelvin)))* FeMolarMass / 1000
        
        #Solubility at M-P: CMPsat assuming that the solubility at the M_P interface is summation of what is diffused at this interface and the equilibrium concentration considering the developed Ecorr
        
        CMP_Diffused = ((0.476 * (1.101 + PhiOX) * ChiOX * DeltaOX) * Corrosionrate / (PhiOX * DFe * RhoOX * (1 - PhiOX))) * FeMolarMass / 1000
        
        CSat_MPNew = CMP_Diffused +  Ce_MPNew
        
        
        #Dissolution of magnetite at O_B interface Fe3O4 + 2H2O + 2H+ + 2e- = 3 Fe(OH)2 
        
        EeFe2p_Fe3O4_OB = -1 * DeltaGFe2p_Fe3O4_OB * 1000/ (n * F) - math.log(10) * 2 * R * TKelvin / (n * F) * pH - 3 * R * TKelvin / (n * F) * math.log(Ce_OB * 1000 / FeMolarMass); # Equilibrium potential for dissolution reaction at O_B interface
        i0Fe2p_Fe3O4_OB = F * kboltzman * TKelvin / h * (math.exp(-1 * self.ActivationEFe2p_Fe3O4_OB  / (R * TKelvin )) * math.exp(-1 * Beta * n * F * EeFe2p_Fe3O4_OB / (R * TKelvin)));
        
        
        #Hydrogen consumption at O_B interface 2H+ + 2e- =  H2 
        #Standard DeltaGH2_Hp is zero at any temperature
        
        EeH2_Hp_OB = -1 * DeltaGH2_Hp_OB * 1000 / (n * F ) - math.log(10) * 2 * R * TKelvin / (n * F) * pH - R * TKelvin / (n * F) * math.log((CeH2_OB * 1000 / H2MolarMass) / kHenry); # At O_B interface
        
        i0H2_Hp_OB = F * kboltzman * TKelvin / h * math.exp(-1 * self.ActivationEH2_Hp_OB / (R * TKelvin)) * math.pow(CeHp/1000 , .67) * math.exp(-1 * Beta * n * F * EeH2_Hp_OB/(R * TKelvin));
        
        
        # Mixed Potential at the O_B interface
        Emixed = R * TKelvin / ( n * F) * math.log((i0H2_Hp_OB * math.exp(Beta * n * F * EeH2_Hp_OB / (R * TKelvin)) + i0Fe2p_Fe3O4_OB * math.exp(Beta * n * F * EeFe2p_Fe3O4_OB / (R * TKelvin))) / (i0H2_Hp_OB * math.exp(-1 * (1 - Beta) * n * F * EeH2_Hp_OB / (R * TKelvin)) + i0Fe2p_Fe3O4_OB * math.exp(-1 * (1 - Beta) * n * F * EeFe2p_Fe3O4_OB / (R * TKelvin))));
        
        
        # Current calculations at O_B interface
        iFe3O4_Fe2p_OB = i0Fe2p_Fe3O4_OB * (math.exp(Beta * n * F / (R * TKelvin) * (EeFe2p_Fe3O4_OB - Emixed)) - math.exp(-1 * (1 - Beta) * n * F / (R * TKelvin) * (EeFe2p_Fe3O4_OB - Emixed)));
        iH2_Hp_OB = i0H2_Hp_OB * (math.exp(Beta * n * F / (R * TKelvin) * (Emixed - EeH2_Hp_OB)) - math.exp(-1 * (1 - Beta) * n * F / (R * TKelvin) * (Emixed - EeH2_Hp_OB)));
        
        
        #Adjusting the concentration of iron species at O_B interface (in mol/l should be calculated back to g/cm^3)
        Ce_OBNew = ( math.exp((Emixed + DeltaGFe2p_Fe3O4_OB * 1000/ (n * F) + 2 * R * TKelvin * math.log(10) * pH / (n * F)) * (-1) * n * F / (3 * R * TKelvin)))* FeMolarMass / 1000
        CSat_OBNew = Ce_OBNew + CMP_Diffused
        
    
        #Adjusting the dissolution rate constant of magnetite
        ked_OB = kd_OB * math.exp(-1 * (1 - Beta) * n * F * (Emixed - EeFe2p_Fe3O4_OB) / (R * TKelvin))
    
        #H2 concentration
        
        CH2generated = 0.005659 * Corrosionrate * (7.32 - PhiOX )
        
        #Diffusivity of H2
        DH2 = 0.0000222 * TKelvin * math.exp(-12400 / (R * TKelvin))
        
        # H2 concentration at M_P interface assuming high mass transfer to the bulk
        CH2coolant = CH2coolant * 1e-6 * H2MolarMass * 101325 * RhoH2O / (8.314 * 298.15 * 1000)  # cm^3/kg water *1e-6 * H2MolarMass *P/(RT) * rowater/1000 = g H2/cm^3 H2O
        
        CeH2_MPNew = CH2coolant + (0.005659 * Corrosionrate * (7.32 - PhiOX) * (hH2 * DeltaOX * ChiOX / (RhoOX * (1 - PhiOX)) + (DH2 * PhiOX))) / (DH2 * hH2 * PhiOX)
        
        # H2 concentration at O_B interface
        
        CeH2_OBNew = CH2coolant + (0.005659 * Corrosionrate * (7.32 - PhiOX ) / hH2)
         
        return CSat_MPNew, CSat_OBNew, CeH2_MPNew, CeH2_OBNew, Ce_OBNew , Ce_MPNew, ked_OB, Emixed, Ecorr