'''
gibbsEnergy.py
The is a class to calculate the Gibbs Free Engery at a given temperature
Author: Khatereh Maaaaaaaaaaa
Date: Aug 13, 2014
'''
import math

class gibbs(object):
    
    def __init__(self, Temperature):
        
        self.Temperature = Temperature + 273.15
        
    def setTemperature(self,Temperature):
        self.Temperature = Temperature
        
    def getDeltaGfFe(self):
            # Fe
        # Using the CP value at 298 and assuming constant CP over the temperatures of 298-573 K
        CP0Fe  = 24.98 #  from Beverskog
        GF0Fe = 0.0 # G of formation at standard temperature and pressure at 298.15K
        S0Fe = 27.28 # entropy at standard temperature and pressure from Beverskog
        HF0Fe = 0.0
        # Standard G formation of elements is  0 at any temperature by definition
        DeltaGfFe = 0.0
        return DeltaGfFe

    def getDeltaGfFe304(self,TKelvin = None):        
        # Fe3O4
        if TKelvin == None:
            TKelvin = self.Temperature
                
        CP0Fe3O4 = 150.73  #From Ziemniak
        GF0Fe3O4 = -1012.57# from Beverskog
        S0Fe3O4 = 146.14 # from Beverskog
        HF0Fe3O4 = -1118.38
        #G formation at T = DeltaGf = G0 + (CP0 - S0)*(TKelvin - 298.15) - T*CP0*Ln(T/298.15) assuming an average CP
        DeltaGfFe3O4 = GF0Fe3O4 + ((CP0Fe3O4 - S0Fe3O4) * (TKelvin - 298.15) - TKelvin * CP0Fe3O4 * math.log(TKelvin / 298.15))/1000
        #DeltaGfFe3O4 = HF0Fe3O4 + (-298 * CP0Fe3O4- S0Fe3O4 * TKelvin + CP0Fe3O4 * TKelvin- CP0Fe3O4 * TKelvin * (math.log(TKelvin/298) ))/1000
        return DeltaGfFe3O4

    def getDeltaGfH2O(self,TKelvin=None):
        # H2O
        # J/mol.K; Used the value of Cp at 25 C, small change over the temperature of 25-300 C
        # All values from Ziemniak
        if TKelvin == None:
            TKelvin = self.Temperature
        CP0H2O = 75.29 
        GF0H2O = -237.1
        S0H2O = 69.95
        HF0H2O = -285.85
        DeltaGfH2O = GF0H2O + ((CP0H2O - S0H2O) * (TKelvin - 298.15) - TKelvin * CP0H2O * math.log(TKelvin / 298.15))/1000
        #DeltaGfH2O = HF0H2O + (-298 * CP0H2O - S0H2O * TKelvin + CP0H2O * TKelvin- CP0H2O * TKelvin * (math.log(TKelvin/298) ))/1000
        return DeltaGfH2O
    
    def getDeltaGfH2(self):
        # H2 (g)
        # All values from Ziemniak
        CP0H2 = 28.84 
        GF0H2 = 0
        S0H2 = 130.68
        HF0H2 = 0
        DeltaGfH2 = 0
        return DeltaGfH2 
        
    def getDeltaGfHPlus(self,TKelvin=None):
        # H+ 
        # All values from Ziemniak
        if TKelvin == None:
            TKelvin = self.Temperature
        CP0Hp = -71.0 
        GF0Hp = 0
        S0Hp = -22.20
        HF0Hp = 0
        DeltaGfHp = GF0Hp + ((CP0Hp - S0Hp) * (TKelvin - 298.15) - TKelvin * CP0Hp * math.log(TKelvin / 298.15))/1000
        #DeltaGfHp = HF0Hp + (-298 * CP0Hp- S0Hp * TKelvin + CP0Hp * TKelvin- CP0Hp * TKelvin * (math.log(TKelvin/298) ))/1000
        return DeltaGfHp 
    
    #Fe(OH)2
    def getDeltaGfFeOH2(self,TKelvin=None):
        if TKelvin == None:
            TKelvin = self.Temperature
        CP0FeOH2 = 133.0 # From Tremaine
        GF0FeOH2 = -447  # From Beverskog
        S0FeOH2 = -29.20 # From Ziemniak
        HF0FeOH2 = -563.08
        DeltaGfFeOH2 = GF0FeOH2 + ((CP0FeOH2 - S0FeOH2) * (TKelvin - 298.15) - TKelvin * CP0FeOH2 * math.log(TKelvin / 298.15))/1000
        #DeltaGfFeOH2 = HF0FeOH2 + (-298 * CP0FeOH2- S0FeOH2 * TKelvin + CP0FeOH2 * TKelvin- CP0FeOH2 * TKelvin * (math.log(TKelvin/298) ))/1000
        return DeltaGfFeOH2
    
    def getDeltaGfFe2Plus(self,TKelvin=None):
        if TKelvin == None:
            TKelvin = self.Temperature
        #Fe2+
        CP0Fe2p = 90.0  #83
        GF0Fe2p = -88.92 #-91.88 # from Beverskog
        S0Fe2p = -107.1 #-105.6
        HF0Fe2p = -88.69
        DeltaGfFe2p = GF0Fe2p + ((CP0Fe2p - S0Fe2p) * (TKelvin - 298.15) - TKelvin * CP0Fe2p * math.log(TKelvin / 298.15))/1000
        #DeltaGfFe2p = HF0Fe2p + (-298 * CP0Fe2p- S0Fe2p * TKelvin + CP0Fe2p * TKelvin- CP0Fe2p * TKelvin * (math.log(TKelvin/298) ))/1000
        return DeltaGfFe2p
 
    
    def getDeltaGFe_Fe2(self,TKelvin=None):
        #Corrosion reaction: Fe + 2H2O = 2H+ +2e- + Fe(OH)2
        if TKelvin == None:
            TKelvin = self.Temperature
        DeltaGFe_Fe2 = 2*self.getDeltaGfHPlus(TKelvin) + self.getDeltaGfFeOH2(TKelvin) - self.getDeltaGfFe() - 2 * self.getDeltaGfH2O(TKelvin)
        return DeltaGFe_Fe2
        
    def getDeltaGFe2Plus_Fe304(self,TKelvin=None):
        #Magnetite precipitation and dissolution  3 Fe(OH)2 = Fe3O4 + 2H2O + 2H+ + 2e-
        if TKelvin == None:
            TKelvin = self.Temperature
        DeltaGFe2p_Fe3O4 = self.getDeltaGfFe304(TKelvin) + 2 * self.getDeltaGfH2O(TKelvin) + 2*self.getDeltaGfHPlus(TKelvin) - 3*self.getDeltaGfFeOH2(TKelvin) 
        return DeltaGFe2p_Fe3O4 
    
    def getDeltaGfH2_Hp(self,TKelvin=None):
        #Hydrogen Production or consumption:   H2 = 2H+ + 2e
        if TKelvin == None:
            TKelvin = self.Temperature
        DeltaGH2_Hp = 2 * self.getDeltaGfHPlus(TKelvin) - self.getDeltaGfH2() 
        #not getting similar delta cp and delta G at higher temperatures????????????????????????????????????????
        return DeltaGH2_Hp

    




