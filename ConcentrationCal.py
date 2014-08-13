import numpy as np
import matplotlib.pyplot as plt
import math


n = 2;
F = 96485;                   # C/mol
R = 8.314;                   # (j/mol.K)
TKelvin = T + 273.15; 
h = 6.62606957e-34           # Plank's constant J.s
kboltzman = 1.3806488e-23;   # Boltzman constant J/K
Beta =0.5;
DeltaGFe_Fe2_MP =
DeltaGFe2p_Fe3O4_MP = 
DeltaGFe2p_Fe3O4_OB =
DeltaGFe3O4_Fe2_OB = 
pH =  
kd_OB = 

def Concentrations(Ce_MP, CeH2_MP,Ce_OB, CeH2_OB):
# At the oxide-Bulk interface O_B
    # Corrosion reaction at the M_P interface; All the reaction are considered in the following direction: R -> O + ne-
    # Fe + 2H2O = 2H+ +2e- + Fe(OH)2
    
    EeFe_Fe2p_MP = -1*DeltaGFe_Fe2_MP/(n*F) - math.log(10)*2*R*TKelvin/(n*F)*pH + R*TKelvin/(n*F)*math.log(Ce_MP); # Equilibrium potential for the oxidation reaction at M_P interface
    i0Fe_Fe2p_MP = F*kboltzman *TKelvin/h*(math.exp(-1*DeltaGFe_Fe2_MP/(R*TKelvin))*Ce_MP*math.exp(Beta*n*F*EeFe_Fe2p_MP/(R*TKelvin)));
    
    #Precipitation of magnetite at M_P interface  3Fe(OH)2 = Fe3O4 + 2H2O + 2H+ + 2e-
    
    EeFe2p_Fe3O4_MP = -1*DeltaGFe2p_Fe3O4_MP/(n*F) - math.log(10)*2*R*TKelvin/(n*F)*pH - 3*R*TKelvin/(n*F)*math.log(Ce_MP);
    
    #Hydrogen production 2H+ + 2e- = H2 at M_P; The potential and current equations are written for the R => O + ne- ( H2 = 2H+ + 2e-)
    #Standard DeltaGH2_Hp is zero at any temperature
    
    EeH2_Hp_MP =  - math.log(10)*2*R*TKelvin/(n*F)*pH - R*TKelvin/(n*F)*math.log(CeH2_MP); # At M_P interface
    i0H2_Hp_MP = F*kboltzman *TKelvin/h*CeH2_MP*math.exp(Beta*n*F*EeH2_Hp_MP/(R*TKelvin));
    
    
    # Mixed Potential at the M_P interface
    Ecorr = R*Tkelvin/(Beta-(1-Beta)*n*F)*math.log((i0H2_Hp_MP*math.exp(Beta*n*F*EeH2_Hp_MP/(R*TKelvin))+i0Fe_Fe2p_MP*math.exp(Beta*n*F*EeFe_Fe2p_MP/(R*TKelvin)))/(i0H2_Hp_MP*math.exp(-1*(1-Beta)*n*F*EeH2_Hp_MP/(R*TKelvin))+i0Fe_Fe2p_MP*math.exp(-1*(1-Beta)*n*F*EeFe_Fe2p_MP/(R*TKelvin))));
    
    # Current calculations at M_P interface
    iFe_Fe2p_MP = i0Fe_Fe2p_MP *(math.exp(Beta*n*F/(R*TKelvin)*(Ecorr-EeFe_Fe2p_MP)) - math.exp(-1*(1-Beta)*n*F/(R*TKelvin)*(Ecorr-EeFe_Fe2p_MP)));
    iH2_Hp_MP = i0H2_Hp_MP *(math.exp(Beta*n*F/(R*TKelvin)*(Ecorr-EeH2_Hp_MP)) - math.exp(-1*(1-Beta)*n*F/(R*TKelvin)*(Ecorr-EeH2_Hp_MP)));
    
    Corrosionrate = iFe_Fe2p_MP * MWFe/ (n*F)
    
    #Adjusting Ce_MP for the Ecorr
    Ce_MPNew = math.exp((Ecorr + DeltaGFe_Fe2_MP/(n*F) + math.log(10)*2*R*TKelvin/(n*F)*pH)*(n*F)/(R*TKelvin))
    
    #Solubility at M-P: CMPsat assuming that the solubility at the M_P interface is summation of what is diffused at this interface and the equilibrium concentration considering the developed Ecorr
    
    CMP_Diffused = (0.476 * (1.101 + PhiOX) * ChiOX * DeltaOX) * Corrosionrate /(PhiOX * DFe * RhoOX * (1 - PhiOX))
    
    CSat_MPNew = CMP_Diffused +  Ce_MPNew;
    
    
    #Dissolution of magnetite at O_B interface Fe3O4 + 2H2O + 2H+ + 2e- = 3 Fe(OH)2 
    
    EeFe2p_Fe3O4_OB = -1*DeltaGFe2p_Fe3O4_OB/(n*F) - math.log(10)*2*R*TKelvin/(n*F)*pH - 3*R*TKelvin/(n*F)*math.log(Ce_OB); # Equilibrium potential for dissolution reaction at O_B interface
    i0Fe2p_Fe3O4_OB = F*kboltzman *TKelvin/h*(math.exp(-1*DeltaGFe2p_Fe3O4_OB/(R*TKelvin))*Ce_OB*math.exp(Beta*n*F*EeFe2p_Fe3O4_OB/(R*TKelvin)));
    
    
    #Hydrogen consumption at O_B interface H2 = 2H+ + 2e- 
    #Standard DeltaGH2_Hp is zero at any temperature
    
    EeH2_Hp_OB =  - math.log(10)*2*R*TKelvin/(n*F)*pH - R*TKelvin/(n*F)*math.log(CeH2_OB); # At O_B interface
    i0H2_Hp_OB = F*kboltzman *TKelvin/h*CeH2_OB*math.exp(Beta*n*F*EeH2_Hp_OB/(R*TKelvin));
    
    # Mixed Potential at the O_B interface
    Emixed = R*Tkelvin/(Beta-(1-Beta)*n*F)*math.log((i0H2_Hp_OB*math.exp(Beta*n*F*EeH2_Hp_OB/(R*TKelvin))+i0Fe3O4_Fe2p_OB*math.exp(Beta*n*F*EeFe2p_Fe3O4_OB/(R*TKelvin)))/(i0H2_Hp_OB*math.exp(-1*(1-Beta)*n*F*EeH2_Hp_OB/(R*TKelvin))+i0Fe3O4_Fe2p_OB*math.exp(-1*(1-Beta)*n*F*EeFe3O4_Fe2p_OB/(R*TKelvin))));
    
    
    # Current calculations at O_B interface
    iFe3O4_Fe2p_OB = i0Fe2p_Fe3O4_OB *(math.exp(Beta*n*F/(R*TKelvin)*(Emixed-EeFe2p_Fe3O4_OB)) - math.exp(-1*(1-Beta)*n*F/(R*TKelvin)*(Emixed-EeFe2p_Fe3O4_OB)));
    iH2_Hp_OB = i0H2_Hp_OB *(math.exp(Beta*n*F/(R*TKelvin)*(Emixed-EeH2_Hp_OB)) - math.exp(-1*(1-Beta)*n*F/(R*TKelvin)*(Emixed-EeH2_Hp_OB)));
    
    
    #Adjusting the concentration of iron species at O_B interface 
    Ce_OBNew = math.exp((Emixed + DeltaGFe3O4_Fe2_OB / (n * F) + 2*R*TKelvin*math.log(10)*pH/(n*F))*(-1)*n*F/(3*R*TKelvin))
    CSat_OBNew = Ce_OBNew;
    

#Adjusting the dissolution rate constant of magnetite
    ked_OB = kd_OB *math.exp((1-Beta)*n*F*Emixed/(R*TKelvin))

    #H2 concentration
    
    CH2generated = 0.005659 * Corrosionrate * (7.32 - PhiOX )
    
     #Diffusivity of H2
    DH2 = 0.0000222 * TKelvin *math.exp(-12400 / (R * T))
    
    # H2 concentration at M_P interface assuminghigh mass transfer to the bulk
    
    CeH2_MPNew = CH2coolant + (0.005659 * Corrosionrate * (7.32 - PhiOX )*(h*DeltaOX * ChiOX / (RhoOX * (1 - PhiOX)) + (DH2 * PhiOX))) / (DH2 * hH2 * PhiOX)
    
    # H2 concentration at O_B interface
    
    CeH2_OBNew = CH2coolant + (0.005659 * Corrosionrate * (7.32 - PhiOX ) / hH2)
     
    return  {CSat_MPNew, CSat_OBNew, CeH2_MPNew, CeH2_OBNew, Ce_OBNew, ked_OB, Emixed, Ecorr}