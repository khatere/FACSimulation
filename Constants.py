
'''
Constants.py
The is a class that assign values for the constants used in the FAC modeling plus initial values for the concentrations at different 
interfaces.
Author: Khatereh Mohajery
Date: Aug 15, 2014
'''
import numpy as np

class Const(object):
    
    def __init__(self): 
        self.CSatMP = 1.16E-08            # CSat_MP: Initial value for solubility of dissolved iron at metal-passivation layer interface    g/cm^3
        self.CSatPO = 1.16E-07            # CSat_ PO: solubility of dissolved iron at passivation-oxide layer interface    g/cm^3
        self.kd0 = 2.00E-02              # kd: Dissolution Rate Constant of Magnetite    cm/s
        self.km0 = 5.41E-01               # km:  Static Mass Transfer Coefficient    cm/s 
        self.f0 = 1.732050808             # f: Surface Area Factor (F* in thesis)    no unit
        self.S0 = 1.10E+00                # S: supersaturation factor of small crystallites  no unit  
        self.CB0 = 4.640E-8                 # CB: concentration of iron in bulk solution    g/cm^3
        self.DeltaP0 = 9.78E-06           # DeltaP: passivating layer thickness    g/cm^2 (10 nm thickness*density of passivating layer)
        self.ChiP0 = 1.2                  # ChiP: tortuosity of passivating layer    
        self.PhiP0 = 0.199989634          # PhiP: porosity of passivating layer    
        self.DFe0 = 0.0002735             # DFe: diffusivity of iron in water    cm^2/s
        self.RhoP0 = 5.15                 # RhoP: Density of passivating layer    g/cm^3
        self.MFe0 = 0.723569              # MFe: mass fraction of iron in oxide     
        self.MP0 = 1.0                    # MP: mass fraction of iron in passivating layer    
        self.RhoFe0 = 7.9                 # RhoFe: Density of Iron    g/cm^3
        self.PhiOX0 = 0.25                # PhiOX: porosity of oxide    
        self.ChiOX0 = 1.8                 # ChiOX: tortuosity of oxide 
        self.RhoOX0 = 5.2                 # RhoOX: Density of Oxide    g/cm^3
        self.CPO0 = 1.16e-8#1.16E-07              # CPO: concentration at passivating-oxide interface    g/cm^3     
        self.DCORs0 = 0.0                 # DCORs: initial concentration at oxide bulk   
        self.DOXs0 = 0.0                 # DOXs: initial oxide layer thickness 
        self.t0 = 2000.0                   # t: time interval for simulation    seconds
        self.Tstep0 = 7000.0                # Tstep: number of steps    
        self.DP0 = 0.000273486            # DP: diffusivity of iron in passivating layer    cm^2/s
        self.CSatOB = 1.16E-7            # CSat_OB: Solubility of dissolved iron at oxide-bulk solution interface    g/cm^3    
        self.CeOB = 6.96e-8#1.0E-10              # Ce_OB: Concentration of dissolved iron at Oxide-Bulk interface    g/cm^3
        self.CrContent0 = 0.001#0.019           # CrContent: Chromium Content    %
        self.Paverage0 = 1.00E-06         # Paverage: minimum size of oxide particle    cm
        self.SDP0 = 1.00E-05              # SDP: maximum size of oxide particle    cm
        self.TS0 = 1.6                    # TS: Supersaturation factor  
        self.FF0 = 1.0                    # FF: Fricton Factor   
        self.U0 = 9.51e2 #2.5e2#188.9965                # U: Linear Velocity of Coolant    cm/s   
        self.Tphi0 = 0.25                 # Tphi: porosity of oxide (for spalling calculation)  
        self.TKd0 = 0.2                   # TKd: dissolution rate constant (for spalling calculation)   
        self.TFSA0 = 1.7321               # TFSA: Surface Area Factor (F* in thesis)    
        self.TFeSAT0 = 1.16E-07           # TFeSAT: solubility of dissolved iron at oxide-bulk layer interface    
        self.Tfe0 = 0.00E+00              # Tfe: bulk iron concentration in solution
        self.RhoH2O0 = 0.92579            # RhoH2O: Density of Water    g/cm^3  
        self.ksp0 = 15000.0              # Ksp: Spalling time constant   
        self.Cell0 = 200.0                # Cell: Number of nodes
        self.Temperature0 = 140.0           # Temperature: Temperature    C    
        self.ConcC4H9ONTotal0 = 0.000516  # ConcC4H9ONTotal: Morpholine concentration     mol/L   
        self.n0 = 2.0                     # n: number of electrons transferd in a reaction
        self.F0 = 96485.0                 # F: Faraday constant    C/mol (s.A/mol) 
        self.R0 = 8.314                   # R: Gas constant    j/mol.K
        self.h0 = 6.63E-34                # h: Plank's constant     J.s
        self.kboltzman0 = 1.38E-23        # kboltzman: Boltzman constant    J/K     
        self.Beta0 = 0.5                  # Beta: symmetry factor or charge transfer coefficient 
        self.FeMolarMass0 = 90.0         # Fe(OH)2 molar mass g/mol
        self.CH2coolant0 = 10.0             # Bulk concentration of hydrogen cm^3/Kg water
        self.H2MolarMass0 = 2.016         # Hydrogen molar mass g/mol
        self.EcorrInitial0 =  -0.55      #Corrosion potential from experimental data V
        self.RhoSteel0 =  7.86            # Density of Steel g/cm^3
        self.CorrosionRateExp0 = 300.0      # Experimental ( Actual measured) Corrosion rate microm/year
        self.iH2_Hp_OB0 = 1.18e-12    # The current at the oxide-solution interface, is used to get the current density and activation energy of magnetite dissolution at this interface, should be adjusted and taken from similar works and data
        self.EmixedInitial0 = -0.50
    def CSat_MP(self, CSatMP = None):
        if CSatMP == None: 
            CSatMP = self.CSatMP 
        return CSatMP
    def CSat_PO(self, CSatPO = None):
        if CSatPO == None: 
            CSatPO = self.CSatPO 
        return CSatPO  
    def kd(self, kd0 = None):
        if kd0 == None: 
            kd0 = self.kd0 
        return kd0 
    def km(self, km0 = None):
        if km0 == None: 
            km0 = self.km0 
        return km0     
    def f(self, f0 = None):
        if f0 == None: 
            f0 = self.f0 
        return f0 
    def S(self, S0 = None):
        if S0 == None: 
            S0 = self.S0 
        return S0 
    def CB(self, CB0 = None):
        if CB0 == None: 
            CB0 = self.CB0 
        return CB0        
    def DeltaP(self, DeltaP0 = None):
        if DeltaP0 == None: 
            DeltaP0 = self.DeltaP0 
        return DeltaP0
    def ChiP(self, ChiP0 = None):
        if ChiP0 == None: 
            ChiP0 = self.ChiP0 
        return ChiP0
    def PhiP(self, PhiP0 = None):
        if PhiP0 == None: 
            PhiP0 = self.PhiP0 
        return PhiP0
    def DFe(self, DFe0 = None):
        if DFe0 == None: 
            DFe0 = self.DFe0 
        return DFe0
    def RhoP(self, RhoP0 = None):
        if RhoP0 == None: 
            RhoP0 = self.RhoP0 
        return RhoP0
    def MFe(self, MFe0 = None):
        if MFe0 == None: 
            MFe0 = self.MFe0 
        return MFe0
    def MP(self, MP0 = None):
        if MP0 == None: 
            MP0 = self.MP0 
        return MP0
    def RhoFe(self, RhoFe0 = None):
        if RhoFe0 == None: 
            RhoFe0 = self.RhoFe0 
        return RhoFe0
    def PhiOX(self, PhiOX0 = None):
        if PhiOX0 == None: 
            PhiOX0 = self.PhiOX0 
        return PhiOX0
    def ChiOX(self, ChiOX0 = None):
        if ChiOX0 == None: 
            ChiOX0 = self.ChiOX0 
        return ChiOX0
    def RhoOX(self, RhoOX0 = None):
        if RhoOX0 == None: 
            RhoOX0 = self.RhoOX0 
        return RhoOX0
    def CPO(self, CPO0 = None):
        if CPO0 == None: 
            CPO0 = self.CPO0 
        return CPO0
    def DCORs(self, DCORs0 = None):
        if DCORs0 == None: 
            DCORs0 = self.DCORs0 
        return S0
    def DOXs(self, DOXs0 = None):
        if DOXs0 == None: 
            DOXs0 = self.DOXs0 
        return DOXs0
    def t(self, t0 = None):
        if t0 == None: 
            t0 = self.t0 
        return t0
    def Tstep(self, Tstep0 = None):
        if Tstep0 == None: 
            Tstep0 = self.Tstep0 
        return Tstep0 
    def DP(self, DP0 = None):
        if DP0 == None: 
            DP0 = self.DP0 
        return DP0 
    def CSat_OB(self, CSatOB = None):
        if CSatOB == None: 
            CSatOB = self.CSatOB 
        return CSatOB 
    def Ce_OB(self, CeOB = None):
        if CeOB == None: 
            CeOB = self.CeOB 
        return CeOB 
    def CrContent(self, CrContent0 = None):
        if CrContent0 == None: 
            CrContent0 = self.CrContent0 
        return CrContent0 
    def Paverage(self, Paverage0 = None):
        if Paverage0 == None: 
            Paverage0 = self.Paverage0 
        return Paverage0 
    def SDP(self, SDP0 = None):
        if SDP0 == None: 
            SDP0 = self.SDP0 
        return SDP0 
    def TS(self, TS0 = None):
        if TS0 == None: 
            TS0 = self.TS0 
        return TS0 
    def FF(self, FF0 = None):
        if FF0 == None: 
            FF0 = self.FF0 
        return FF0 
    def U(self, U0 = None):
        if U0 == None: 
            U0 = self.U0 
        return U0 
    def Tphi(self, Tphi0 = None):
        if Tphi0 == None: 
            Tphi0 = self.Tphi0 
        return Tphi0 
    def TKd(self, TKd0 = None):
        if TKd0 == None: 
            TKd0 = self.TKd0 
        return TKd0 
    def TFSA(self, TFSA0 = None):
        if TFSA0 == None: 
            TFSA0 = self.TFSA0 
        return TFSA0 
    def TFeSAT(self, TFeSAT0 = None):
        if TFeSAT0 == None: 
            TFeSAT0 = self.TFeSAT0 
        return TFeSAT0 
    def Tfe(self, Tfe0 = None):
        if Tfe0 == None: 
            Tfe0 = self.Tfe0 
        return Tfe0 
    def RhoH2O(self, RhoH2O0 = None):
        if RhoH2O0 == None: 
            RhoH2O0 = self.RhoH2O0 
        return RhoH2O0 
    def ksp(self, ksp0 = None):
        if ksp0 == None: 
            ksp0 = self.ksp0 
        return ksp0 
    def Cell(self, Cell0 = None):
        if Cell0 == None: 
            Cell0 = self.Cell0 
        return Cell0 
    def Temperature(self, Temperature0 = None):
        if Temperature0 == None: 
            Temperature0 = self.Temperature0 
        return Temperature0 
    def ConcC4H9ONTotal(self, ConcC4H9ONTotal0 = None):
        if ConcC4H9ONTotal0 == None: 
            ConcC4H9ONTotal0 = self.ConcC4H9ONTotal0 
        return ConcC4H9ONTotal0 
    def n(self, n0 = None):
        if n0 == None: 
            n0 = self.n0 
        return n0 
    def F(self, F0 = None):
        if F0 == None: 
            F0 = self.F0 
        return F0 
    def R(self, R0 = None):
        if R0 == None: 
            R0 = self.R0 
        return R0 
    def h(self, h0 = None):
        if h0 == None: 
            h0 = self.h0 
        return h0 
    def kboltzman(self, kboltzman0 = None):
        if kboltzman0 == None: 
            kboltzman0 = self.kboltzman0 
        return kboltzman0 
    def Beta(self, Beta0 = None):
        if Beta0 == None: 
            Beta0 = self.Beta0 
        return Beta0 
    def FeMolarMass(self, FeMolarMass0 = None):
        if FeMolarMass0 == None: 
            FeMolarMass0 = self.FeMolarMass0 
        return FeMolarMass0 
    def CH2coolant(self, CH2coolant0 = None):
        if CH2coolant0 == None: 
            CH2coolant0 = self.CH2coolant0 
        return CH2coolant0 
    def H2MolarMass(self, H2MolarMass0 = None):
        if H2MolarMass0 == None: 
            H2MolarMass0 = self.H2MolarMass0 
        return H2MolarMass0  
    def EcorrInitial(self, EcorrInitial0 = None):
        if EcorrInitial0 == None: 
            EcorrInitial0 = self.EcorrInitial0 
        return EcorrInitial0 
    def EmixedInitial(self, EmixedInitial0 = None):
        if EmixedInitial0 == None: 
            EmixedInitial0 = self.EmixedInitial0 
        return EmixedInitial0 
    def RhoSteel(self, RhoSteel0 = None):
        if RhoSteel0 == None: 
            RhoSteel0 = self.RhoSteel0 
        return RhoSteel0 
    def CorrosionRateExp(self, CorrosionRateExp0 = None):
        if CorrosionRateExp0 == None: 
            CorrosionRateExp0 = self.CorrosionRateExp0 
        return CorrosionRateExp0  
    def iH2_Hp_OB(self, iH2_Hp_OB0 = None):
        if iH2_Hp_OB0 == None: 
            iH2_Hp_OB0 = self.iH2_Hp_OB0 
        return iH2_Hp_OB0  
    
    