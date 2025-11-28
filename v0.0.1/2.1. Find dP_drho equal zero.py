# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import iapws


def xmass_to_xmol(xmass,mass_weight1,mass_weight2):
    xmol_converted = (xmass/mass_weight1)/(xmass/mass_weight1 + (1-xmass)/mass_weight2)
    return xmol_converted

def Patm_to_bar(Patm):
    return Patm * 1.01325

N_i = 3000
Tol_d = pow(10,-1)

Patm_target = 9 # [atm]
P_target = Patm_to_bar(Patm_target) #[bar]
T_target = 273.15 + 180 #[K]
xmass_target = 0.45 # [kg NH3/kg H2O]
MW_NH3 = 17.03052 # [g/mol]
MW_H2O = 18.01528 # [g/mol]
xmol_target = xmass_to_xmol(xmass_target,MW_NH3,MW_H2O) #[mol/mol]

rho_initial = 1
step_rho = 0.1 # [kg/m³]

i = 1
while i <= N_i:
    print("i = " + str(i))
    rho_minusStep = rho_initial - step_rho
    print("rho_minus = "+ str(rho_minusStep)+".")
    rho_plusStep = rho_initial + step_rho
    print("rho_plus = "+ str(rho_plusStep)+".")
    P_minusStep = iapws.ammonia.H2ONH3()._prop(rho_minusStep, T_target, xmol_target).get("P") * 10 # [bar]
    P_plusStep = iapws.ammonia.H2ONH3()._prop(rho_plusStep, T_target, xmol_target).get("P") * 10 # [bar]
    #P_guess = iapws.ammonia.H2ONH3()._prop(rho_guess, T_target, xmol_target).get("P") * 10 # [bar]
    differential_P_rho = (P_plusStep-P_minusStep)/(rho_plusStep-rho_minusStep)
    print("derivative = " + str(differential_P_rho))
    
    if differential_P_rho > 0 and abs(differential_P_rho) >= Tol_d:
        rho_initial += step_rho
        i += 1
    elif differential_P_rho < 0 and abs(differential_P_rho) >= Tol_d:
        rho_initial -= step_rho
        i += 1
    else:
        print("P = " + str(iapws.ammonia.H2ONH3()._prop(rho_initial, T_target, xmol_target).get("P") * 10) + " bar.")
        break
if i >= N_i:
    print("Iterations exceeded.")
else:
    print("Successful!")
    
    
    
    
    