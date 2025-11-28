# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

def xmass_to_xmol(xmass,mass_weight1,mass_weight2):
    xmol_converted = (xmass/mass_weight1)/(xmass/mass_weight1 + (1-xmass)/mass_weight2)
    return xmol_converted

def Patm_to_bar(Patm):
    return Patm * 1.01325

import iapws


N_0i = 3000
N_0j = 3000
Tol = pow(10,-7)

Patm_target = 100 # [atm]
P_target = Patm_to_bar(Patm_target) # [bar]
T_target = 273.15 + 132 # [K]
xmass_target = 0.99999 # [kg NH3/kg H2O]
MW_NH3 = 17.03026 # [g/mol] referência do iapws, no docuemnto de nh3h2o
MW_H2O = 18.015268 # [g/mol] referência do iapws, no docuemnto de nh3h2o
xmol_target = xmass_to_xmol(xmass_target,MW_NH3,MW_H2O) # [mol/mol]

#rho_guess_hand = 0.01 # [kg/m³]
rho_guess = 0.01
step_rho = 0.00001 # [kg/m³]

print("For:")
print("P = " + str(P_target) + " bar")
print("T = " + str(T_target) + " K")
print("x = " + str(xmass_target) + " (mass fraction of ammonia)")
print("")
print("we have:")


# Curva de rho é ascendente de quase zero até um certo ponto, onde começa a ser descendente.
# f(rho) = P0 - P(rho), por isso o differential é com sinal trocado em cima.

j = 1
while j <= N_0j:
    #print("j = " + str(j))
    rho_minusStep = rho_guess - step_rho
    #print("rho = "+ str(rho_minusStep)+".")
    rho_plusStep = rho_guess + step_rho
    #print("rho = "+ str(rho_plusStep)+".")
    if rho_minusStep <= 0 or rho_plusStep <= 0:
        print("Error! rho is negative (diverged) after " + str(j) + " iterations. Check the input conditions.")
        break
    else:
        P_minusStep = iapws.ammonia.H2ONH3()._prop(rho_minusStep, T_target, xmol_target).get("P") * 10 # [bar]
        P_plusStep = iapws.ammonia.H2ONH3()._prop(rho_plusStep, T_target, xmol_target).get("P") * 10 # [bar]
        #P_guess = iapws.ammonia.H2ONH3()._prop(rho_guess, T_target, xmol_target).get("P") * 10 # [bar]
    differential_P_rho = (-P_plusStep+P_minusStep)/(rho_plusStep-rho_minusStep)

    P_guess = iapws.ammonia.H2ONH3()._prop(rho_guess, T_target, xmol_target).get("P") * 10 # [bar]
    rho_now = rho_guess - (P_target - P_guess) / differential_P_rho
    if abs(rho_now - rho_guess) < Tol:
        print("rho = "+ str(rho_now)+", after "+str(j)+" iterations.")
        #print("P = " + str(iapws.ammonia.H2ONH3()._prop(rho_now, T_target, xmol_target).get("P") * 10) + " bar.")
        break
    else:
        rho_guess = rho_now
        j += 1
        #print("j+1 = " + str(j))
        
if j > N_0j:
    print("Method failed after " + str(j-1) + " iteration(s).")












