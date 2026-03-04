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


N_guess = 3000
N_0i = 3000
Tol_BS = pow(10,-7) #Tolerance bisection method

Patm_target = 100 # [atm]
P_target = Patm_to_bar(Patm_target) # [bar]
T_target = 273.15 + 132 # [K]
xmass_target = 0.99999 # [kg NH3/kg H2O]
MW_NH3 = 17.03026 # [g/mol] referência do iapws, no docuemnto de nh3h2o
MW_H2O = 18.015268 # [g/mol] referência do iapws, no docuemnto de nh3h2o
xmol_target = xmass_to_xmol(xmass_target,MW_NH3,MW_H2O) # [mol/mol]

#rho_guess_hand = 0.01 # [kg/m³]
rho_guessA = 0.1 # [kg/m³]
rho_guessB = 50.0 # [kg/m³]
step_rho = 1 # [kg/m³]
step_diff_rho = 0.0001

print("For:")
print("P = " + str(P_target) + " bar")
print("T = " + str(T_target) + " K")
print("x = " + str(xmass_target) + " (mass fraction of ammonia)")
print("")
print("we have:")


# Curva de rho é ascendente de quase zero até um certo ponto, onde começa a ser descendente.
# f(rho) = P0 - P(rho), por isso o differential é com sinal trocado em cima.
i1 = 1
while i1 < N_guess:
    #print("i = " + str(i))
    rho_minusStep1 = rho_guessA - step_diff_rho
    #print("rho_minus = "+ str(rho_minusStep)+".")
    rho_plusStep1 = rho_guessA + step_diff_rho
    #print("rho_plus = "+ str(rho_plusStep)+".")
    P_minusStep1 = iapws.ammonia.H2ONH3()._prop(rho_minusStep1, T_target, xmol_target).get("P") # [bar]
    P_plusStep1 = iapws.ammonia.H2ONH3()._prop(rho_plusStep1, T_target, xmol_target).get("P") # [bar]
    #P_guess = iapws.ammonia.H2ONH3()._prop(rho_guess, T_target, xmol_target).get("P") * 10 # [bar]
    differential_P_rho1 = (P_plusStep1-P_minusStep1)/(rho_plusStep1-rho_minusStep1)
    #print("derivative = " + str(differential_P_rho))
    if differential_P_rho1 <= 0:
        rho_guessA += step_rho
        i1 += 1
    else:
        print("guess one = " + str(rho_guessA) + " kg/m³.")
        print("guess one is valid after " + str(i1-1) + " iterations.")
        break

i2 = 1
while i2 < N_guess:
    #print("i2 = " + str(i2))
    rho_minusStep2 = rho_guessB - step_diff_rho
    #print("rho_minus = "+ str(rho_minusStep)+".")
    rho_plusStep2 = rho_guessB + step_diff_rho
    #print("rho_plus = "+ str(rho_plusStep)+".")
    P_minusStep2 = iapws.ammonia.H2ONH3()._prop(rho_minusStep2, T_target, xmol_target).get("P") # [bar]
    P_plusStep2 = iapws.ammonia.H2ONH3()._prop(rho_plusStep2, T_target, xmol_target).get("P") # [bar]
    #P_guess = iapws.ammonia.H2ONH3()._prop(rho_guess, T_target, xmol_target).get("P") * 10 # [bar]
    differential_P_rho2 = (P_plusStep2-P_minusStep2)/(rho_plusStep2-rho_minusStep2)
    #print("derivative 2 = " + str(differential_P_rho2))
    if differential_P_rho2 <= 0:
        rho_guessB -= step_rho
        i2 += 1
    elif differential_P_rho2 > pow(10,-3):
        rho_guessB += step_rho
        i2 += 1
    else:
        print("guess two = " + str(rho_guessB) + " kg/m³.")
        print("guess two is valid after " + str(i2-1) + " iterations.")
        break

# Bisection Method
#-----------------------------------
#-----------------------------------
#-----------------------------------


i = 1
while i <= N_0i:
    avg_rho_guess = rho_guessA + (rho_guessB - rho_guessA)/2
    function_guess1 = P_target - iapws.ammonia.H2ONH3()._prop(rho_guessA, T_target, xmol_target).get("P") * 10 # [bar]
    function_avg_guess = P_target - iapws.ammonia.H2ONH3()._prop(avg_rho_guess, T_target, xmol_target).get("P") * 10 # [bar]
    if function_avg_guess == 0 or abs((rho_guessB - rho_guessA)/2) < Tol_BS:
        break
    elif function_guess1*function_avg_guess > 0:
        rho_guessA = avg_rho_guess
    else:
        rho_guessB = avg_rho_guess
    i += 1

if i > N_0i:
    print("Method failed after " + str(i-1) + " iteration(s).")
else:
    print("rho_average = "+ str(avg_rho_guess)+", after "+str(i)+" iterations.")
    print("Delta_P = " + str(function_avg_guess))
    print("")
    #print("Take care! Guess two could be under estimated.")

#-----------------------------------
#-----------------------------------
#-----------------------------------

"""















