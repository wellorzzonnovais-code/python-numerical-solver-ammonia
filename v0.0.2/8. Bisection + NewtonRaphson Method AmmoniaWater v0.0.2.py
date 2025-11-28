# -*- coding: utf-8 -*-
"""

Find zero point using Bisection (first) and Newton-Raphson (second) Method.

"""


import iapws
import useful_tools as ut

import time
start_time = time.time() # Measuring above 0.1 seconds


P_target = 15.0*ut.pressure_bar_to_Pa # [Pa]
T_target = 273.15+120 # [K]
xmass_target = 0.42 # [kg NH3/kg H2O]
MW_NH3 = 17.03026 # [g/mol] referência do iapws, no docuemnto de nh3h2o
MW_H2O = 18.015268 # [g/mol] referência do iapws, no docuemnto de nh3h2o
xmol_target = ut.xmass_to_xmol(xmass_target,MW_NH3,MW_H2O) # [mol/mol]

rho_guess = [0.1,25.1] # [kg/m³] [6.369055814002477,20.536754868000003] zeros of 11 bar, 132°C and xmass = 0.42

n = 10 # number to divide guess range
tol_BS = pow(10,-1) # Tolerance to zero for Bisection Method
tol_NR = pow(10,-7) # Tolerance to zero for Newton-Raphson Method
N = 3000 # limit for numbers of iterations


print("For:")
print("P = " + str(P_target) + " Pa")
print("T = " + str(T_target) + " K")
print("x_mass = " + str(xmass_target) + " (mass fraction of ammonia)")
print("rho guess = " + str(rho_guess) + " kg/m³")
print("")



def function_P_rho(rho):
    try:
        function = (iapws.ammonia.H2ONH3()._prop(rho, T_target, xmol_target).get("P") * ut.pressure_MPa_to_Pa) - P_target
        return function
    except:
        return print("Error in guess = " + str(rho) + ". \n")


x_ranges = ut.divide_into_ranges(rho_guess[0],rho_guess[1],n)
# print(x_ranges)
function_ranges = []
for x_min, x_max in x_ranges:
    function_ranges.append([function_P_rho(x_min), function_P_rho(x_max)])
# print(function_ranges)
x_zero_ranges = ut.x_with_zero_ranges(x_ranges, function_ranges, tol_NR)
# print(ut.x_with_zero_ranges(x_ranges, function_ranges, tol))



# Bisection Method
#-----------------------------------

def find_zero_bisection_rho(guess_min, guess_max, tolerance, limit_iterations):
    i = 1
    while i <= limit_iterations:
        guess_avg = guess_min + (guess_max - guess_min)/2
        function_guess_min = function_P_rho(guess_min) # [bar]
        function_guess_avg = function_P_rho(guess_avg) # [bar]
        if function_guess_avg == 0 or abs((guess_min - guess_max)/2) < tolerance:
            return guess_avg
        elif function_guess_min*function_guess_avg > 0:
            guess_min = guess_avg
        else:
            guess_max = guess_avg
        i += 1
        
    if i > N:
        return print("Method failed after " + str(i-1) + " iteration(s).")


# Check if I have a range or zero points, and apply bisection method just within ranges
x_with_zero = []
for x in x_zero_ranges:
    if not(isinstance(x, list)):
        x_with_zero.append(x)
    else:
        x_zero = find_zero_bisection_rho(x[0],x[1],tol_BS,N)
        x_with_zero.append(x_zero)


def find_zero_newtonraphson_rho(guess, tolerance, limit_iterations):
    i = 1
    while i <= N:
        guess_minusStep = guess*(1-tolerance) # tolerance used to estimate differential curve
        guess_plusStep = guess*(1+tolerance)
        if guess_minusStep <= 0 or guess_plusStep <= 0:
            return print("Error! Guess diverged after " + str(i) + " iterations. Check the input conditions.")
        
        differential_function = (function_P_rho(guess_plusStep)-function_P_rho(guess_minusStep))/(guess_plusStep-guess_minusStep)
        guess_now = guess-function_P_rho(guess)/differential_function
        if abs(guess_now-guess) < tolerance:
            return guess_now #print("guess = "+ str(guess_now)+", after "+str(i)+" iterations.")
        else:
            guess = guess_now
            i += 1
            
    if i > N:
        return print("Method failed after " + str(i-1) + " iteration(s).")

     
x_zeros = []
for x in x_with_zero:
    x_zero = find_zero_newtonraphson_rho(x, tol_NR, N)
    x_zeros.append(x_zero)
    
# Print answers
if x_zeros == []:
    print("Those guesses does not have a zero point or diverged. Try refine the range or change values for guesses.")
else:
    print("we have zero with guess = ")
    print(x_zeros)


executionTime = (time.time() - start_time)
print("--- %s seconds ---" % executionTime)
















"""
def xmass_to_xmol(xmass,mass_weight1,mass_weight2):
    xmol_converted = (xmass/mass_weight1)/(xmass/mass_weight1 + (1-xmass)/mass_weight2)
    return xmol_converted

def Patm_to_bar(Patm):
    return Patm * 1.01325

import iapws

N_guess = 3000
N_0i = 3000
N_0j = 3000
Tol_BS = pow(10,-0) #Tolerance bisection method
Tol_NR = pow(10,-7) #Tolerance Newton-Raphson method

Patm_target = 100 # [atm]
P_target = Patm_to_bar(Patm_target) # [bar]
T_target = 273.15 + 132 # [K]
xmass_target = 0.99999 # [kg NH3/kg H2O]
MW_NH3 = 17.03026 # [g/mol] referência do iapws, no docuemnto de nh3h2o
MW_H2O = 18.015268 # [g/mol] referência do iapws, no docuemnto de nh3h2o
xmol_target = xmass_to_xmol(xmass_target,MW_NH3,MW_H2O) # [mol/mol]

#rho_guess_hand = 0.01 # [kg/m³]
rho_guess1 = 0.1 # [kg/m³]
rho_guess2 = 50.0 # [kg/m³]
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
    rho_minusStep1 = rho_guess1 - step_diff_rho
    #print("rho_minus = "+ str(rho_minusStep)+".")
    rho_plusStep1 = rho_guess1 + step_diff_rho
    #print("rho_plus = "+ str(rho_plusStep)+".")
    P_minusStep1 = iapws.ammonia.H2ONH3()._prop(rho_minusStep1, T_target, xmol_target).get("P") # [bar]
    P_plusStep1 = iapws.ammonia.H2ONH3()._prop(rho_plusStep1, T_target, xmol_target).get("P") # [bar]
    #P_guess = iapws.ammonia.H2ONH3()._prop(rho_guess, T_target, xmol_target).get("P") * 10 # [bar]
    differential_P_rho1 = (P_plusStep1-P_minusStep1)/(rho_plusStep1-rho_minusStep1)
    #print("derivative = " + str(differential_P_rho))
    if differential_P_rho1 <= 0:
        rho_guess1 += step_rho
        i1 += 1
    else:
        print("guess one is valid after " + str(i1-1) + " iterations.")
        break

i2 = 1
while i2 < N_guess:
    #print("i = " + str(i))
    rho_minusStep2 = rho_guess2 - step_diff_rho
    #print("rho_minus = "+ str(rho_minusStep)+".")
    rho_plusStep2 = rho_guess2 + step_diff_rho
    #print("rho_plus = "+ str(rho_plusStep)+".")
    P_minusStep2 = iapws.ammonia.H2ONH3()._prop(rho_minusStep2, T_target, xmol_target).get("P") # [bar]
    P_plusStep2 = iapws.ammonia.H2ONH3()._prop(rho_plusStep2, T_target, xmol_target).get("P") # [bar]
    #P_guess = iapws.ammonia.H2ONH3()._prop(rho_guess, T_target, xmol_target).get("P") * 10 # [bar]
    differential_P_rho2 = (P_plusStep2-P_minusStep2)/(rho_plusStep2-rho_minusStep2)
    #print("derivative = " + str(differential_P_rho))
    if differential_P_rho2 <= 0:
        rho_guess2 -= step_rho
        i2 += 1
    else:
        print("guess two is valid after " + str(i2-1) + " iterations.")
        break

# Bisection Method
#-----------------------------------
#-----------------------------------
#-----------------------------------


i = 1
while i <= N_0i:
    avg_rho_guess = rho_guess1 + (rho_guess2 - rho_guess1)/2
    function_guess1 = P_target - iapws.ammonia.H2ONH3()._prop(rho_guess1, T_target, xmol_target).get("P") * 10 # [bar]
    function_avg_guess = P_target - iapws.ammonia.H2ONH3()._prop(avg_rho_guess, T_target, xmol_target).get("P") * 10 # [bar]
    if function_avg_guess == 0 or abs((rho_guess2 - rho_guess1)/2) < Tol_BS:
        break
    elif function_guess1*function_avg_guess > 0:
        rho_guess1 = avg_rho_guess
    else:
        rho_guess2 = avg_rho_guess
    i += 1

if i > N_0i:
    print("Method failed after " + str(i-1) + " iteration(s).")
else:
    print("rho_average = "+ str(avg_rho_guess)+", after "+str(i)+" iterations.")
    print("Delta_P = " + str(function_avg_guess))

#-----------------------------------
#-----------------------------------
#-----------------------------------

rho_guess = avg_rho_guess

j = 1
while j <= N_0j:
    #print("j = " + str(j))
    rho_minusStep = rho_guess - step_rho
    #print("rho = "+ str(rho_minusStep)+".")
    #rho_plusStep = rho_guess + step_rho
    #print("rho = "+ str(rho_plusStep)+".")
    if rho_minusStep <= 0 or rho_guess <= 0:
        print("Error! rho is negative (diverged) after " + str(j) + " iterations. Check the input conditions.")
        break
    else:
        P_minusStep = iapws.ammonia.H2ONH3()._prop(rho_minusStep, T_target, xmol_target).get("P") * 10 # [bar]
        #P_plusStep = iapws.ammonia.H2ONH3()._prop(rho_plusStep, T_target, xmol_target).get("P") * 10 # [bar]
        P_guess = iapws.ammonia.H2ONH3()._prop(rho_guess, T_target, xmol_target).get("P") * 10 # [bar]
    differential_P_rho = (-P_guess+P_minusStep)/(rho_guess-rho_minusStep)

    P_guess = iapws.ammonia.H2ONH3()._prop(rho_guess, T_target, xmol_target).get("P") * 10 # [bar]
    rho_now = rho_guess - (P_target - P_guess) / differential_P_rho
    if abs(rho_now - rho_guess) < Tol_NR:
        print("rho = "+ str(rho_now)+", after "+str(j)+" iterations.")
        #print("P = " + str(iapws.ammonia.H2ONH3()._prop(rho_now, T_target, xmol_target).get("P") * 10) + " bar.")
        break
    else:
        rho_guess = rho_now
        j += 1
        #print("j+1 = " + str(j))
        
if j > N_0j:
    print("Method failed after " + str(j-1) + " iteration(s).")

"""