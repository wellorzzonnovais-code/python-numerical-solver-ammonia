//# -*- coding: utf-8 -*-
"""

Find zero point using exclusively Newton-Raphson Method.

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
tol = pow(10,-7) # Tolerance to zero
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
x_zero_ranges = ut.x_with_zero_ranges(x_ranges, function_ranges, tol)
# print(x_zero_ranges)
# print()

# print((x_zero_ranges[0][0]+x_zero_ranges[0][1])/2)


def listoflists_into_avgnumbers(list_of_guesses_or_ranges):
    list_avg = []
    j=0
    for n in list_of_guesses_or_ranges:
        # print(list_of_guesses_or_ranges)
        if not(isinstance(n, list)):
            list_avg.append(n)
        else:
            # print(j)
            list_avg.append((n[0]+n[1])/2)
    return list_avg
    

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

     

x_guess = listoflists_into_avgnumbers(x_zero_ranges)
# print(x_guess)
x_zeros = []
for x in x_guess:
    x_zero = find_zero_newtonraphson_rho(x, tol, N)
    x_zeros.append(x_zero)
    
# Print answers
if x_zeros == []:
    print("Those guesses does not have a zero point or diverged. Try refine the range or change values for guesses.")
else:
    print("we have zero with guess = ")
    print(x_zeros)


executionTime = (time.time() - start_time)
print("--- %s seconds ---" % executionTime)








