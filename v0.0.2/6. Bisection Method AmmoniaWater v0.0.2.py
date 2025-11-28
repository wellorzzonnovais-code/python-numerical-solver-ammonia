# -*- coding: utf-8 -*-
"""

Find zero point using exclusively Bisection Method.

"""


import iapws
import useful_tools as ut/

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
        x_zero = find_zero_bisection_rho(x[0],x[1],tol,N)
        x_with_zero.append(x_zero)


# Print answers
if x_with_zero == []:
    print("Those intervals does not have a zero point. Try refine the range or change values for guesses.")
    # print("Function ranges: "+ str(function_ranges))
else:
    print("we have zero with guess = ")
    print(x_with_zero)

executionTime = (time.time() - start_time)
print("--- %s seconds ---" % executionTime)















