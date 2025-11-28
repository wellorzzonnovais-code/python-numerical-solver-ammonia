# -*- coding: utf-8 -*-
"""
This tool divides a functions into n parts and find a zero point between those.
Can be used as a first guess to others functions

# Restrictoins or To-do:
    1. Does not cover two zeros between the guesses range;
    2. Does not address outside functions, yet;
    3. Does not return more than just the range (could find the zero point itself).
    
"""

import useful_tools
import iapws


Tol_BS = pow(10,-3) #Tolerance bisection method

P_target = 11.0*useful_tools.pressure_bar_to_Pa # [Pa]
#print(P_target)
T_target = 273.15+130 # [K]
xmass_target = 0.42 # [kg NH3/kg H2O]
MW_NH3 = 17.03026 # [g/mol] referência do iapws, no docuemnto de nh3h2o
MW_H2O = 18.015268 # [g/mol] referência do iapws, no docuemnto de nh3h2o
xmol_target = useful_tools.xmass_to_xmol(xmass_target,MW_NH3,MW_H2O) # [mol/mol]

#rho_guess_hand = 0.01 # [kg/m³]
rho_guess = [20, 50] # [kg/m³]
#rho_guessBi = 500 # [kg/m³]
step_rho = 1.0 # [kg/m³]
step_diff_rho = 0.0001

n = 100


def divide_into_ranges(x_min, x_max, n_divisions):
    diff = (x_max - x_min)/n_divisions
    list = []
    i = 1
    while i <= n_divisions:
        min_range = x_min + diff*(i-1)
        max_range = min_range + diff
        list.append([min_range, max_range])
        i +=1
    return list


x_ranges = divide_into_ranges(rho_guess[0],rho_guess[1],n) # (7.427516, 7.427518, 100)
# for x in x_ranges:
#     print(x)

def function_P_rho(rho):
    try:
        function = (iapws.ammonia.H2ONH3()._prop(rho, T_target, xmol_target).get("P") * useful_tools.pressure_MPa_to_Pa) - P_target
        return function
    except:
        return print("Error in guess = " + str(rho) + ". \n")

function_ranges = []

for x_min, x_max in x_ranges:
    function_ranges.append([function_P_rho(x_min), function_P_rho(x_max)])
    

# print()
# print(function_ranges)

# for x in function_ranges:
#     print(x)

# def test_zeros(list_guesses, list_functions):

def test_function_zero(x_ranges, function_ranges, tol): 
    
    function_ranges_with_zero = []
    i = 0
    for f in function_ranges:
        i += 1
        if (f[0]*f[1] < 0):
            function_ranges_with_zero.append(f) # ATTENTION: cannot add multiple rho ranges, even using append list
            print()
            print("There is a zero in " + str(x_ranges[i-1]) + ", function between " + str(function_ranges_with_zero) + ".")
        elif abs(f[0]) <= tol:
            guess_with_zero = x_ranges[i-1][0]
            print()
            print("Guess with function equal to zero: " + str(guess_with_zero) + ".")
            break
        elif abs(f[1]) <= tol:
            guess_with_zero = x_ranges[i-1][1]
            print()
            print("Guess with function equal to zero: " + str(guess_with_zero) + ".")
            break
    
    if 'guess_with_zero' in locals():
        return print("\n" + str(guess_with_zero))
    else:
        return print("\n" + str(function_ranges_with_zero))


test_function_zero(x_ranges, function_ranges, Tol_BS)


