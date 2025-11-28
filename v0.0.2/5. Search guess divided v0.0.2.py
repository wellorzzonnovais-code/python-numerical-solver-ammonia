# -*- coding: utf-8 -*-
"""
This tool divides a functions into n parts and find a zero point between those.
Can be used as a first guess to others functions

# Restrictoins or To-do:
    2. Does not address outside functions, yet;
    3. Does not return more than just the range (could find the zero point itself).
    
"""

import useful_tools
import iapws


P_target = 15.0*useful_tools.pressure_bar_to_Pa # [Pa]
T_target = 273.15+120 # [K]
xmass_target = 0.42 # [kg NH3/kg H2O]
MW_NH3 = 17.03026 # [g/mol] referência do iapws, no docuemnto de nh3h2o
MW_H2O = 18.015268 # [g/mol] referência do iapws, no docuemnto de nh3h2o
xmol_target = useful_tools.xmass_to_xmol(xmass_target,MW_NH3,MW_H2O) # [mol/mol]

rho_guess = [0.1,25.1] # [kg/m³] [6.369055814002477,20.536754868000003] zeros of 11 bar, 132°C and xmass = 0.42

n = 100 # number to divide guess range
tol = pow(10,-7) # Tolerance to zero


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


def function_P_rho(rho):
    try:
        function = (iapws.ammonia.H2ONH3()._prop(rho, T_target, xmol_target).get("P") * useful_tools.pressure_MPa_to_Pa) - P_target
        return function
    except:
        return print("Error in guess = " + str(rho) + ". \n")


def function_zeros_range(x_ranges, function_ranges, tolerance): 
    function_ranges_with_zero = []
    guess_with_zero = []
    i = 0
    for f in function_ranges:
        i += 1
        if (f[0]*f[1] < 0):
            function_ranges_with_zero.append(f) # ATTENTION: cannot add multiple rho ranges, even using append list
            print()
            print("There is a zero in " + str(x_ranges[i-1]))
            print("with function between " + str(f) + ".")
        elif abs(f[0]) <= tol:
            guess_with_zero.append(x_ranges[i-1][0])
            print()
            print("Guess with function equal to zero: " + str(guess_with_zero) + ".")
        elif abs(f[1]) <= tol:
            guess_with_zero.append(x_ranges[i-1][1])
            print()
            print("Guess with function equal to zero: " + str(guess_with_zero) + ".")
    
    if 'guess_with_zero' in locals():
        return guess_with_zero
    else:
        return function_ranges_with_zero


x_ranges = divide_into_ranges(rho_guess[0],rho_guess[1],n)
function_ranges = []
for x_min, x_max in x_ranges:
    function_ranges.append([function_P_rho(x_min), function_P_rho(x_max)])
function_zeros_range(x_ranges, function_ranges, tol)


