# -*- coding: utf-8 -*-
"""

Find zero point using exclusively Bisection Method, exclusively Newton-Raphson Method, or with an hybrid method.

It was based on version 0.0.2 of each method, but with some changes on callings to functions using T_target, xmol_target and P_target.

"""


import iapws
import useful_tools as ut

# import time

   
# start_time = time.time() # Measuring above 0.1 seconds
"""
P_target = 15.0*ut.pressure_bar_to_Pa # [Pa]
T_target = 273.15+120 # [K]
xmass_target = 0.42 # [kg NH3/kg H2O]

rho_guess_min = 0.1 # [kg/m³]
rho_guess_max = 25.1 # [kg/m³]

n = 10 # number to divide guess range
tol = pow(10,-7) # Tolerance to zero
N = 3000 # limit for numbers of itera
"""


def function_P_rho(rho, T_target, xmol_target, P_target):
    try:
        function = (iapws.ammonia.H2ONH3()._prop(rho, T_target, xmol_target).get("P") * ut.pressure_MPa_to_Pa) - P_target
        return function
    except:
        return print("Error in guess = " + str(rho) + ". \n")


def find_zero_bisection_rho(guess_min, guess_max, T_target, xmol_target, P_target, tolerance, limit_iterations):
    i = 1
    while i <= limit_iterations:
        guess_avg = guess_min + (guess_max - guess_min)/2
        function_guess_min = function_P_rho(guess_min, T_target, xmol_target, P_target) # [bar]
        function_guess_avg = function_P_rho(guess_avg, T_target, xmol_target, P_target) # [bar]
        if function_guess_avg == 0 or abs((guess_min - guess_max)/2) < tolerance:
            return guess_avg
        elif function_guess_min*function_guess_avg > 0:
            guess_min = guess_avg
        else:
            guess_max = guess_avg
        i += 1
        
    if i > N:
        return print("Method failed after " + str(i-1) + " iteration(s).")


# -----------------------------------------------------
# Use SI on each variable, as prof. Simões observed

# P_target = pressure in Pa
# T_target = temperature in K 
# xmass_target = ammonia mass fraction
# rho_guess_min = density in kg/m³
# rho_guess_max = density in kg/m³

# n = number to divide guess range
# tol = tolerance to zero
# N = limit for numbers of iterations

def my_prop_bissection(P_target, T_target, xmass_target, rho_guess_min, rho_guess_max, n, tol, N):
    MW_NH3 = 17.03026 # [g/mol] referência do iapws, no docuemnto de nh3h2o
    MW_H2O = 18.015268 # [g/mol] referência do iapws, no docuemnto de nh3h2o
    xmol_target = ut.xmass_to_xmol(xmass_target,MW_NH3,MW_H2O) # [mol/mol]
    
    x_ranges = ut.divide_into_ranges(rho_guess_min,rho_guess_max,n)
    # print(x_ranges)
    function_ranges = []
    for x_min, x_max in x_ranges:
        function_ranges.append([function_P_rho(x_min, T_target, xmol_target, P_target), function_P_rho(x_max, T_target, xmol_target, P_target)])
    # print(function_ranges)
    x_zero_ranges = ut.x_with_zero_ranges(x_ranges, function_ranges, tol)
    # print(ut.x_with_zero_ranges(x_ranges, function_ranges, tol))

    # Check if I have a range or zero points, and apply bisection method just within ranges
    x_with_zero = []
    for x in x_zero_ranges:
        if not(isinstance(x, list)):
            x_with_zero.append(x)
        else:
            x_zero = find_zero_bisection_rho(x[0],x[1],T_target,xmol_target,P_target,tol,N)
            x_with_zero.append(x_zero)   

    # prop["rho"] = x_with_zero[0]
    
    # Print answers
    if x_with_zero == []:
        print("Those intervals does not have a zero point. Try refine the range or change values for guesses.")
        # print("Function ranges: "+ str(function_ranges))
    else:
        rho_found = x_with_zero[0]   # Here i am assuming only the first zero point is valid.
    prop = {}
    prop["M"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("M")
    prop["rho"] = rho_found
    prop["u"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("u")
    prop["s"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("s")
    prop["h"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("h")
    prop["g"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("g")
    prop["a"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("a")
    prop["cv"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("cv")
    prop["cp"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("cp")
    prop["w"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("w")
    prop["fugH2O"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("fugH2O")
    prop["fugNH3"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("fugNH3")
    
    return prop

    
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
    

def find_zero_newtonraphson_rho(guess, T_target, xmol_target, P_target, tolerance, limit_iterations):
    i = 1
    while i <= limit_iterations:
        guess_minusStep = guess*(1-tolerance) # tolerance used to estimate differential curve
        guess_plusStep = guess*(1+tolerance)
        if guess_minusStep <= 0 or guess_plusStep <= 0:
            return print("Error! Guess diverged after " + str(i) + " iterations. Check the input conditions.")
        
        differential_function = (function_P_rho(guess_plusStep, T_target, xmol_target, P_target)-function_P_rho(guess_minusStep, T_target, xmol_target, P_target))/(guess_plusStep-guess_minusStep)
        guess_now = guess-function_P_rho(guess, T_target, xmol_target, P_target)/differential_function
        if abs(guess_now-guess) < tolerance:
            return guess_now #print("guess = "+ str(guess_now)+", after "+str(i)+" iterations.")
        else:
            guess = guess_now
            i += 1
            
    if i > limit_iterations:
        return print("Method failed after " + str(i-1) + " iteration(s).")
    
    
def my_prop_newraph(P_target, T_target, xmass_target, rho_guess_min, rho_guess_max, n, tol, N):
    MW_NH3 = 17.03026 # [g/mol] referência do iapws, no docuemnto de nh3h2o
    MW_H2O = 18.015268 # [g/mol] referência do iapws, no docuemnto de nh3h2o
    xmol_target = ut.xmass_to_xmol(xmass_target,MW_NH3,MW_H2O) # [mol/mol]
    
    x_ranges = ut.divide_into_ranges(rho_guess_min,rho_guess_max,n)
    function_ranges = []
    for x_min, x_max in x_ranges:
        function_ranges.append([function_P_rho(x_min, T_target, xmol_target, P_target), function_P_rho(x_max, T_target, xmol_target, P_target)])
    x_zero_ranges = ut.x_with_zero_ranges(x_ranges, function_ranges, tol)

    x_guess = listoflists_into_avgnumbers(x_zero_ranges)
    # print(x_guess)
    x_with_zero = []
    for x in x_guess:
        x_zero = find_zero_newtonraphson_rho(x, T_target, xmol_target, P_target, tol, N)
        x_with_zero.append(x_zero)


    # Print answers
    if x_with_zero == []:
        print("Those intervals does not have a zero point or diverged. Try refine the range or change values for guesses.")
        # print("Function ranges: "+ str(function_ranges))
    else:
        rho_found = x_with_zero[0]   # Here i am assuming only the first zero point is valid.

    prop = {}
    prop["M"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("M")
    prop["rho"] = rho_found
    prop["u"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("u")
    prop["s"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("s")
    prop["h"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("h")
    prop["g"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("g")
    prop["a"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("a")
    prop["cv"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("cv")
    prop["cp"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("cp")
    prop["w"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("w")
    prop["fugH2O"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("fugH2O")
    prop["fugNH3"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("fugNH3")
    
    return prop
    


def my_prop_hybrid(P_target, T_target, xmass_target, rho_guess_min, rho_guess_max, n, tol_BS, tol_NR, N):
    MW_NH3 = 17.03026 # [g/mol] referência do iapws, no docuemnto de nh3h2o
    MW_H2O = 18.015268 # [g/mol] referência do iapws, no docuemnto de nh3h2o
    xmol_target = ut.xmass_to_xmol(xmass_target,MW_NH3,MW_H2O) # [mol/mol]
    
    x_ranges = ut.divide_into_ranges(rho_guess_min,rho_guess_max,n)
    # print(x_ranges)
    function_ranges = []
    for x_min, x_max in x_ranges:
        function_ranges.append([function_P_rho(x_min, T_target, xmol_target, P_target), function_P_rho(x_max, T_target, xmol_target, P_target)])
    # print(function_ranges)
    x_zero_ranges = ut.x_with_zero_ranges(x_ranges, function_ranges, tol_NR)
    # print(ut.x_with_zero_ranges(x_ranges, function_ranges, tol))


    x_with_zero_BS = []
    for x in x_zero_ranges:
        if not(isinstance(x, list)):
            x_with_zero_BS.append(x)
        else:
            x_zero_BS = find_zero_bisection_rho(x[0],x[1],T_target,xmol_target,P_target,tol_BS,N)
            x_with_zero_BS.append(x_zero_BS)  


    x_with_zero_NR = []
    for x in x_with_zero_BS:
        x_zero_NR = find_zero_newtonraphson_rho(x, T_target, xmol_target, P_target, tol_NR, N)
        x_with_zero_NR.append(x_zero_NR)

   
    # Print answers
    if x_with_zero_NR == []:
        print("Those intervals does not have a zero point or diverged. Try refine the range or change values for guesses.")
        # print("Function ranges: "+ str(function_ranges))
    else:
        rho_found = x_with_zero_NR[0]   # Here i am assuming only the first zero point is valid.

    prop = {}
    prop["M"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("M")
    prop["rho"] = rho_found
    prop["u"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("u")
    prop["s"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("s")
    prop["h"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("h")
    prop["g"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("g")
    prop["a"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("a")
    prop["cv"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("cv")
    prop["cp"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("cp")
    prop["w"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("w")
    prop["fugH2O"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("fugH2O")
    prop["fugNH3"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("fugNH3")
    
    return prop
    
    
    
    
    
    
    