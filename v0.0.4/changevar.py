# -*- coding: utf-8 -*-
"""
Description
--
Find zero point using exclusively Bisection Method, exclusively Newton-Raphson Method, or with an hybrid method.

Methods based on version 0.0.2 of each one, but with some changes on callings to functions using T_target, xmol_target and P_target.
"""

"""
Version 0.0.4
--
Added description into each line for some general function and for Bissection Method

TO DO:
    - Introduce relative errors for all methods;
    - Description of each line in Newton-Raphson and Hybrid Methods;
    - Create a checker to determine if the function variables are within the specified range from IAPWS;
    

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


# --------------------------------------------------------------------------
# Use SI on each variable, as prof. Simões observed

# P_target = pressure in Pa
# T_target = temperature in K 
# xmass_target = ammonia mass fraction
# rho_guess_min = density in kg/m³
# rho_guess_max = density in kg/m³

# n = number to divide guess range
# tol = tolerance to zero
# N = limit for numbers of iterations


# --------------------------------------------------------------------------
# Defining function f(rho) to find a zero when pressure target is found. All the other variables remains constant when the methods to find zero were applied.
# --------------------------------------------------------------------------
def function_P_rho(rho, T_target, xmol_target, P_target): # Definition of function
    try:       # if the function fails, move to the other line "except"
        function = (iapws.ammonia.H2ONH3()._prop(rho, T_target, xmol_target).get("P") * ut.pressure_MPa_to_Pa) - P_target       # Function that gives zero when target pressure is acchieved
        return function     # return function
    except:     # When the function has an error, goes to this line
        return print("Error in guess = " + str(rho) + ". \n")   # Expression to inform that it has an error.


# --------------------------------------------------------------------------
# Divide a range into n ranges/divisions and create a list of them.
def divide_into_ranges(x_min, x_max, n_divisions):  # Definition of divisions to make
    diff = (x_max - x_min)/n_divisions      # size of each division
    list = []       # creates a list
    i = 1   # variable created just to iterate with each division
    while i <= n_divisions:     # while loop for each iteration/division
        x_min_new = x_min + diff*(i-1)      # minimum value 
        x_max_new = x_min_new + diff        # maximum value
        list.append([x_min_new, x_max_new])     # put the values into a list with two numbers
        i +=1   # iteration variable increases to the next step
    return list     # returns the list with the ranges divided


# --------------------------------------------------------------------------
# Check if the the list of ranges and list of results (from funtion) has a zero, and return a list with just the ranges within a zero, or even a zero point when those are within the list in input.
# --------------------------------------------------------------------------
def x_with_zero_ranges(x_ranges, function_ranges, tolerance):   # Define the function
    x_ranges_with_zero = []     # Creates the output with the x ranges with zero(s) inside.
    i = 0   # iteration number
    for f in function_ranges:   # for each function range (list) in the list specified (list inside a list)
        if (f[0]*f[1] < 0): # Check if the x_min and x_máx has different signals (to cross a zero)
            x_ranges_with_zero.append(x_ranges[i])  # Add this zero range
        elif abs(f[0]) <= tolerance:    # Check if x_min is a zero (within tolerance)
            x_ranges_with_zero.append(x_ranges[i][0])  # Add this zero range
        elif abs(f[1]) <= tolerance:    # Check if x_máx is a zero (within tolerance)
            x_ranges_with_zero.append(x_ranges[i][1])  # Add this zero range
        i += 1  # Increases the iteration variable to go to the next value
    return x_ranges_with_zero   # Returnonly the x ranges with zero(s) in the middle.


# --------------------------------------------------------------------------
# Defines the x necessary to reach the zero in the function above (i.e. function_P_rho), according to some parameters needed to introduced already.
# --------------------------------------------------------------------------

# EDITING TO REACH RELATIVE ERROR ANALYSIS
def find_zero_bisection_rho(x_min, x_max, T_target, xmol_target, P_target, tolerance, limit_iterations):    # Defines this function 
    i = 1   # variable just to iterate
    x_avg_before = 0
    while i <= limit_iterations:    # check if the number of iteration at the moment is above the limit number of iterations
        x_avg = x_min + (x_max - x_min)/2   # define the middle point (x average) between the x_min and x_max values
        function_x_min = function_P_rho(x_min, T_target, xmol_target, P_target) # [bar]     # Define the function solution at x_min
        function_x_avg = function_P_rho(x_avg, T_target, xmol_target, P_target) # [bar]     # Define the function solution at x_máx
        if i > 1 and (function_x_avg == 0 or abs((x_avg-x_avg_before)/x_avg_before) < tolerance):   # Check if function in x_avg is zero or error in x difference is less than tolerance
        # if function_x_avg == 0 or abs((x_max - x_min)/2) < tolerance:   # Check if function in x_avg is zero or error in x difference is less than tolerance
            return x_avg    # return x average
        elif function_x_min*function_x_avg > 0:     # Check if the function solution in x_avg has the same signal than in x_min
            x_min = x_avg   # Change x_min to x_avg
        else:
            x_max = x_avg   # Change x_máx to x_avg
        i += 1  # Increase the number of iteration variable
        x_avg_before = x_avg
        function_x_avg_before = function_x_avg
        
        
    if i > limit_iterations:   # If the number of iterations exceed the number specified by the user, then shows an error (just printed)
        return print("Method failed after " + str(i-1) + " iteration(s).")  # Error message shown to the user.


# --------------------------------------------------------------------------
# Defines the ammonia-water properties using pressure instead of density (rho), based on IAWPS library, requiring guesses of rho and other variables.
# --------------------------------------------------------------------------
def my_prop_bissection(P_target, T_target, xmass_target, rho_guess_min, rho_guess_max, n_divisions, tolerance, limit_iterations):  # Defines the function
    MW_NH3 = 17.03026 # [g/mol]     Referenced number from iapws, inside nh3h2o documentation.
    MW_H2O = 18.015268 # [g/mol]    Referenced number from iapws, inside nh3h2o documentation.
    xmol_target = ut.xmass_to_xmol(xmass_target,MW_NH3,MW_H2O) # [mol/mol]  # Transformation of mass fraction into molar fraction
    
    x_ranges = divide_into_ranges(rho_guess_min,rho_guess_max,n_divisions)  # Divide rho range into small ones
    function_ranges = []    # Define a list of function ranges
    for x_min, x_max in x_ranges:   # For each x_min and x_max in the list of list from x_ranges, do:
        function_ranges.append([function_P_rho(x_min, T_target, xmol_target, P_target), function_P_rho(x_max, T_target, xmol_target, P_target)])    # Add a list of function solution pair in x_min and x_max (range)
    x_zero_ranges = x_with_zero_ranges(x_ranges, function_ranges, tolerance)  # Check if there are zeros between the ranges

# Check if I have a range or zero points, and apply bisection method just within ranges
    x_with_zero = []    # Create a list of this variable
    for x in x_zero_ranges:     # Take each range into x_zero_ranges
        if not(isinstance(x, list)):    # If the "range" is not a list (i.e. is a number, for example), then is added directly into the answer,m because it is necessarily a zero.
            x_with_zero.append(x) # Added into the answer list
        else:   # If is not a zero directly
            x_zero = find_zero_bisection_rho(x[0],x[1],T_target,xmol_target,P_target,tolerance,limit_iterations)     # Apply the bissection method to find a zero
            x_with_zero.append(x_zero)  # Add the found zero into the list of x with solution equal to zero in function_P_rho
            
   
# Printing answers
    if x_with_zero == []:   # If the list is empty means it was not found a zero, so a message appears
        print("Those intervals does not have a zero point. Try refine the range or change values for guesses.")     # Message that a zero has not found.
    else:   # If it has something inside the list, then:
        rho_found = x_with_zero[0]   # IMPORTANT: Here i am assuming only the first zero point is valid. This can be changed if necessary.
        
    prop = {}   # Creating the properties dictionary
    prop["M"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("M")     # Property Mixture molecular mass, [g/mol]
    prop["rho"] = rho_found     # Property Density, [kg/m³]
    prop["u"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("u")     # Property Specific internal energy, [kJ/kg]
    prop["s"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("s")     # Property Specific entropy, [kJ/kgK]
    prop["h"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("h")     # Property Specific enthalpy, [kJ/kg]
    prop["g"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("g")     # Property Specific Helmholtz energy, [kJ/kg]
    prop["a"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("a")     # Property Specific gibbs energy, [kJ/kg]
    prop["cv"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("cv")     # Property Specific isochoric heat capacity, [kJ/kgK]
    prop["cp"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("cp")     # Property Specific isobaric heat capacity, [kJ/kgK]
    prop["w"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("w")     # Property Speed of sound, [m/s]
    prop["fugH2O"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("fugH2O")     # Property Fugacity of water, [-]
    prop["fugNH3"] = iapws.ammonia.H2ONH3()._prop(rho_found, T_target, xmol_target).get("fugNH3")     # Property Fugacity of ammonia, [-]
    
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
    
    
    
    
    
    
    