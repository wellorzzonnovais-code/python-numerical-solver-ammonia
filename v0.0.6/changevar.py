# -*- coding: utf-8 -*-
"""
Description
--
Find zero point using exclusively Bisection Method, exclusively Newton-Raphson Method, or with an hybrid method.

Methods based on version 0.0.2 of each one, but with some changes on callings to functions using T_target, xmol_target and P_target.
"""

"""
Log

--------------------------------------
Version 0.0.4
--
Added description into each line for some general function and for Bissection Method

TO DO:
    - Introduce relative errors for all methods;
    - Description of each line in Newton-Raphson and Hybrid Methods;
    - Create a checker to determine if the function variables are within the specified range from IAPWS;

--------------------------------------
Version 0.0.5
--
Added description into each line for Newton-Raphson and Hybrid Method
Introduced relative errors (for x values) in Bissection methods;
Updated main variables name and description in header only for general functions and Bissection Method.

TO DO:
    - Update main variables name and description in header
    - In "x_with_zero_ranges", I must define relative error based on x values, not fucntion solutions. I can delete the lines that verify if the values are within tolerance, because it is a relative error based on last iteration, which is not applied inside this function "x_with_zero_ranges".
    - Create a checker to determine if the function variables are within the specified range from IAPWS;
    - Create an automated error test within these ranges;
    - Introduce error description for each function/definition.

--------------------------------------
Version 0.0.6
--
In function "x_with_zero_ranges", it was eliminated checking zeros with function solution and tolerance: it makes no sense to test here;
Reorganized some variables in functions for Newton-Raphson and Hybrid methods;
Introduced relative errors (for function to be zero - function_P_rho - solutions) in Newton-Raphson method;
Updated relative errors function to be based on function solutoins instead of density values;


TO DO:
    - Create a checker to determine if the function variables are within the specified range from IAPWS;
    - Introduce error description for each function/definition.
    - Create an automated error test within these ranges;

"""

import iapws
import useful_tools as ut


# --------------------------------------------------------------------------
# Using SI on each variable, as prof. Simões observed
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# Defining function f(rho) to find a zero when pressure target is found. All the other variables remains constant when the methods to find zero were applied.
# --------------------------------------------------------------------------
# rho = density [kg/m³]
# T = temperature [K]
# xmol = molar fraction [mol ammonia/mol solution]
# P_target = pressure target [Pa]

def function_P_rho(rho, T, xmol, P_target): # Definition of function
    try:       # if the function fails, move to the other line "except"
        function = (iapws.ammonia.H2ONH3()._prop(rho, T, xmol).get("P") * ut.pressure_MPa_to_Pa) - P_target       # Function that gives zero when target pressure is acchieved
        return function     # return function
    except:     # When the function has an error, goes to this line
        return print("Error in guess = " + str(rho) + ". \n")   # Expression to inform that it has an error.


# --------------------------------------------------------------------------
# Divide a range into n ranges/divisions and create a list of them.
# --------------------------------------------------------------------------
# x_min = lower value of defined range
# x_max = higher value of defined range
# n_divisions = number of divisions

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
# Check if the list of ranges and list of solutions from these ranges (from funtion) has a zero, and return a list with just the ranges with a zero.
# --------------------------------------------------------------------------
# x_ranges = list containing a list of lower and higher values of a range.
# function_ranges = list containing a list of solutions for lower and higher values of a range, respectively.

def x_with_zero_ranges(x_ranges, function_ranges):   # Define the function
    x_ranges_with_zero = []     # Creates the output with the x ranges with zero(s) inside.
    i = 0   # iteration number
    for f in function_ranges:   # for each function range (list) in the list specified (list inside a list)
        if (f[0]*f[1] < 0): # Check if the x_min and x_máx has different signals (to cross a zero)
            x_ranges_with_zero.append(x_ranges[i])  # Add this zero range
        i += 1  # Increases the iteration variable to go to the next value
    return x_ranges_with_zero   # Returnonly the x ranges with zero(s) in the middle.


# --------------------------------------------------------------------------
# Defines the x necessary to reach the zero in the function above (i.e. function_P_rho), according to some parameters needed to introduced already, using the Bissection zero finding method.
# --------------------------------------------------------------------------
# x_min = lower value of defined x range (for density range, [kg/m³])
# x_max = higher value of defined x range (for density range, [kg/m³])
# T = temperature [K]
# xmol = molar fraction [mol ammonia/mol solution]
# P_target = pressure target [Pa]
# tolerance = specified tolerance for x to reach zero solution [-]
# limit_iterations = limit number of iterations [-]

def find_zero_bisection_rho(P_target, T, xmol, x_min, x_max, tolerance, limit_iterations):    # Defines this function 
    i = 1   # variable just to iterate
    # x_avg_before = 0
    # function_x_avg_before = 0
    while i <= limit_iterations:    # check if the number of iteration at the moment is above the limit number of iterations
        x_avg = x_min + (x_max - x_min)/2   # define the middle point (x average) between the x_min and x_max values
        function_x_min = function_P_rho(x_min, T, xmol, P_target) # [bar]     # Define the function solution at x_min
        function_x_avg = function_P_rho(x_avg, T, xmol, P_target) # [bar]     # Define the function solution at x_máx
        if i > 1 and (abs(function_x_avg/P_target) < tolerance):   # Check if function in x_avg is zero or relative error in function solution is less than tolerance
            return x_avg    # return x average
        elif function_x_min*function_x_avg > 0:     # Check if the function solution in x_avg has the same signal than in x_min
            x_min = x_avg   # Change x_min to x_avg
        else:
            x_max = x_avg   # Change x_máx to x_avg
        i += 1  # Increase the number of iteration variable
        # x_avg_before = x_avg
        # function_x_avg_before = function_x_avg
        
        
    if i > limit_iterations:   # If the number of iterations exceed the number specified by the user, then shows an error (just printed)
        return print("Method failed after " + str(i-1) + " iteration(s).")  # Error message shown to the user.


# --------------------------------------------------------------------------
# Defines the ammonia-water properties using pressure instead of density (rho), based on IAWPS library, requiring guesses of rho and other variables to apply Bissection Method to find zeros.
# --------------------------------------------------------------------------
# P_target = pressure target [Pa]
# T = temperature [K]
# xmass = mass fraction [kg ammonia/kg solution]
# density_guess_min = lower value of guessed solution density range [kg/m³]
# density_guess_max = higher value of guessed solution density range [kg/m³]
# n_divisions = number of divisions to divide initial rho range [-]
# tolerance = relative error specified for function to reach zero solution [-]
# limit_iterations = limit number of iterations [-]

def my_prop_bissection(P_target, T, xmass, density_guess_min, density_guess_max, n_divisions, tolerance, limit_iterations):  # Defines the function
    MW_NH3 = 17.03026 # [g/mol]     Referenced number from iapws, inside nh3h2o documentation.
    MW_H2O = 18.015268 # [g/mol]    Referenced number from iapws, inside nh3h2o documentation.
    xmol = ut.xmass_to_xmol(xmass,MW_NH3,MW_H2O) # [mol/mol]  # Transformation of mass fraction into molar fraction

    x_ranges = divide_into_ranges(density_guess_min,density_guess_max,n_divisions)  # Divide rho range into small ones
    function_ranges = []    # Define a list of function ranges
    for x_min, x_max in x_ranges:   # For each x_min and x_max in the list of list from x_ranges, do:
        function_ranges.append([function_P_rho(x_min, T, xmol, P_target), function_P_rho(x_max, T, xmol, P_target)])    # Add a list of function solution pair in x_min and x_max (range)
    x_zero_ranges = x_with_zero_ranges(x_ranges, function_ranges)  # Check if there are zeros between the ranges

# Check if I have a range or zero points, and apply bisection method just within ranges
    x_with_zero = []    # Create a list of this variable
    for x in x_zero_ranges:     # Take each range into x_zero_ranges
        if not(isinstance(x, list)):    # If the "range" is not a list (i.e. is a number, for example), then is added directly into the answer,m because it is necessarily a zero.
            x_with_zero.append(x) # Added into the answer list
        else:   # If is not a zero directly
            x_zero = find_zero_bisection_rho(P_target,T,xmol,x[0],x[1],tolerance,limit_iterations)     # Apply the bissection method to find a zero           
            x_with_zero.append(x_zero)  # Add the found zero into the list of x with solution equal to zero in function_P_rho
            
   
# Printing answers
    if x_with_zero == []:   # If the list is empty means it was not found a zero, so a message appears
        print("Those intervals does not have a zero point. Try refine the range or change values for guesses.")     # Message that a zero has not found.
    else:   # If it has something inside the list, then:
        rho_found = x_with_zero[0]   # IMPORTANT: Here i am assuming only the first zero point is valid. This can be changed if necessary.
        
    prop = {}   # Creating the properties dictionary
    prop["M"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("M")     # Property Mixture molecular mass, [g/mol]
    prop["rho"] = rho_found     # Property Density, [kg/m³]
    prop["u"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("u")     # Property Specific internal energy, [kJ/kg]
    prop["s"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("s")     # Property Specific entropy, [kJ/kgK]
    prop["h"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("h")     # Property Specific enthalpy, [kJ/kg]
    prop["g"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("g")     # Property Specific Helmholtz energy, [kJ/kg]
    prop["a"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("a")     # Property Specific gibbs energy, [kJ/kg]
    prop["cv"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("cv")     # Property Specific isochoric heat capacity, [kJ/kgK]
    prop["cp"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("cp")     # Property Specific isobaric heat capacity, [kJ/kgK]
    prop["w"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("w")     # Property Speed of sound, [m/s]
    prop["fugH2O"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("fugH2O")     # Property Fugacity of water, [-]
    prop["fugNH3"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("fugNH3")     # Property Fugacity of ammonia, [-]
    
    return prop

# --------------------------------------------------------------------------
# Defines the transformation from range lists (list of "lists with two numbers") into list of average numbers in the middle of the range.
# --------------------------------------------------------------------------
# list_of_guesses_or_ranges = a list containing ranges (two numbers inside a list) or even values.

def listoflists_into_avgnumbers(list_of_guesses_or_ranges):     # Define the transformation function
    list_avg = []   # Creates a list named list_avg
    for n in list_of_guesses_or_ranges:     # For each range specified in the list, do:
        if not(isinstance(n, list)):    # If n is not a list (which means it is not a range, for this case)
            list_avg.append(n)  # Include already in the average list
        else:   # If not:
            list_avg.append((n[0]+n[1])/2)  # Includes in the list the average number of range
    return list_avg     # Return the list with average values
    

# --------------------------------------------------------------------------
# Defines the x necessary to reach the zero in the function above (i.e. function_P_rho), according to some parameters needed to introduced already, using the Newton-Raphson zero finding method.
# --------------------------------------------------------------------------
# P_target = pressure target [Pa]
# T = temperature [K]
# xmol = molar fraction [mol ammonia/mol solution]
# x_value = value of guessed solution for function_P_rho (for solution density, [kg/m³])
# n_divisions = number of divisions to divide initial x range [-]
# tolerance = relative error specified for function to reach zero solution [-]
# limit_iterations = limit number of iterations [-]

def find_zero_newtonraphson_rho(P_target, T, xmol, x_value, tolerance, limit_iterations):     # Define this functions with those parameters
    i = 1   # Define the iteration variable
    while i <= limit_iterations:    # Definition of loop according to number of iterations
        x_value_minusStep = x_value*(1-tolerance)   # Calculate a nearby x (lower value) to estimate the differential funciton numerically with "tolerance" variation in x. NOTE: other values than tolerance can be used.
        x_value_plusStep = x_value*(1+tolerance)   # Calculate a nearby x (lower value) to estimate the differential funciton numerically with "tolerance" variation in x. NOTE: other values than tolerance can be used.
        if x_value_minusStep <= 0 or x_value_plusStep <= 0:     # When these values are negative, FOR DENSITY, it has diverged.
            return print("Error! Guess diverged after " + str(i) + " iterations, and a negative number for x was found. Check the input conditions.")
        differential_function = (function_P_rho(x_value_plusStep, T, xmol, P_target)-function_P_rho(x_value_minusStep, T, xmol, P_target))/(x_value_plusStep-x_value_minusStep)     # According to Newton-Raphson Method the numerical differential was created.
        if i==1:      
            function_x_value = function_P_rho(x_value, T, xmol, P_target)   # Calculates the solution for x_value
        x_value_new = x_value-function_x_value/differential_function    # According to Newton-Raphson Method, a new x_value was created.
        function_x_value_new = function_P_rho(x_value_new, T, xmol, P_target)   # Calculates the solution for x_value_new
        if abs(function_x_value_new/P_target) < tolerance:    # Check if the specified number difference reached the relative tolerance.
            return x_value_new  # Return x_value value found
        else:   # Otherwise:
            x_value = x_value_new   # Update x_value to be the new one found in last iteration loop.
            function_x_value = function_x_value_new   # Update function_x_value to be the new one found in last iteration loop.
            i += 1  # Update the iteration tracker variable
 
    if i > limit_iterations:    # If the number of iteration has exceeded, an error message appears
        return print("Method failed after " + str(i-1) + " iteration(s).")  # Error message when limit_iterations has acchieved.
    


# --------------------------------------------------------------------------
# Defines the ammonia-water properties using pressure instead of density (rho), based on IAWPS library, requiring guesses of rho and other variables to apply Newton-Raphson Method to find zeros.
# --------------------------------------------------------------------------
# P_target = pressure target [Pa]
# T = temperature [K]
# xmass = mass fraction [kg ammonia/kg solution]
# density_guess_min = lower value of guessed solution density range [kg/m³]
# density_guess_max = higher value of guessed solution density range [kg/m³]
# n_divisions = number of divisions to divide initial rho range [-]
# tolerance = relative error specified for function to reach zero solution [-]
# limit_iterations = limit number of iterations [-]

def my_prop_newraph(P_target, T, xmass, density_guess_min, density_guess_max, n_divisions, tolerance, limit_iterations):  # Defines the function
    MW_NH3 = 17.03026 # [g/mol]     Referenced number from iapws, inside nh3h2o documentation.
    MW_H2O = 18.015268 # [g/mol]    Referenced number from iapws, inside nh3h2o documentation.
    xmol = ut.xmass_to_xmol(xmass,MW_NH3,MW_H2O) # [mol/mol]  # Transformation of mass fraction into molar fraction
    
    x_ranges = divide_into_ranges(density_guess_min,density_guess_max,n_divisions)  # Divide rho range into small ones
    function_ranges = []    # Define a list of function ranges
    for x_min, x_max in x_ranges:   # For each x_min and x_max in the list of list from x_ranges, do:
        function_ranges.append([function_P_rho(x_min, T, xmol, P_target), function_P_rho(x_max, T, xmol, P_target)])    # Add a list of function solution pair in x_min and x_max (range)
    x_zero_ranges = x_with_zero_ranges(x_ranges, function_ranges)  # Check if there are zeros between the ranges


# Convert ranges into a number (average of range) and apply Newton-Raphson method just with these values.
    x_guess = listoflists_into_avgnumbers(x_zero_ranges)    # Converting x with zero solution ranges with function_P_rho
    x_with_zero = []    # Create a variable to store x with zero solutions
    for x in x_guess:   # For each average zero number in x_guess, do:
        x_zero = find_zero_newtonraphson_rho(P_target, T, xmol, x, tolerance, limit_iterations)    # Find a zero using this number as a guess with Newton-Raphson method
        x_with_zero.append(x_zero)  # Include the zero point into the list with x resulting in function equal to zero.


# Printing answers
    if x_with_zero == []:    # If the list is empty means it was not found a zero, so a message appears
        print("Those intervals does not have a zero point or diverged. Try refine the range or change values for guesses.")     # Message that a zero has not found.
    else:   # If it has something inside the list, then:
        rho_found = x_with_zero[0]   # IMPORTANT: Here i am assuming only the first zero point is valid. This can be changed if necessary.

    prop = {}   # Creating the properties dictionary
    prop["M"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("M")     # Property Mixture molecular mass, [g/mol]
    prop["rho"] = rho_found     # Property Density, [kg/m³]
    prop["u"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("u")     # Property Specific internal energy, [kJ/kg]
    prop["s"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("s")     # Property Specific entropy, [kJ/kgK]
    prop["h"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("h")     # Property Specific enthalpy, [kJ/kg]
    prop["g"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("g")     # Property Specific Helmholtz energy, [kJ/kg]
    prop["a"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("a")     # Property Specific gibbs energy, [kJ/kg]
    prop["cv"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("cv")     # Property Specific isochoric heat capacity, [kJ/kgK]
    prop["cp"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("cp")     # Property Specific isobaric heat capacity, [kJ/kgK]
    prop["w"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("w")     # Property Speed of sound, [m/s]
    prop["fugH2O"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("fugH2O")     # Property Fugacity of water, [-]
    prop["fugNH3"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("fugNH3")     # Property Fugacity of ammonia, [-]
    
    return prop
   
    

# --------------------------------------------------------------------------
# Defines the ammonia-water properties using pressure instead of density (rho), based on IAWPS library, requiring guesses of rho and other variables to apply Bissection (at first) and Newton-Raphson (second) Methods to find zeros. This is called Hybrid Method by the authors.
# --------------------------------------------------------------------------
# P_target = pressure target [Pa]
# T = temperature [K]
# xmass = mass fraction [kg ammonia/kg solution]
# density_guess_min = lower value of guessed solution density range [kg/m³]
# density_guess_max = higher value of guessed solution density range [kg/m³]
# n_divisions = number of divisions to divide initial rho range [-]
# tolerance = relative error specified for function to reach zero solution [-]
# limit_iterations = limit number of iterations [-]

def my_prop_hybrid(P_target, T, xmass, density_guess_min, density_guess_max, n_divisions, tolerance_BS, tolerance_NR, limit_iterations):
    MW_NH3 = 17.03026 # [g/mol]     Referenced number from iapws, inside nh3h2o documentation.
    MW_H2O = 18.015268 # [g/mol]    Referenced number from iapws, inside nh3h2o documentation.
    xmol = ut.xmass_to_xmol(xmass,MW_NH3,MW_H2O) # [mol/mol]  # Transformation of mass fraction into molar fraction


    x_ranges = divide_into_ranges(density_guess_min,density_guess_max,n_divisions)  # Divide rho range into small ones
    function_ranges = []    # Define a list of function ranges
    for x_min, x_max in x_ranges:   # For each x_min and x_max in the list of list from x_ranges, do:
        function_ranges.append([function_P_rho(x_min, T, xmol, P_target), function_P_rho(x_max, T, xmol, P_target)])    # Add a list of function solution pair in x_min and x_max (range)
    x_zero_ranges = x_with_zero_ranges(x_ranges, function_ranges)  # Check if there are zeros between the ranges


# Check if I have a range or zero points and apply bisection method just within ranges
    x_with_zero_BS = []    # Create a list of this variable for bissection method
    for x in x_zero_ranges:     # Take each range into x_zero_ranges
        if not(isinstance(x, list)):    # If the "range" is not a list (i.e. is a number, for example), then is added directly into the answer,m because it is necessarily a zero.
            x_with_zero_BS.append(x) # Added into the answer list
        else:   # If is not a zero directly
            x_zero_BS = find_zero_bisection_rho(P_target,T,xmol,x[0],x[1],tolerance_BS,limit_iterations)     # Apply the bissection method to find a zero        
            x_with_zero_BS.append(x_zero_BS)  # Add the found zero into the list of x with solution equal to zero in function_P_rho


# Apply Newton-Raphson method just within the values found with bissection method.
    x_with_zero_NR = []    # Create a variable to store x with zero solutions
    for x in x_with_zero_BS:   # For each average zero number in x_guess, do:
        x_zero_NR = find_zero_newtonraphson_rho(P_target, T, xmol, x, tolerance_NR, limit_iterations)    # Find a zero using this number as a guess with Newton-Raphson method
        x_with_zero_NR.append(x_zero_NR)  # Include the zero point into the list with x resulting in function equal to zero.


# Printing answers
    if x_with_zero_NR == []:    # If the list is empty means it was not found a zero, so a message appears
        print("Those intervals does not have a zero point or diverged. Try refine the range or change values for guesses.")     # Message that a zero has not found.
    else:   # If it has something inside the list, then:
        rho_found = x_with_zero_NR[0]   # IMPORTANT: Here i am assuming only the first zero point is valid. This can be changed if necessary.

    prop = {}   # Creating the properties dictionary
    prop["M"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("M")     # Property Mixture molecular mass, [g/mol]
    prop["rho"] = rho_found     # Property Density, [kg/m³]
    prop["u"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("u")     # Property Specific internal energy, [kJ/kg]
    prop["s"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("s")     # Property Specific entropy, [kJ/kgK]
    prop["h"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("h")     # Property Specific enthalpy, [kJ/kg]
    prop["g"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("g")     # Property Specific Helmholtz energy, [kJ/kg]
    prop["a"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("a")     # Property Specific gibbs energy, [kJ/kg]
    prop["cv"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("cv")     # Property Specific isochoric heat capacity, [kJ/kgK]
    prop["cp"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("cp")     # Property Specific isobaric heat capacity, [kJ/kgK]
    prop["w"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("w")     # Property Speed of sound, [m/s]
    prop["fugH2O"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("fugH2O")     # Property Fugacity of water, [-]
    prop["fugNH3"] = iapws.ammonia.H2ONH3()._prop(rho_found, T, xmol).get("fugNH3")     # Property Fugacity of ammonia, [-]
    
    return prop
    
    
    
    
    
    
    