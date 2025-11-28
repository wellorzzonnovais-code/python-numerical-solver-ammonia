# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# Tools for conversions

def xmass_to_xmol(xmass,mass_weight1,mass_weight2):
    xmol_converted = (xmass/mass_weight1)/(xmass/mass_weight1 + (1-xmass)/mass_weight2)
    return xmol_converted

pressure_MPa_to_Pa = pow(10,6)

pressure_bar_to_Pa = pow(10,5)

pressure_atm_to_Pa = 101325


# Tools to work with searching zeros in functions

# divide a range into n ranges and make a list of it them.
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

# elaborate a list of guesses range with a zero function point, or even a zero point when those are within the list in input.
def x_with_zero_ranges(x_ranges, function_ranges, tolerance): 
    # function_ranges_with_zero = []
    x_ranges_with_zero = []
    i = 0
    for f in function_ranges:
        i += 1
        if (f[0]*f[1] < 0):
            # function_ranges_with_zero.append(f)
            x_ranges_with_zero.append(x_ranges[i-1])
            # print()
            # print("There is a zero in " + str(x_ranges[i-1]))
            # print("with function between " + str(f) + ".")
        elif abs(f[0]) <= tolerance:
            x_ranges_with_zero.append(x_ranges[i-1][0])
            # print()
            # print("Guess with function equal to zero: " + str(guess_with_zero) + ".")
        elif abs(f[1]) <= tolerance:
            x_ranges_with_zero.append(x_ranges[i-1][1])
            # print()
            # print("Guess with function equal to zero: " + str(guess_with_zero) + ".")
    
    return x_ranges_with_zero
