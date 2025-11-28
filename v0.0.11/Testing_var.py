#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 16:43:01 2021

@author: unknown
"""

"""
Version 0.0.1 to 0.0.6
Nothing done.

Version 0.0.7
Try triple point automatically

Version 0.0.8
-- nothing done --
"""


# -----------------------------------------------------------------------------

# import iapws.ammonia as ammonia
# import changevar
# import useful_tools as tools

# # MW_NH3 = 17.03026 # [g/mol]     Referenced number from iapws, inside nh3h2o documentation.
# # MW_H2O = 18.015268 # [g/mol]    Referenced number from iapws, inside nh3h2o documentation.

# filewrite = open("test_Ttr.csv", 'w')
# filewrite.write("x [mol ammonia/mol];T [K]\n")
# filewrite.close()

# x = 0
# fileappen = open("test_Ttr.csv", 'a')


# while x >= 0 and x <= 1:
#     T = ammonia.Ttr(x)
#     # P = changevar.my_prop_hybrid(P, T, tools.xmol_to_xmass(x, MW_NH3, MW_H2O), 0.0001, 150.0001, 100, pow(10,-1), pow(10,-6), 3000)
#     fileappen.write(str(x) + ";" + str(T) + "\n")
#     x+=0.001
# fileappen.close()

# fileread = open("test_Ttr.csv", 'r')
# print(fileread.read())
# filewrite.close()

# -----------------------------------------------------------------------------

# import changevar

# P = 10000
# T = 100
# xmass = 0 

# rho_min = 0.0001
# rho_max = 150.0001
# init_divisions = 100
# tolerance = pow(10,-4)
# limit_iterations = 3000

# delta_P = 10000
# delta_T = 100
# delta_xmass = 0 


# filewrite = open("test_changevar.txt", 'w')
# filewrite.write("P [Pa];T [K];xmass [kg ammonia/kg]\n")
# filewrite.close()

# fileappen = open("test_changevar.txt", 'a')
# while P >= 10000 and P <= 40*pow(10,6):
#     try:
#         Properties = changevar.my_prop_bissection(P, T, xmass, rho_min, rho_max, init_divisions, tolerance, limit_iterations)
#         fileappen.write(str(Properties) + "\n\n\n")
#         P += delta_P
#     except:
#         fileappen.write("-\n\n\n")
#         P += delta_P
# fileappen.close()

# fileread = open("test_changevar.txt", 'r')
# print(fileread.read())
# filewrite.close()

# -----------------------------------------------------------------------------

import changevar
import time

start_time_BS = time.time() # Measuring above 0.1 seconds

print(changevar.my_prop_bissection(1500000, 393.25, 0.42, 0.1, 25.1, 10, pow(10,-7), 30))

print()
executionTime = (time.time() - start_time_BS)
print("--- %s seconds ---" % executionTime)



start_time_NR = time.time() # Measuring above 0.1 seconds

print()
print(changevar.my_prop_newraph(1500000, 393.25, 0.42, 0.1, 25.1, 10, pow(10,-7), 30))

print()
executionTime = (time.time() - start_time_NR)
print("--- %s seconds ---" % executionTime)



start_time_HB = time.time() # Measuring above 0.1 seconds

print()
print(changevar.my_prop_hybrid(1500000, 393.25, 0.42, 0.1, 25.1, 10, pow(10,-1), pow(10,-7), 30))

print()
executionTime = (time.time() - start_time_HB)
print("--- %s seconds ---" % executionTime)








