#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 16:43:01 2021

@author: unknown
"""

import changevar
import time


start_time_BS = time.time() # Measuring above 0.1 seconds

print(changevar.my_prop_bissection(1500000, 393.25, 0.42, 0.1, 25.1, 10, pow(10,-7), 3000))

print()
executionTime = (time.time() - start_time_BS)
print("--- %s seconds ---" % executionTime)




start_time_NR = time.time() # Measuring above 0.1 seconds

print()
print(changevar.my_prop_newraph(1500000, 393.25, 0.42, 0.1, 25.1, 10, pow(10,-7), 3000))

print()
executionTime = (time.time() - start_time_NR)
print("--- %s seconds ---" % executionTime)




start_time_HB = time.time() # Measuring above 0.1 seconds

print()
print(changevar.my_prop_hybrid(1500000, 393.25, 0.42, 0.1, 25.1, 10, pow(10,-1), pow(10,-7), 3000))

print()
executionTime = (time.time() - start_time_HB)
print("--- %s seconds ---" % executionTime)

