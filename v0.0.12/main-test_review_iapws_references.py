#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 12:23:16 2021

@author: unknown
"""

import iapws
import changevar
import useful_tools as ut

rho = 0.4521
T = 460
xmol = 1.0

# print(iapws.ammonia.H2ONH3()._prop(rho, T, xmol))


P = 101325
MW_NH3 = 17.03026 # [g/mol]     Referenced number from iapws, inside nh3h2o documentation.
MW_H2O = 18.015268 # [g/mol]    Referenced number from iapws, inside nh3h2o documentation.

xmass = ut.xmol_to_xmass(xmol, MW_NH3, MW_H2O)

print("\n")
print(changevar.my_prop_bissection(P, T, xmass, -0.1, 25.1, 10, pow(10,-7), 30))
