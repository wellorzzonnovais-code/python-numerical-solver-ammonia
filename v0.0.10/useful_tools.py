# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""
Version 0.0.1 to 0.0.6
Some changes (not described)

Version 0.0.7

Added xmol_to_xmass function
"""

# Tools for conversions

def xmass_to_xmol(xmass,mass_weight1,mass_weight2):
    xmol_converted = (xmass/mass_weight1)/(xmass/mass_weight1 + (1-xmass)/mass_weight2)
    return xmol_converted

def xmol_to_xmass(xmol,mass_weight1,mass_weight2):
    xmass_converted = xmol*mass_weight1/(xmol*mass_weight1 + (1-xmol)*mass_weight2)
    return xmass_converted

pressure_MPa_to_Pa = pow(10,6)

pressure_bar_to_Pa = pow(10,5)

pressure_atm_to_Pa = 101325
