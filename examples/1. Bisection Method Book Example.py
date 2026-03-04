# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

def f(num1):
    return pow(num1,3)+4*pow(num1,2)+num1-10

N_0 = 30000
Tol = pow(10,-15)
a = 1
b = 2

i = 1
while i <= N_0:
    p = a + (b - a)/2
    if f(p) == 0 or abs((b-a)/2) < Tol:
        break
    elif f(a)*f(p) > 0:
        a = p
    else:
        b = p
    i += 1

if i > N_0:
    print("Method failed after " + str(i-1) + " iteration(s).")
else:
    print("p = "+ str(p)+", after "+str(i)+" iterations.")
    print("f(p) = "+str(f(p)))










