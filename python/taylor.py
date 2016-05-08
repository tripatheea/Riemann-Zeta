from __future__ import division

import math
import matplotlib as mpl
import numpy as np



import sympy as sp
from mpmath import *

import inspect

from scipy.misc import comb
from scipy.misc import factorial 
from scipy.special import gamma

from scipy.special import zeta
from scipy.misc import derivative

from sympy.abc import x
from sympy.utilities.lambdify import lambdify, implemented_function
from sympy import Function

from sympy.core.mul import Mul
from sympy.core.add import Add
from mpmath.ctx_mp_python import mpf

from operator import mul
from time import time


precision = 53

mp.prec = precision

mp.pretty = True



def calculate_factorial_ratio(k, m):
	# This function calculates (k + m)! / k!
	# This is much faster than directly computing (k + m)! / k! because this function computes only (m - 1) products, 
	# compared to (k + m)! / k!, which has to comput 2k + m products.
	
	mp.dps = 50

	result = 1
	for i in range(k + m, k, -1):
		result = fmul(result, i)

	return result


def get_coefficients_of_derivative(m, original_coeffs):
	# original_coeffs holds the Taylor series coefficients of our original function.
	# We want the kth coefficient of the mth deriviative of our original function.
	# That is, we want the coefficient of the x^k term of the d^m f(x) / dx^m if f(x) is our original function.

	# We want the nth coefficient of the original function. 
	# This is well-defined because the least values of k and m are both 0 so n >= 0 and so, n here is always well-defined.

	coefficients_of_derivative = []

	terms_in_original_function = len(original_coeffs)	# Say this is N.
	
	# We can get coefficients of only N - m  terms.
	for i in range(0, terms_in_original_function - m):
		coeff = fmul(mpf(original_coeffs[i + m]), mpf(calculate_factorial_ratio(i, m)))

		coefficients_of_derivative.append( coeff )

	return coefficients_of_derivative


def get_coefficients(function, p, x0, n):
	taylor_series = sp.series(expr=function, x=p, x0=x0, n=n).removeO()
	taylor_series_polynomial =  sp.Poly(taylor_series, (p))
	coeffs = [ sp.N(c, precision) for c in taylor_series_polynomial.all_coeffs() ][::-1]

	return coeffs

def write_coeffs(filename, coeffs, m):
	f = open(filename, "a")

	f.write("m={}:".format(m))
	for i in range(len(coeffs)):
		if i < len(coeffs) - 1:
			f.write("{},".format(coeffs[i]))
		else:
			f.write("{}\n".format(coeffs[i]))

	f.close()



def get_taylor_series(n, filename):


	f = open(filename, "w")
	f.close()

	p = sp.symbols('p')

	psi_0 = sp.cos( 2*mp.pi * (p**2 - p - 1/16)) / sp.cos( 2*mp.pi*p )

	C_0 = psi_0
	coeffs =  get_coefficients( C_0, p, 0, n)
	write_coeffs( filename=filename, coeffs=coeffs, m=0 )


	psi_factor = get_coefficients_of_derivative(m=3, original_coeffs=coeffs)
	C_1 = [ mpf( - (mpf(1) / mpf( (2**5) * 3 * (mp.pi**2) )) * coeff) for coeff in psi_factor ]
	write_coeffs( filename=filename, coeffs=C_1, m=1 )


	psi_factor_1 = get_coefficients_of_derivative(m=2, original_coeffs=coeffs)
	psi_factor_2 = get_coefficients_of_derivative(m=6, original_coeffs=coeffs)
	C_2 = [ (1 / ((2**6) * (mp.pi**2))) * psi_1 + ( 1 / ((2**11) * (3**2) * (mp.pi**4))) * psi_2 for (psi_1, psi_2) in zip( psi_factor_1, psi_factor_2) ]
	write_coeffs( filename=filename, coeffs=C_2, m=2 )


	psi_factor_1 = get_coefficients_of_derivative(m=1, original_coeffs=coeffs)
	psi_factor_2 = get_coefficients_of_derivative(m=5, original_coeffs=coeffs)
	psi_factor_3 = get_coefficients_of_derivative(m=9, original_coeffs=coeffs)
	C_3 = [ - (1 / ( (2**6) * (mp.pi**2))) * psi_1 - (1 / ( (2**8) * 3 * 5 * (mp.pi**4))) * psi_2 -  (1 / ( (2**16) * (3**4) * (mp.pi**6))) * psi_3 for (psi_1, psi_2, psi_3) in zip(psi_factor_1, psi_factor_2, psi_factor_3) ]
	write_coeffs( filename=filename, coeffs=C_3, m=3 )

	
	psi_factor_1 = coeffs
	psi_factor_2 = get_coefficients_of_derivative(m=4, original_coeffs=coeffs)
	psi_factor_3 = get_coefficients_of_derivative(m=8, original_coeffs=coeffs)
	psi_factor_4 = get_coefficients_of_derivative(m=12, original_coeffs=coeffs)
	C_4 = [ (1 / ( (2**7) * (mp.pi**2))) * psi_1 + (19 / ( (2**13) * 3 * (mp.pi**4))) * psi_2 + (11 / ( (2**17) * (3**2) * 5 * (mp.pi**6))) * psi_3 + (1 / ( (2**23) * (3**5) * (mp.pi**8))) * psi_4 for (psi_1, psi_2, psi_3, psi_4) in zip(psi_factor_1, psi_factor_2, psi_factor_3, psi_factor_4) ]
	write_coeffs( filename=filename, coeffs=C_4, m=4 )
	

	'''

	elif k == 1:	
		psi_factor = sp.diff(psi_0, p, 3)
		C_1 = - (1 / ( (2**5) * 3 * (mp.pi**2) )) * psi_factor
		# all_coeffs.append( get_coefficients( C_1, p, 0, n) )
		write_coeffs( filename=filename, coeffs=get_coefficients( C_1, p, 0, n), k=1 )
	elif k == 2:
		psi_factor_1 = sp.diff(psi_0, p, 2)
		psi_factor_2 = sp.diff(psi_0, p, 6)
		C_2 = (1 / ((2**6) * (mp.pi**2))) * psi_factor_1 + ( 1 / ((2**11) * (3**2) * (mp.pi**4))) * psi_factor_2
		# all_coeffs.append( get_coefficients( C_2, p, 0, n) )
		write_coeffs( filename=filename, coeffs=get_coefficients( C_2, p, 0, n), k=2 )
	elif k == 3:
		psi_factor_1 = sp.diff(psi_0, p, 1)
		psi_factor_2 = sp.diff(psi_0, p, 5)
		psi_factor_3 = sp.diff(psi_0, p, 9)
		C_3 = - (1 / ( (2**6) * (mp.pi**2))) * psi_factor_1 - (1 / ( (2**8) * 3 * 5 * (mp.pi**4))) * psi_factor_2 -  (1 / ( (2**16) * (3**4) * (mp.pi**6))) * psi_factor_3
		# all_coeffs.append( get_coefficients( C_3, p, 0, n) )
		write_coeffs( filename=filename, coeffs=get_coefficients( C_3, p, 0, n), k=3 )
	elif k == 4:
		psi_factor_1 = psi_0
		psi_factor_2 = sp.diff(psi_0, p,  4)
		psi_factor_3 = sp.diff(psi_0, p,  8)
		psi_factor_4 = sp.diff(psi_0, p, 12)
		C_4 = (1 / ( (2**7) * (mp.pi**2))) * psi_factor_1 + (19 / ( (2**13) * 3 * (mp.pi**4))) * psi_factor_2 + (11 / ( (2**17) * (3**2) * 5 * (mp.pi**6))) * psi_factor_3 + (1 / ( (2**23) * (3**5) * (mp.pi**8))) * psi_factor_4
		# all_coeffs.append( get_coefficients( C_4, p, 0, n) )
		write_coeffs( filename=filename, coeffs=get_coefficients( C_4, p, 0, n), k=4 )
	'''

def get_riemann_siegel_C_s(p_value, order):
	
	C_s = []

	p = sp.symbols('p')
	psi_0 = sp.cos( 2*mp.pi * (p**2 - p - 1/16)) / sp.cos( 2*mp.pi*p )

	if order >= 0:
		C_0 = psi_0.evalf(subs={p:p_value})
		C_s.append(C_0)
	
	if order >= 1:
		psi_factor = sp.diff(psi_0, p, 3).evalf(subs={p:p_value})
		C_1 = - (1 / ( (2**5) * 3 * (mp.pi**2) )) * psi_factor
		C_s.append(C_1)
	
	if order >= 2:
		psi_factor_1 = sp.diff(psi_0, p, 2).evalf(subs={p:p_value})
		psi_factor_2 = sp.diff(psi_0, p, 6).evalf(subs={p:p_value})

		C_2 = (1 / ((2**6) * (mp.pi**2))) * psi_factor_1 + ( 1 / ((2**11) * (3**2) * (mp.pi**4))) * psi_factor_2
		C_s.append(C_2)

	if order >= 3:
		psi_factor_1 = sp.diff(psi_0, p, 1).evalf(subs={p:p_value})
		psi_factor_2 = sp.diff(psi_0, p, 5).evalf(subs={p:p_value})
		psi_factor_3 = sp.diff(psi_0, p, 9).evalf(subs={p:p_value})

		C_3 = - (1 / ( (2**6) * (mp.pi**2))) * psi_factor_1 - (1 / ( (2**8) * 3 * 5 * (mp.pi**4))) * psi_factor_2 -  (1 / ( (2**16) * (3**4) * (mp.pi**6))) * psi_factor_3
		C_s.append(C_3)

	if order >= 4:
		psi_factor_1 = psi_0.evalf(subs={p:p_value})
		psi_factor_2 = sp.diff(psi_0, p,  4).evalf(subs={p:p_value})
		psi_factor_3 = sp.diff(psi_0, p,  8).evalf(subs={p:p_value})
		psi_factor_4 = sp.diff(psi_0, p, 12).evalf(subs={p:p_value})
		
		C_4 = (1 / ( (2**7) * (mp.pi**2))) * psi_factor_1 + (19 / ( (2**13) * 3 * (mp.pi**4))) * psi_factor_2 + (11 / ( (2**17) * (3**2) * 5 * (mp.pi**6))) * psi_factor_3 + (1 / ( (2**23) * (3**5) * (mp.pi**8))) * psi_factor_4
		C_s.append(C_4)
	
	while len(C_s) < 5:
		C_s.append(0)

	return C_s






if __name__ == '__main__':
	get_taylor_series(n=50, filename="data/coeffs.dat")



