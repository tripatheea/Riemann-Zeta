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
		coeff = fmul(mpf(original_coeffs[i + m] * (2**(i + m))), mpf(calculate_factorial_ratio(i, m))) / 2**i

		coefficients_of_derivative.append( coeff )

	return coefficients_of_derivative

def write_coeffs(filename, coeffs, m):
	f = open(filename, "a")

	f.write("m={}:".format(m))
	for i in range(len(coeffs)):
		if i < len(coeffs) - 1:
			f.write("{},".format(coeffs[i]))
		else:
			f.write("{}\n".format(coeffs[i]))

	f.close()

def get_coefficients(function, p, x0, n):

	taylor_series = sp.series(expr=function, x=p, x0=x0, n=n).removeO()
	
	coefficients_dict_sympy_pow_keys = taylor_series.as_coefficients_dict()

	coefficients_dict_string_keys = {}
	for key in coefficients_dict_sympy_pow_keys.keys():
		coefficients_dict_string_keys[str(key)] = coefficients_dict_sympy_pow_keys[key]

	taylor_coefficients = []
	for i in range(n):
		if i == 0:
			key = '1'
		elif i == 1:
			key = 'p'
		else:
			key = "(p - 0.5)**{}".format(i)
		taylor_coefficients.append( coefficients_dict_string_keys[key] / (2**i) )


	coeffs = [ sp.N(c, precision) for c in taylor_coefficients ]

	return coeffs


def get_taylor_series(n, filename):

	f = open(filename, "w")
	f.close()

	p = sp.symbols('p')

	psi_0 = sp.cos( 2*mp.pi * (p**2 - p - 1/16)) / sp.cos( 2*mp.pi*p )

	C_0 = psi_0
	coeffs =  get_coefficients( C_0, p, 1/2, n)
	write_coeffs( filename=filename, coeffs=coeffs, m=0 )

	psi_factor = get_coefficients_of_derivative(m=3, original_coeffs=coeffs)
	C_1 = [ mpf( (mpf(1) / mpf( (2**5) * 3 * (mp.pi**2) )) * coeff) for coeff in psi_factor ]
	write_coeffs( filename=filename, coeffs=C_1, m=1 )

	psi_factor_1 = get_coefficients_of_derivative(m=2, original_coeffs=coeffs)
	psi_factor_2 = get_coefficients_of_derivative(m=6, original_coeffs=coeffs)
	C_2 = [ (1 / ((2**6) * (mp.pi**2))) * psi_1 + ( 1 / ((2**11) * (3**2) * (mp.pi**4))) * psi_2 for (psi_1, psi_2) in zip( psi_factor_1, psi_factor_2) ]
	write_coeffs( filename=filename, coeffs=C_2, m=2 )

	
	psi_factor_1 = get_coefficients_of_derivative(m=1, original_coeffs=coeffs)
	psi_factor_2 = get_coefficients_of_derivative(m=5, original_coeffs=coeffs)
	psi_factor_3 = get_coefficients_of_derivative(m=9, original_coeffs=coeffs)
	C_3 = [ (1 / ( (2**6) * (mp.pi**2))) * psi_1 + (1 / ( (2**8) * 3 * 5 * (mp.pi**4))) * psi_2 + (1 / ( (2**16) * (3**4) * (mp.pi**6))) * psi_3 for (psi_1, psi_2, psi_3) in zip(psi_factor_1, psi_factor_2, psi_factor_3) ]
	write_coeffs( filename=filename, coeffs=C_3, m=3 )

	
	psi_factor_1 = coeffs
	psi_factor_2 = get_coefficients_of_derivative(m=4, original_coeffs=coeffs)
	psi_factor_3 = get_coefficients_of_derivative(m=8, original_coeffs=coeffs)
	psi_factor_4 = get_coefficients_of_derivative(m=12, original_coeffs=coeffs)
	C_4 = [ (1 / ( (2**7) * (mp.pi**2))) * psi_1 + (19 / ( (2**13) * 3 * (mp.pi**4))) * psi_2 + (11 / ( (2**17) * (3**2) * 5 * (mp.pi**6))) * psi_3 + (1 / ( (2**23) * (3**5) * (mp.pi**8))) * psi_4 for (psi_1, psi_2, psi_3, psi_4) in zip(psi_factor_1, psi_factor_2, psi_factor_3, psi_factor_4) ]
	write_coeffs( filename=filename, coeffs=C_4, m=4 )
	





# if __name__ == '__main__':
	# get_taylor_series(n=50, filename="data/coeffs.dat")



