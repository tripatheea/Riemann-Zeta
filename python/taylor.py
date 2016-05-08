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

def get_taylor_series(n):


	n = 5
	# mp.dps = 100
	mp.prec = 500
	mp.pretty = True

	print mp

	
	a = mp.sin
	# print mp.taylor(a, 0, 5)


	p = sp.symbols('p')

	taylor_series = sp.series(expr=sp.cos( 2*np.pi * (p**2 - p - 1/16)) / sp.cos( 2*np.pi*p ), x=p, x0=(0.5), n=n).removeO()
	
	taylor_series_polynomial =  sp.Poly(taylor_series, (p))
	# print taylor_series_polynomial.all_coeffs()
	
	print "\n"*2	

	print taylor_series_polynomial.all_coeffs()[0].evalf(54)

	

	f = implemented_function(Function('f'), lambda x: x+mp.pi)
	lam_f = lambdify(x, f(x))
	# print mpf( lam_f(4) )

	

	psi_0 = implemented_function( Function('psi_0'), lambda x: (mp.cos( 2 *mp.pi*(x*x - x - 1/16) ) / mp.cos(2*mp.pi*x)) )

	lam_psi_0 = lambdify(x, psi_0(x))

	
	taylor_coeffs = mp.taylor( psi_0, 0.5, n)

	print taylor_coeffs[0]


	# print "\n"*2
	# The coefficients in taylor_coeffs will have psi_0's so we need to evaluate them.

	evaluated_coeffs = []
	for coeff in taylor_coeffs:

		
		# print sp.N( coeff )

		if type(coeff) is Mul:
			mpf_coeff = coeff.as_coeff_Mul()
			c = lambdify(x, mpf_coeff)
			# coe = mpf_coeff, c(0)
			coe = reduce(mul, list(c(0)), 1)

			# print coe

		elif type(coeff) is Add:
			mpf_coeff = coeff.as_coeff_Add()
			c = lambdify(x, mpf_coeff)
			# print mpf_coeff
			coe = sum(list(c(0)))
			# print coe
		elif type(coeff) is mpf:
			# evaluated_coeffs.append( coeff )
			coe = coeff
		else:
			raise Exception("Error: Unknown data type found!")
		
		evaluated_coeffs.append( coe )
	

	print evaluated_coeffs


	'''
	a =   taylor_coeffs[0].as_coeff_Mul()[1]

	print a

	b = lambdify(x, a )
	print b

	print b(100)
	'''

def get_riemann_siegel_C_s(p_value, order):
	
	C_s = []

	p = sp.symbols('p')
	psi_0 = sp.cos( 2*np.pi * (p**2 - p - 1/16)) / sp.cos( 2*np.pi*p )

	if order >= 0:
		C_0 = psi_0.evalf(subs={p:p_value})
		C_s.append(C_0)
	
	if order >= 1:
		psi_factor = sp.diff(psi_0, p, 3).evalf(subs={p:p_value})
		C_1 = - (1 / ( (2**5) * 3 * (np.pi**2) )) * psi_factor
		C_s.append(C_1)
	
	if order >= 2:
		psi_factor_1 = sp.diff(psi_0, p, 2).evalf(subs={p:p_value})
		psi_factor_2 = sp.diff(psi_0, p, 6).evalf(subs={p:p_value})

		C_2 = (1 / ((2**6) * (np.pi**2))) * psi_factor_1 + ( 1 / ((2**11) * (3**2) * (np.pi**4))) * psi_factor_2
		C_s.append(C_2)

	if order >= 3:
		psi_factor_1 = sp.diff(psi_0, p, 1).evalf(subs={p:p_value})
		psi_factor_2 = sp.diff(psi_0, p, 5).evalf(subs={p:p_value})
		psi_factor_3 = sp.diff(psi_0, p, 9).evalf(subs={p:p_value})

		C_3 = - (1 / ( (2**6) * (np.pi**2))) * psi_factor_1 - (1 / ( (2**8) * 3 * 5 * (np.pi**4))) * psi_factor_2 -  (1 / ( (2**16) * (3**4) * (np.pi**6))) * psi_factor_3
		C_s.append(C_3)

	if order >= 4:
		psi_factor_1 = psi_0.evalf(subs={p:p_value})
		psi_factor_2 = sp.diff(psi_0, p,  4).evalf(subs={p:p_value})
		psi_factor_3 = sp.diff(psi_0, p,  8).evalf(subs={p:p_value})
		psi_factor_4 = sp.diff(psi_0, p, 12).evalf(subs={p:p_value})
		
		C_4 = (1 / ( (2**7) * (np.pi**2))) * psi_factor_1 + (19 / ( (2**13) * 3 * (np.pi**4))) * psi_factor_2 + (11 / ( (2**17) * (3**2) * 5 * (np.pi**6))) * psi_factor_3 + (1 / ( (2**23) * (3**5) * (np.pi**8))) * psi_factor_4
		C_s.append(C_4)
	
	while len(C_s) < 5:
		C_s.append(0)

	return C_s






if __name__ == '__main__':
	get_taylor_series(n=0)
