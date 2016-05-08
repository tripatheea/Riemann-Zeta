from __future__ import division

import math
import matplotlib as mpl
import numpy as np



import sympy as sp



from scipy.misc import comb
from scipy.misc import factorial 
from scipy.special import gamma

from scipy.special import zeta
from scipy.misc import derivative


def read_riemann_siegel_coeffs(filename):
	f = open(filename, 'r')

	lines = f.read().split("\n")

	all_coeffs = []
	for line in lines:
		if len(line) != 0:
			coeffs = [ float(c) for c in line.split(":")[1].split(",")]

			all_coeffs.append(coeffs)

	return all_coeffs

def evaluate_riemann_siegel_C_s(p):
	all_coeffs = read_riemann_siegel_coeffs('data/coeffs.dat')
	results = []

	for k in range(0, 4):
		result = 0.0
		for i in range(len(all_coeffs[k])):
			result += all_coeffs[k][i] * p**i

		results.append(result)

	return results
			



def riemann_siegel_psi(p):	
	return sp.cos( 2*np.pi * (p**2 - p - 1/16)) / sp.cos( 2*np.pi*p )



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



def numeric_riemann_siegel_remainder(t):
	R = 0.0
	f = np.sqrt(t / (2 * np.pi))
	N, p = int(f), abs( f - int(f) )	# The integer and fractional part of sqrt( t / 2pi ).

	C_s = evaluate_riemann_siegel_C_s(p)

def riemann_siegel_remainder(t, order=0):
	R = 0.0

	f = np.sqrt(t / (2 * np.pi))

	N, p = int(f), abs( f - int(f) )	# The integer and fractional part of sqrt( t / 2pi ).

	C_s = get_riemann_siegel_C_s(p, order)

	all_coeffs = evaluate_riemann_siegel_C_s(p)

	print C_s

	print all_coeffs

	# R_s = []
	# R_s.append( C_s[0] )
	# R_s.append( C_s[1] * ((t / (2*np.pi))**(-1/2)) )
	# R_s.append( C_s[2] * ((t / (2*np.pi))**(-2/2)) )
	# R_s.append( C_s[2] * ((t / (2*np.pi))**(-3/2)) )
	# R_s.append( C_s[2] * ((t / (2*np.pi))**(-4/2)) )

	# R_s = [ ((-1)**(N - 1)) * (f**(-1/2)) * C_s[i] * ((t / (2*np.pi))**(-i/2)) for i in range(0, 5)]

	# print "The R's are:"
	# for i in range(len(R_s)):
	# 	print "R{} = {}".format(i, R_s[i])


	# return sum(R_s)

	return 5.0

def riemann_siegel_theta(t):
	first_term = np.angle( gamma( (2.j*t + 1) / 4) )
	second_term = t * np.log(np.pi) / 2
	
	return first_term - second_term

def z_function(t):

	upper_limit = math.floor( np.sqrt(t / (2 * np.pi)) )

	z = 0.0
	k = 1
	while k**2 < t / (2*np.pi):		
		z += (k**(-0.5)) * np.cos( riemann_siegel_theta(t) - t * np.log(k) )
		k += 1

	z *= 2

	remainder = riemann_siegel_remainder(t, order=4)

	z += remainder

	print "Z({}) = {}".format(t, z)
	
	return z


def zeta_function(t):
	z = z_function(t)
	return z



if __name__ == '__main__':
	zeta_function(t=25)
