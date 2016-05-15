from __future__ import division

import math

import numpy as np

import sympy as sp

from scipy.special import gamma


def read_riemann_siegel_coeffs(filename):
	f = open(filename, 'r')

	lines = f.read().split("\n")

	all_coeffs = []
	for line in lines:
		if len(line) != 0:
			coeffs = [ float(c) for c in line.split(":")[1].split(",")]

			all_coeffs.append(coeffs)

	return all_coeffs

def evaluate_riemann_siegel_C_series(p, order):

	all_coeffs = read_riemann_siegel_coeffs('data/coeffs.dat')
	results = []

	if order >= 0:
		results.append( sum( [ all_coeffs[0][i] * p**i for i in range(len(all_coeffs[0])) ] ) )

	if order >= 1:
		results.append( sum( [ all_coeffs[1][i] * p**i for i in range(len(all_coeffs[1])) ] ) )

	if order >= 2:
		results.append( sum( [ all_coeffs[2][i] * p**i for i in range(len(all_coeffs[2])) ] ) )

	if order >= 3:
		results.append( sum( [ all_coeffs[3][i] * p**i for i in range(len(all_coeffs[3])) ] ) )

	if order >= 4:
		results.append( sum( [ all_coeffs[4][i] * p**i for i in range(len(all_coeffs[4])) ] ) )

	while len(results) < 5:
		results.append(0)

	return results
			

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

def riemann_siegel_remainder(t, order=0, use_series_expansion=True):
	
	# f = sp.sqrt(t / (2 * sp.pi))
	f = np.sqrt(t / (2 * np.pi))

	N, p = int(f), abs( f - int(f) )	# The integer and fractional part of sqrt( t / 2pi ).

	if use_series_expansion:
		C_s = evaluate_riemann_siegel_C_series(1 - 2*p, order)
	else:
		C_s = get_riemann_siegel_C_s(p, order)
	
	R_s = [ (f**(-1/2))*C_s[i] * f**(- i ) for i in range(len(C_s)) ]
	
	if (N - 1) % 2 == 1:
		return - sum(R_s)
	
	return sum(R_s)
	
def riemann_siegel_theta(t):
	first_term = np.angle( gamma( (2.j*t + 1) / 4) )
	second_term = t * np.log(np.pi) / 2
	
	return first_term - second_term

def z_function(t, remainder_order=4, use_series_expansion=True):

	upper_limit = math.floor( np.sqrt(t / (2 * np.pi)) )

	z = 0.0
	k = 1
	while k**2 < t / (2*np.pi):		
		z += (k**(-0.5)) * np.cos( riemann_siegel_theta(t) - t * np.log(k) )
		k += 1

	z *= 2

	remainder = riemann_siegel_remainder(t, order=remainder_order, use_series_expansion=use_series_expansion)

	z += remainder

	# print "Z({}) = {}".format(t, z)
	
	return z

def calculate_z(t):	# Convenient wrapper to use for roots.py
	return z_function(t, remainder_order=4, use_series_expansion=True)

def zeta_function(t):
	z = z_function(t)
	return z



if __name__ == '__main__':
	z_function(t=18, remainder_order=4, use_series_expansion=True)
