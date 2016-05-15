from __future__ import division

import math
import numpy as np

from time import time

import sympy as sp

import mpmath as mp

from mpmath.ctx_mp_python import mpf

from scipy.misc import factorial 
from scipy.special import gamma

precision = 53
mp.prec = precision
mp.pretty = True


def calculate_factorial_ratio(n, i):
	# This function calculates (n + i - 1)! / (n - i)!
	
	mp.dps = 50

	k = (n - i)

	result = 1
	for j in range(k + 2*i - 1, k, -1):
		result = mp.fmul(result, j)

	return result

def dirichlet_eta(s, N):
	'''
		This is essentially the Euler transformation applied to Dirichlet eta. 
	'''

	mp.dps = 50

	def calculate_d_k(k, n):
		d_k_total = mpf(0.0)
		for i in range(k):
			# numerator = factorial( ( n + i - 1) ) * 4**i
			# denominator = factorial( n - i ) * factorial( 2*i )
			# summand = numerator / denominator
			summand = (calculate_factorial_ratio(n, i) * (4**i)) / factorial(2 * i)
			
			d_k_total += summand

		return n * d_k_total

	d_n = calculate_d_k(N, N)

	eta_total = 0.0
	for k in range(N):
		d_k = calculate_d_k(k, N)
		numerator = ( (-1)**k ) * (d_k - d_n)
		denominator = (k + 1)**s
		eta_total += numerator / denominator

	return  - (1 / d_n) * eta_total 

def alternating_series(s, N):
	eta = dirichlet_eta(s, N)
	denominator = 1 - 2**(1 - s)
	zeta = eta / denominator

	return zeta

def riemann_siegel_theta(t):
	first_term = np.angle( gamma( (2.j*t + 1) / 4) )
	second_term = t * np.log(np.pi) / 2
	
	return first_term - second_term

def zeta_function(s, N):
	z = alternating_series(s, N)

	return z

def z_function(t, N=100000):
	zeta = zeta_function(1/2 + (1.j)*t, N)

	return mp.re( np.exp( 1.j * riemann_siegel_theta(t) ) * zeta )

def calculate_z(t):	# Convenient wrapper to use for roots.py
	return z_function(t, N=25)

if __name__ == '__main__':

	start = time()

	# print zeta_function(s=1/2 + 25.j, N=1000)
	# print z_function(t=18, N=100)

	eta = dirichlet_eta(1, N=25)
	print eta
	print abs(eta - np.log(2))

	end = time()

	print "Calculated using alternating series in {:.4f} seconds.".format(float(end - start))
