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

def n_choose_k(n, k):
	
	j = n - k

	numerator = 1
	for i in range(1, k + 1):
		numerator *= (j + i)

	denominator = factorial(k)

	return numerator / denominator

def dirichlet_eta(s, N):

	def calculate_d_n(n):
		total = 0.0
		for k in range(n + 1):
			if k % 2 == 0:
				alternating_factor = 1
			else:
				alternating_factor = -1

			total += alternating_factor * n_choose_k(n, k) / ( k + 1)**s
		return total

	eta = 0.0
	for n in range(N + 1):
		d_n = calculate_d_n(n)
		eta += d_n / (2**(n + 1))

	return eta





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


	# print zeta_function(s=1/2 + 25.j, N=1000)
	# print z_function(t=18, N=100)

	start = time()

	eta = dirichlet_eta(1, N=25)
	print eta
	print abs(eta - np.log(2))

	end = time()

	print "Calculated using alternating series in {:.4f} seconds.".format(float(end - start))
