from __future__ import division

import math
import numpy as np

from time import time

from scipy.misc import comb
from scipy.misc import factorial 
from scipy.special import gamma


def dirichlet_eta(s, N):
	'''
		This is essentially the Euler transformation applied to Dirichlet eta. 
	'''

	def calculate_d_k(k, n):
		d_k_total = 0.0
		for i in range(k):
			numerator = factorial( ( n + i - 1) ) * 4**i
			denominator = factorial( n - i ) * factorial( 2*i )
			summand = numerator / denominator
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

	return np.real( np.exp( 1.j * riemann_siegel_theta(t) ) * zeta )

if __name__ == '__main__':

	start = time()

	# print zeta_function(s=1/2 + 25.j, N=1000)
	print z_function(t=18, N=100)

	end = time()

	print "Calculated using alternating series in {:.4f} seconds.".format(float(end - start))
