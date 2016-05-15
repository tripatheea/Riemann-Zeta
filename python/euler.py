from __future__ import division

import math
import matplotlib as mpl
import numpy as np

from time import time

from scipy.misc import comb
from scipy.misc import factorial 
from scipy.special import gamma

import mpmath as mp


def brute_force_sum(N, s):
	total = 0.0
	for n in range(1, N):
		total += np.power(n, -s)
	return total


def calculate_bernoulli_numbers():
	'''
		Calculates the nth Bernoulli number B_n. Uses the B_1 = 1/2 convention.
	'''

	n_max = 20
	bernoulli_numbers = [1]

	for m in range(1, n_max + 1):
		B_m = 0
		for k in range(0, m):
			B_m += comb(m, k, exact=True) * bernoulli_numbers[k] / (m - k + 1)

		B_m = 1 - B_m
		bernoulli_numbers.append( B_m )

	return bernoulli_numbers


def bernoulli_sum(N, s):
	total = 0.0

	bernoulli_numbers = calculate_bernoulli_numbers()

	v_max = 10

	for v in range(1, v_max + 1):
		B_2v = bernoulli_numbers[ 2 * v ]
		fac = factorial(2 * v)

		s_product = s
		for i in range(1, (2*v - 2) + 1 ):
			s_product *= s + i

		N_term = N**( - s - 2*v + 1)

		total += B_2v * s_product * N_term / fac

	return total



def riem_sieg_gamma_function(s):
	numerator = gamma( s / 2 )
	denominator = gamma( (1 - s) / 2 )
	return (np.pi**(0.5 - s)) * (numerator / denominator)


def riemann_siegel_theta(t):
	first_term = np.angle( gamma( (2.j*t + 1) / 4) )
	second_term = t * np.log(np.pi) / 2
	
	return first_term - second_term




def zeta_function(s, N):

	z = brute_force_sum(N=N, s=s)
	z += bernoulli_sum(N=N, s=s)
	z += (N**(1 - s)) / (s - 1)
	z += 0.5 * N**(-s)
	
	return z


def z_function(t, N=100000):
	zeta = zeta_function(1/2 + (1.j)*t, N)

	return mp.re( np.exp( 1.j * riemann_siegel_theta(t) ) * zeta )

def calculate_z(t):	# Convenient wrapper to use for roots.py
	return z_function(t, N=1000)


if __name__ == '__main__':
	start = time()

	# print zeta_function(s=(1/2 + 25.j), N=100000, method="euler")
	print z_function(t=18, N=100000)

	end = time()

	print "Calculated using Euler's method in {:.4f} seconds.".format(float(end - start))


	
