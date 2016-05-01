from __future__ import division

import math
import matplotlib as mpl
import numpy as np

from scipy.misc import comb
from scipy.misc import factorial 


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


def alternating_series(s, N):
	'''
		This is essentially the Euler transformation applied to Dirichlet eta. 
	'''

	def e_k(k):
		total = 0.0
		for j in range(k, N + 1):
			total += comb(N, j)
		return total


	prefactor = 1 / ( 1 - 2**(1 - s))


	first_term = 0.0	

	for k in range(1, N + 1):
		first_term += ( (-1)**k ) / (k**s)

	second_term = 0.0
	for k in range( N + 1, 2*N + 1):
		second_term +=  (1 / (2**N)) * ( (-1)**k ) * e_k(k) / (k**s)

	zeta = first_term + second_term

	zeta *= prefactor

	return zeta


def zeta_function(s, N, method="euler"):

	z = 0.0

	if method == "euler":
		z = brute_force_sum(N=N, s=s)

		z += bernoulli_sum(N=N, s=s)

		z += (N**(1 - s)) / (s - 1)

		z += 0.5 * N**(-s)

	elif method == "alternating":
		z = alternating_series(s, N)
		return z

	return z
	


print zeta_function(s=1/2, N=4, method="euler")
print zeta_function(s=1/2, N=50, method="alternating")












# def alternating_series(s, N):
# 	'''
# 		This is essentially the Euler transformation applied to Dirichlet eta. 
# 	'''

# 	def e_k(k):
# 		total = 0.0
# 		for j in range(k, N + 1):
# 			total += comb(N, j)
# 		return total


# 	prefactor = 1 / ( 1 - 2**(1 - s))

# 	# first_term = 0.0

# 	zeta = 0.0

# 	for k in range(1, N + 1):
# 		first_term = ( (-1)**k ) / (k**s)

# 		k2 = k + 1
# 		second_term = ( (-1)**k2 ) * e_k(k2) / (k2**s)
# 		second_term *= 1 / (2**N)

# 		zeta += first_term + second_term

# 	zeta *= prefactor

# 	return zeta