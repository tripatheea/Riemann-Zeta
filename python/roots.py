from __future__ import division

import math
import matplotlib as mpl
import numpy as np

from time import time

from scipy.misc import comb
from scipy.misc import factorial 
from scipy.special import gamma

import dirichlet, euler, riemann_siegel

import mpmath as mp

from mpmath.ctx_mp_python import mpf

precision = 53
mp.prec = precision
mp.pretty = True


def next_upper_bracket(t, z_function):

	current_sign = np.sign( mp.re( z_function(t) ) )
	next_sign = current_sign

	while current_sign == next_sign:
		t += 1
		next_sign = np.sign( mp.re( z_function(t) ) )

	return t

def find_root(t, method="riemann_siegel"):
	# Finds root at the point nearest >= t.
	n_max = 10000

	eps = 1e-9
	lower = t

	if method == "euler":
		z_function = lambda t: euler.calculate_z(t)
	elif method == "dirichlet":
		z_function = lambda t: dirichlet.calculate_z(t)
	else:
		z_function = lambda t: riemann_siegel.calculate_z(t)

	# Need to find the lowest upper bracketing value.
	upper = next_upper_bracket(t=t, z_function=z_function)

	n = 0
	while n <= n_max:
		mid = (upper + lower) / 2
		value = z_function(mid)

		# print "Trying between {} and {}".format(lower, upper)
		# print "Got Z({}) = {}".format(mid, value)

		if abs(value) < eps: # Found a root.
			print "The root closest to t={} in the upward direction is {}.".format(t, mid)
			return mid
		else:
			n += 1
			if np.sign( mp.re(z_function(lower)) ) == np.sign( mp.re(z_function(mid)) ):
				lower = mid
			else:
				upper = mid

	return None



def write_roots(filename):
	f = open(filename,'w')
	t = 1

	for i in range(5):
		root = find_root(t)
		f.write(str(root) + "\n")
		t = np.ceil(root)

	f.close()



find_root(5, method="euler")

# write_roots("roots.dat")


# Profiling. This later moves to profile.py

start = time()
find_root(10, method="euler")
end = time()

print "Found root with Euler's method in {} seconds.\n\n".format(end - start)

start = time()
find_root(10, method="dirichlet")
end = time()

print "Found root with Dirichlet's eta function in {} seconds.\n\n".format(end - start)

start = time()
find_root(10, method="riemann_siegel")
end = time()

print "Found root with Riemann Siegel formula in {} seconds.\n\n".format(end - start)





# Profiling Ends.