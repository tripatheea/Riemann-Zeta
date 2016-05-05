from __future__ import division

import math
import matplotlib as mpl
import numpy as np

from time import time

from scipy.misc import comb
from scipy.misc import factorial 
from scipy.special import gamma

import zeta



def next_upper_bracket(t):
	current_sign = np.sign( zeta.z_function(t, N=100) )
	next_sign = current_sign

	while current_sign == next_sign:
		t += 1
		next_sign = np.sign( zeta.z_function(t, N=100) )

	return t

def find_root(t):
	# Finds root at the point nearest >= t.
	n_max = 10000

	eps = 1e-9
	lower = t

	# Need to find the lowest upper bracketing value.
	upper = next_upper_bracket(t)

	n = 0
	while n <= n_max:
		mid = (upper + lower) / 2
		value = zeta.z_function(mid)

		# print "Trying between {} and {}".format(lower, upper)
		# print "Got Z({}) = {}".format(mid, value)

		if abs(value) < eps: # Found a root.
			print "The root closest to t={} in the upward direction is {}.".format(t, mid)
			return mid
		else:
			n += 1
			if np.sign(zeta.z_function(lower)) == np.sign(zeta.z_function(mid)):
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



# find_root(5)

write_roots("roots.dat")