from __future__ import division

import math
import matplotlib as mpl
import numpy as np


from sympy import Symbol
import mpmath as mp


from scipy.misc import comb
from scipy.misc import factorial 
from scipy.special import gamma

from scipy.special import zeta
from scipy.misc import derivative



def riemann_siegel_psi(p):	
	return mp.cos( 2*np.pi * (p**2 - p - 1/16)) / mp.cos( 2*np.pi*p )



def get_riemann_siegel_C_s(p, order):
	
	C_s = []


	if order >= 0:
		C_0 = riemann_siegel_psi(p)
		C_s.append(C_0)
	
	if order >= 1:
		psi_factor = mp.diff(riemann_siegel_psi, x=p, n=3)
		C_1 = - (1 / ((2**5) * 3 * (np.pi**2))) * psi_factor
		C_s.append(C_1)

	if order >= 2:
		psi_factor_1 = mp.diff(riemann_siegel_psi, x=p, n=6)
		psi_factor_2 = mp.diff(riemann_siegel_psi, x=p, n=2)

		C_2 = ( 1 / ((2**11) * (3**2) * (np.pi**4))) * psi_factor_1 + (1 / ((2**6) * (np.pi**2))) * psi_factor_2
		C_s.append(C_2)

	if order >= 3:
		psi_factor_1 = mp.diff(riemann_siegel_psi, x=p, n=9)
		psi_factor_2 = mp.diff(riemann_siegel_psi, x=p, n=5)
		psi_factor_3 = mp.diff(riemann_siegel_psi, x=p, n=1)

		C_3 = (1 / ( (2**16) * (3**4) * (np.pi**6))) * psi_factor_1 - (1 / ( (2**8) * 3 * 5 * (np.pi**4))) * psi_factor_2 + (1 / ( (2**6) * (np.pi**2))) * psi_factor_3
		C_s.append(C_3)

	if order >= 4:
		psi_factor_1 = mp.diff(riemann_siegel_psi, x=p, n=12)
		psi_factor_2 = mp.diff(riemann_siegel_psi, x=p, n=8)
		psi_factor_3 = mp.diff(riemann_siegel_psi, x=p, n=4)
		psi_factor_4 = riemann_siegel_psi(p)

		C_4 = (1 / ( (2**23) * (3**5) * (np.pi**8))) * psi_factor_1 + (11 / ( (2**17) * (3**2) * 5 * (np.pi**6))) * psi_factor_2 + (19 / ( (2**13) * 3 * (np.pi**4))) * psi_factor_3 + (1 / ( (2**7) * (np.pi**2))) * psi_factor_4
		C_s.append(C_4)

	while len(C_s) < 5:
		C_s.append(0)


	print C_s[0]
	return C_s



def riemann_siegel_remainder(t, order=0):
	R = 0.0

	f = np.sqrt(t / (2 * np.pi))

	N, p = int(f), abs( f - int(f) )	# The integer and fractional part of sqrt( t / 2pi ).

	C_s = get_riemann_siegel_C_s(p, order)

	R_s = []
	R_s.append( C_s[0] )
	R_s.append( C_s[1] * ((t / (2*np.pi))**(-1/2)) )
	R_s.append( C_s[2] * ((t / (2*np.pi))**(-2/2)) )
	R_s.append( C_s[2] * ((t / (2*np.pi))**(-3/2)) )
	R_s.append( C_s[2] * ((t / (2*np.pi))**(-4/2)) )

	R_s = [ C_s[i] * ((t / (2*np.pi))**(-i/2)) for i in range(0, 5)]

	print
	print "The C's are:"
	i = 0
	for term in C_s:
		print "C{} = {}".format(i, float(term))
		i += 1

	print
	

	print "The corrected C's are:"
	i = 0
	for term in C_s:
		print "C{} = {}".format(i, float( ((t / (2*np.pi))**(-1/4)) * term))
		i += 1

	print

	print "The R's are: "
	i = 0
	for term in R_s:
		print "R{} = {}".format(i, float(term))

	print 

	factor = (f**(-1/2)) * sum(R_s)
	
	print "\nCorrection = {}".format(factor)

	R = ((-1)**(N - 1)) * factor

	return R

def riemann_siegel_theta(t):
	first_term = np.angle( gamma( (2.j*t + 1) / 4) )
	second_term = t * np.log(np.pi) / 2
	
	return first_term - second_term

def z_function(t):

	upper_limit = math.floor( np.sqrt(t / (2 * np.pi)) )

	
	z = 0.0
	k = 1
	
	while k**2 < t / (2*np.pi):
		# print "doing"
		z += (k**(-0.5)) * np.cos( riemann_siegel_theta(t) - t * np.log(k) )
		k += 1

	z = 2*z

	zeta = z / np.exp( riemann_siegel_theta(t) * 1.j )

	remainder = riemann_siegel_remainder(t, order=4)

	print "Without R, z = {}".format(z)

	z = z + remainder

	print "With R = {}, z = {}".format(remainder, z)
	
	return z


def zeta_function(t):
	z = z_function(t)
	return z


	
print zeta_function(t=25)


ps =  riemann_siegel_psi(0.692569)
cor = ((18/(2*np.pi))**(-1/4)) * ps

print ps, cor

# psi = lambda p: mp.cos( 2*mp.pi * (p**2 - p - 1/16)) / mp.cos( 2*mp.pi*p )


# print psi(5)


# print mp.diff(psi, x=2, n=1)
# print mp.diff(riemann_siegel_psi, x=2, n=1)



# psi = lambda p: mp.cos( 2*np.pi * (p**2 - p - 1/16)) / mp.cos( 2*np.pi*p )
# psi = riemann_siegel_psi
# print mp.diff(psi, x=0.692, n=3)