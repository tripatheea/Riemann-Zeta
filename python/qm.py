from __future__ import division

import math
import matplotlib as mpl
import numpy as np


from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FormatStrFormatter

from sets import Set

import sys
import math
from collections import defaultdict

# matplotlib
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm


# Stuff for calculating areas.
from scipy.integrate import simps
from scipy import interpolate
from scipy import optimize

from numpy import trapz


from matplotlib import gridspec

from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png


from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

from scipy.stats import norm
from scipy.stats import gamma
from scipy import arange, array, exp

from scipy.stats import binned_statistic


import sympy as sp



from scipy.misc import comb
from scipy.misc import factorial 
from scipy.special import gamma

from scipy.special import zeta
from scipy.misc import derivative


mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)






def get_eigenvalues(n):
	# Construct creation operator first.
	a_creation = np.zeros((n, n))
	for i in range(n):
		a_creation[i, i - 1] = np.sqrt(i)

	a_anhilation = np.zeros((n, n))
	for i in range(n):
		a_anhilation[i - 1, i] = np.sqrt(i)

	position = a_creation + a_anhilation
	momentum = 1.j * ( a_creation - a_anhilation)

	berry_H = (position * momentum + momentum * position) / 2

	# print berry_H

	eigenvalues = sorted(np.linalg.eigvalsh(berry_H))


	return eigenvalues

def get_eigenvalue_differences(n):
	
	eigenvalues = get_eigenvalues(n)



	normalized_differences = np.diff(eigenvalues)
	normalized_differences *= 1 / np.mean(normalized_differences)

	return normalized_differences


def construct_operators(n):
	
	normalized_differences = get_eigenvalue_differences(n)

	
	plt.hist(normalized_differences, color="red", bins=100, lw=5, histtype='step', edgecolor="red", normed=1)

	plt.autoscale()
	plt.xlim(0, 3)

	plt.gca().xaxis.set_minor_locator(MultipleLocator(0.1))
	plt.gca().yaxis.set_minor_locator(MultipleLocator(0.1))
	
	plt.tick_params(which='major', width=5, length=25, labelsize=70)
	plt.tick_params(which='minor', width=3, length=15)



	plt.xlabel("Normalized Zero Difference", labelpad=50)
	plt.ylabel("Normalized Frequency", labelpad=50)

	plt.gcf().set_size_inches(30, 24, forward=1)


	plt.savefig("plots/qm.pdf")
	plt.clf()
	
	return normalized_differences

def write_eigenvalues(filename, eigenvalues):
	f = open(filename, "w")
	for eigenvalue in eigenvalues:
		f.write(str(eigenvalue) + "\n")

	f.close()



def write_min_eigenvalue_diff_vs_N():

	n_range = range(1000, 2000)

	f = open("data/min_difference.dat", "a")

	for N in n_range:
		minimum = min(get_eigenvalue_differences(N))
		f.write("{},{}\n".format(N, minimum))

	f.close()





def read_min_eigenvalues_differences_vs_N():
	f = open("data/min_differences.dat", "r")

	lines = f.read().split("\n")

	N_s, mins = [], []
	for line in lines:
		if len(line) != 0:
			content = line.split(",")
			N_s.append( int( content[0] ) )
			mins.append( float( content[1] ) )

	
	return N_s, mins


def plot_min_eigenvalues_differences_vs_N():

	max_N = 200

	N_s, mins = read_min_eigenvalues_differences_vs_N()

	N_s, mins = N_s[:max_N], mins[:max_N]

	plt.plot(N_s, mins, color="orange", lw=5)
	# plt.hist(N_s, weights=mins, bins=100, color="purple", histtype='step', lw=5, normed=True)

	# plt.xscale('log')

	# plt.autoscale()

	plt.xlabel("Matrix Size, $N$", labelpad=30, fontsize=70)
	plt.ylabel("Min. Eigenvalue Difference", labelpad=30, fontsize=70)

	plt.xlim(0, max_N)
	plt.ylim(0, plt.ylim()[1])

	plt.gca().xaxis.set_minor_locator(MultipleLocator(10))
	plt.gca().yaxis.set_minor_locator(MultipleLocator(0.02))
	
	plt.tick_params(which='major', width=5, length=25, labelsize=70)
	plt.tick_params(which='minor', width=3, length=15)

	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.grid()
	plt.savefig("plots/qm_min_eigenvalues_differences.pdf")
	plt.clf()


def plot_max_eigenvalue_vs_N():
	N_s = []

	max_s = []

	for i in range(5, 100):
		N_s.append(i)
		max_s.append( max(get_eigenvalues(i)) )

	plt.plot(N_s, max_s, lw=5, color="green")

	plt.autoscale()

	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.savefig("plots/max_eigenvalues.pdf")
	plt.clf()


def plot_max_eigenvalue_diff_vs_N():
	N_s = []

	max_s = []

	for i in range(5, 100):
		N_s.append(i)
		max_s.append( max(get_eigenvalue_differences(i)) )

	plt.plot(N_s, max_s, lw=5, color="green")

	plt.autoscale()

	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.savefig("plots/max_eigenvalue_diff.pdf")
	plt.clf()


def plot_min_eigenvalue_vs_N():
	N_s = []

	minimum_eigenvalues = []

	for i in range(5, 101):
		N_s.append(i)
		minimum_eigenvalues.append( min([ l for l in np.abs(get_eigenvalues(i)) if l > 1e-5  ] ) )

	plt.plot(N_s, minimum_eigenvalues, lw=5, color="green")
	
	plt.xlabel("Matrix Size, $N$", labelpad=30, fontsize=70)
	plt.ylabel("Min. Eigenvalue", labelpad=30, fontsize=70)

	plt.autoscale()

	plt.gca().xaxis.set_minor_locator(MultipleLocator(10))
	plt.gca().yaxis.set_minor_locator(MultipleLocator(0.05))
	
	plt.tick_params(which='major', width=5, length=25, labelsize=70)
	plt.tick_params(which='minor', width=3, length=15)

	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.savefig("plots/qm_min_eigenvalues.pdf")
	plt.clf()






def plot_qm_eigenvalues():



	plt.hist(get_eigenvalues(100), label="N = 100", bins=20, color="red", edgecolor='red', histtype='step', lw=5, normed=1)
	plt.hist(get_eigenvalues(200), label="N = 200", bins=50, color="blue", edgecolor='blue', histtype='step', lw=5, normed=1)
	plt.hist(get_eigenvalues(500), label="N = 500", bins=50, color="green", edgecolor='green', histtype='step', lw=5, normed=1)
	

	# plt.xscale('log')

	# plt.autoscale()

	plt.xlabel("Eigenvalues", labelpad=30, fontsize=70)
	
	plt.legend()

	# plt.xlim(0, max_N)
	plt.ylim(0, plt.ylim()[1] * 1.2)

	plt.gca().xaxis.set_minor_locator(MultipleLocator(100))
	plt.gca().yaxis.set_minor_locator(MultipleLocator(0.00025))
	
	plt.tick_params(which='major', width=5, length=25, labelsize=70)
	plt.tick_params(which='minor', width=3, length=15)

	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.grid()
	plt.savefig("plots/qm_eigenvalues.pdf")
	plt.clf()





if __name__ == '__main__':
	# eigenvalues = get_eigenvalue_differences(n=1000)
	# write_eigenvalues("data/qm.dat", eigenvalues)

	# plot_
	
	# write_min_eigenvalue_diff_vs_N()

	# plot_min_vs_N()

	# plot_max_eigenvalue_vs_N()

	# plot_max_eigenvalue_diff_vs_N()
	# plot_min_eigenvalue_vs_N()

	# plot_min_eigenvalues_differences_vs_N()

	# plot_qm_eigenvalues()

	pass