from __future__ import division

from subprocess import call

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

from scipy.special import zeta

from scipy.stats import binned_statistic

import euler
import random_matrices
import qm
import dirichlet
import riemann_siegel


mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)

def parse_file(f):
	f = open(f, 'r')
	lines = f.read().split("\n")
	zeros = [float(line) for line in lines if len(line) != 0]
	return zeros

def plot_zeros():
	zeros = parse_file("data/zeros_2M.dat")

	plt.hist(zeros, color="red", histtype='step', bins=25, lw=5, normed=1)

	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)


	plt.savefig("plots/zeta_zeros.pdf")
	plt.clf()


def plot_normalized_differences():
	zeros = np.array( parse_file("data/zeros_2M.dat") )

	normalized_differences = np.diff(zeros) * np.log(max(zeros) / (2 * np.pi)) / (2 * np.pi)


	plt.hist(normalized_differences, label="Successive Difference of the Roots of $\zeta(t)$", color="red", bins=500, lw=8, histtype='step', edgecolor="red", normed=1)

	plt.legend()

	plt.autoscale()
	plt.xlim(0, 3)
	plt.ylim(0, 1.1)

	plt.xlabel("Successive Difference", labelpad=30, fontsize=70)
	plt.ylabel("Normalized Frequency", labelpad=30, fontsize=70)

	plt.gca().xaxis.set_minor_locator(MultipleLocator(0.1))
	plt.gca().yaxis.set_minor_locator(MultipleLocator(0.02))
	
	plt.tick_params(which='major', width=5, length=25, labelsize=70)
	plt.tick_params(which='minor', width=3, length=15)


	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.grid()
	plt.savefig("plots/zeta_zeros_normalized_differences.pdf")
	plt.clf()


def plot_random_matrices_eigenvalues_differences():
	GOE_differences, GUE_differences = random_matrices.get_eigenvalues_differences(25, 100000)

	plt.hist(GOE_differences, color="blue", label="Gaussian Orthogonal Ensemble", bins=100, histtype='step', lw=5, edgecolor="blue", normed=1)
	plt.hist(GUE_differences, color="orange", label="Gaussian Unitrary Ensemble", bins=100, histtype='step', lw=5, edgecolor="orange", normed=1)

	plt.legend()

	plt.autoscale()

	plt.xlabel("Successive Difference", labelpad=30, fontsize=70)
	plt.ylabel("Normalized Frequency", labelpad=30, fontsize=70)

	plt.xlim(0, 3)
	plt.ylim(0, 1.2*plt.gca().get_ylim()[1])
	
	plt.gca().xaxis.set_minor_locator(MultipleLocator(0.1))
	plt.gca().yaxis.set_minor_locator(MultipleLocator(0.04))
	
	plt.tick_params(which='major', width=5, length=25, labelsize=70)
	plt.tick_params(which='minor', width=3, length=15)

	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.grid()
	plt.savefig("plots/random_matrices_eigenvalues_differences.pdf")
	plt.clf()


def plot_GUE_eigenvalues_and_zeros_differences():
	GOE_differences, GUE_differences = random_matrices.get_eigenvalues_differences(25, 100)

	zeros = np.array( parse_file("data/zeros_2M.dat") )
	normalized_zeros_differences = np.diff(zeros) * np.log(max(zeros) / (2 * np.pi)) / (2 * np.pi)
	normalized_zeros_differences = normalized_zeros_differences / np.mean(normalized_zeros_differences)
	
	print min(normalized_zeros_differences)

	'''
	zeta_bin_count = int( ( max(normalized_zeros_differences) - min(normalized_zeros_differences) ) / 0.1 )
	GUE_bin_count = int( ( max(GUE_differences) - min(GUE_differences) ) / 0.1 )


	plt.hist(normalized_zeros_differences, color="blue", label="Non-trivial Roots of $\zeta(t)$", bins=zeta_bin_count, histtype='step', lw=5, edgecolor="green", normed=1)
	plt.hist(GUE_differences, color="orange", label="Gaussian Unitrary Ensemble", bins=GUE_bin_count, histtype='step', lw=5, edgecolor="orange", normed=1)

	plt.legend()

	plt.autoscale()

	plt.xlabel("Successive Difference", labelpad=30, fontsize=70)
	plt.ylabel("Normalized Frequency", labelpad=30, fontsize=70)

	plt.xlim(0, 3)
	plt.ylim(0, 1.2*plt.gca().get_ylim()[1])
	
	plt.gca().xaxis.set_minor_locator(MultipleLocator(0.1))
	plt.gca().yaxis.set_minor_locator(MultipleLocator(0.04))
	
	plt.tick_params(which='major', width=5, length=25, labelsize=70)
	plt.tick_params(which='minor', width=3, length=15)

	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.grid()
	plt.savefig("plots/GUE_eigenvalues_and_zeros_differences.pdf")
	plt.clf()

	'''




def plot_zeros_and_eigenvalue_differences():
	GOE_differences, GUE_differences = random_matrices.get_eigenvalues_differences(25, 100000)

	zeros = np.array( parse_file("data/zeros_2M.dat") )

	qm_eigenvalues = np.array( parse_file("data/qm.dat") )

	normalized_zeros_differences = np.diff(zeros) * np.log(max(zeros) / (2 * np.pi)) / (2 * np.pi)

	normalized_zeros_differences = normalized_zeros_differences / np.mean(normalized_zeros_differences)
	
	zeta_bin_count = int( ( max(normalized_zeros_differences) - min(normalized_zeros_differences) ) / 0.05 )
	GOE_bin_count = int( ( max(GOE_differences) - min(GOE_differences) ) / 0.05 )
	qm_bin_count = int( ( max(qm_eigenvalues) - min(qm_eigenvalues) ) / 0.05 )

	plt.hist(normalized_zeros_differences, color="red", label="Non-trivial Roots of $\zeta(t)$", bins=zeta_bin_count, histtype='step', lw=5, edgecolor="red", normed=1)
	plt.hist(GOE_differences, color="blue", label="Gaussian Orthogonal Ensemble", bins=GOE_bin_count, histtype='step', lw=5, edgecolor="blue", normed=1)
	plt.hist(qm_eigenvalues, color="green", label="Berry's Hamiltonian", bins=qm_bin_count, histtype='step', lw=5, edgecolor="green", normed=1)

	plt.legend()

	plt.xlabel("Successive Difference", labelpad=30, fontsize=70)
	plt.ylabel("Normalized Frequency", labelpad=30, fontsize=70)

	plt.autoscale()

	plt.xlim(0, 3)
	plt.ylim(0, 1.2*plt.gca().get_ylim()[1])
	
	plt.gca().xaxis.set_minor_locator(MultipleLocator(0.1))
	plt.gca().yaxis.set_minor_locator(MultipleLocator(0.04))
	
	plt.tick_params(which='major', width=5, length=25, labelsize=70)
	plt.tick_params(which='minor', width=3, length=15)

	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.grid()
	plt.savefig("plots/everything.pdf")
	plt.clf()



# plot_zeros()

# plot_normalized_differences()

# plot_random_matrices_eigenvalues_differences()

# plot_GUE_eigenvalues_and_zeros_differences()

# plot_zeros_and_eigenvalue_differences()