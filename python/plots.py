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


# RootPy
from rootpy.plotting import Hist, HistStack, Legend
import rootpy.plotting.root2matplotlib as rplt
from rootpy.plotting import Hist2D


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

import rootpy.plotting.views


from random_matrices import *



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
	zeros = parse_file("data/zeros_external.dat")

	plt.hist(zeros, color="red", histtype='step', bins=25, lw=5, normed=1)


	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.tight_layout(pad=1.08, h_pad=1.08, w_pad=1.08)


	plt.savefig("plots/zeros_external.pdf")
	plt.clf()


def plot_normalized_differences():
	zeros = np.array( parse_file("data/zeros_2M.dat") )

	normalized_differences = np.diff(zeros) * np.log(max(zeros) / (2 * np.pi)) / (2 * np.pi)


	plt.hist(normalized_differences, color="red", bins=500, lw=5, edgecolor="red", normed=1)

	plt.autoscale()
	plt.xlim(0, 3)

	plt.xlabel("Normalized Zero Difference", labelpad=50)
	plt.ylabel("Normalized Frequency", labelpad=50)

	plt.gcf().set_size_inches(30, 24, forward=1)


	plt.savefig("plots/normalized_differences.png")
	plt.clf()

def plot_zeros_and_eigenvalue_differences():
	eigenvalue_differences = calculate_eigenvalues_differences(N=5, number_of_matrices=1000000)

	print np.mean(eigenvalue_differences)

	zeros = np.array( parse_file("data/zeros_2M.dat") )

	normalized_zeros_differences = np.diff(zeros) * np.log(max(zeros) / (2 * np.pi)) / (2 * np.pi)


	plt.hist(normalized_zeros_differences, color="red", label="Zeros Differences", bins=500, histtype='step', lw=5, edgecolor="red", normed=1)

	plt.hist(eigenvalue_differences, color="blue", label="Eigenvalues Differences", bins=500, histtype='step', lw=5, edgecolor="blue", normed=1)

	plt.legend()

	plt.autoscale()

	plt.xlim(0, 3)
	
	plt.gcf().set_size_inches(30, 24, forward=1)

	plt.grid()
	plt.savefig("plots/everything.pdf")
	plt.clf()


# plot_zeros()

plot_zeros_and_eigenvalue_differences()