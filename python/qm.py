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


def construct_operators(n):
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

	eigenvalues = sorted( list(set(np.absolute(np.array(sorted(np.linalg.eigvalsh(berry_H)))))) )
	# print eigenvalues


	normalized_differences = np.diff(eigenvalues)
	normalized_differences *= 1 / np.mean(normalized_differences)

	filtered_differences = [f for f in normalized_differences if f > 1e-6]

	normalized_differences = filtered_differences


	plt.hist(normalized_differences, color="red", bins=100, lw=5, edgecolor="red", normed=1)

	plt.autoscale()
	plt.xlim(0, 3)

	plt.xlabel("Normalized Zero Difference", labelpad=50)
	plt.ylabel("Normalized Frequency", labelpad=50)

	plt.gcf().set_size_inches(30, 24, forward=1)


	plt.savefig("plots/qm.pdf")
	plt.clf()
	


if __name__ == '__main__':
	construct_operators(n=10000)
