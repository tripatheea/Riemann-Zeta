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



mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)




def generate_random_matrix(N):
	

	random_matrix = np.random.randn(N, N)
	transpose = np.transpose(random_matrix)


	return random_matrix + transpose

	

def generate_random_matrices(N, number_of_matrices):
	random_matrices = []
	for i in range(number_of_matrices):
		random_matrices.append( generate_random_matrix(N) )

	return random_matrices


def probability_distribution(H):
	N = np.shape(H)[0]
	C_N = ((2 * np.pi)**( - (N * (N + 1)) / 2 )) * N**( (N**2) / 2 )
	trace_H2 = np.trace( H*H )

	weight = C_N * np.exp( - (N / 2) * trace_H2 )

	return weight

def weight_matrices(matrices):
	weights = []
	for matrix in matrices:
		weight = probability_distribution(matrix)
		weights.append(weight)

	return weights

def compute_eigenvalues(matrices):
	all_eigenvalues = []
	for matrix in matrices:
		eigenvalues = sorted( np.linalg.eigvalsh(matrix) )
		all_eigenvalues.append( eigenvalues )

	return all_eigenvalues



def calculate_eigenvalues_differences(N, number_of_matrices):
	random_matrices = generate_random_matrices(N, number_of_matrices)
	# weights = weight_matrices(random_matrices)
	all_eigenvalues = compute_eigenvalues(random_matrices)

	# flatten_eigenvalues = sorted( np.ndarray.flatten( np.array( all_eigenvalues ) ) )
	# flatten_eigenvalues = sorted([eigenvalues[int(N / 2)] for eigenvalues in all_eigenvalues] )
	# eigenvalues_differences = ( np.diff(flatten_eigenvalues) )
	eigenvalues_differences = [ eigenvalues[int(N / 2) + 1] - eigenvalues[int(N / 2)] for eigenvalues in all_eigenvalues ]

	# expanded_weights = np.ndarray.flatten( np.array( [ [np.real(weight)]*N for weight in weights ] ) )



	# normalized_eigenvalues_differences = eigenvalues_differences
	normalized_eigenvalues_differences = eigenvalues_differences / np.mean(eigenvalues_differences)



	return normalized_eigenvalues_differences
	


def plot_normalized_differences(N, number_of_matrices):
	differences = calculate_eigenvalues_differences(N, number_of_matrices)


	plt.hist(differences, color="red", lw=5, edgecolor="red", bins=500, normed=1)

	plt.autoscale()

	plt.xlim(0, 3)
	
	plt.gcf().set_size_inches(30, 24, forward=1)


	plt.savefig("plots/eigenvalues_differences.pdf")
	plt.clf()


if __name__ == '__main__':
	plot_normalized_differences(N=50, number_of_matrices=100000)
