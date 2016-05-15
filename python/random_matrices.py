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

from scipy.stats import binned_statistic



mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)



def generate_GOE_matrix(matrix_size):
	

	random_matrix = np.random.randn(matrix_size, matrix_size)
	transpose = np.transpose(random_matrix)

	return ( random_matrix + transpose ) / 2


def generate_GUE_matrix(matrix_size):
	

	random_matrix_real = np.matrix( np.random.randn(matrix_size, matrix_size) )
	random_matrix_im = np.matrix( 1.j * np.random.randn(matrix_size, matrix_size) )

	random_matrix = random_matrix_real + random_matrix_im

	random_matrix_H = random_matrix.getH()


	return ( random_matrix + random_matrix_H ) / 2



def generate_matrix_ensemble(matrix_size, number_of_matrices, ensemble_type):
	random_matrices = []

	for i in range(number_of_matrices):
		
		if ensemble_type == "GOE":
			random_matrix = generate_GOE_matrix(matrix_size)
		elif ensemble_type == "GUE":
			random_matrix = generate_GUE_matrix(matrix_size)

		random_matrices.append( random_matrix )

	return random_matrices




def compute_eigenvalues(matrices):
	all_eigenvalues = []
	for matrix in matrices:
		eigenvalues = sorted( np.linalg.eigvalsh(matrix) )
		all_eigenvalues.append( eigenvalues )

	return all_eigenvalues



def calculate_eigenvalues_differences(matrices):
	
	matrix_size = np.sqrt( matrices[0].size )

	
	all_eigenvalues = compute_eigenvalues(matrices)

	# flatten_eigenvalues = sorted( np.ndarray.flatten( np.array( all_eigenvalues ) ) )
	# flatten_eigenvalues = sorted([eigenvalues[int(matrix_size / 2)] for eigenvalues in all_eigenvalues] )
	# eigenvalues_differences = ( np.diff(flatten_eigenvalues) )
	eigenvalues_differences = [ eigenvalues[int(matrix_size / 2) + 1] - eigenvalues[int(matrix_size / 2)] for eigenvalues in all_eigenvalues ]

	# expanded_weights = np.ndarray.flatten( np.array( [ [np.real(weight)]*matrix_size for weight in weights ] ) )



	# normalized_eigenvalues_differences = eigenvalues_differences
	normalized_eigenvalues_differences = eigenvalues_differences / np.mean(eigenvalues_differences)



	return normalized_eigenvalues_differences
	

def plot_normalized_differences(matrix_size, number_of_matrices, ensemble_type="GOE"):

	random_matrices = generate_matrix_ensemble(matrix_size, number_of_matrices, ensemble_type)


	differences = calculate_eigenvalues_differences(random_matrices)


	plt.hist(differences, color="red", lw=5, edgecolor="red", bins=500, normed=1)

	plt.autoscale()

	plt.xlim(0, 3)
	
	plt.gcf().set_size_inches(30, 24, forward=1)


	plt.savefig("plots/" + str(ensemble_type) + "_eigenvalues_differences.pdf")
	plt.clf()



def get_eigenvalues_differences(matrix_size, number_of_matrices):

	GOE_differences, GUE_differences = [], []

	for i in range(number_of_matrices):
		random_GOE_matrix = generate_GOE_matrix(matrix_size)
		random_GUE_matrix = generate_GUE_matrix(matrix_size)

		# Calculate eigenvalues.
		GOE_eigenvalues = sorted( np.linalg.eigvalsh(random_GOE_matrix) )
		GUE_eigenvalues = sorted( np.linalg.eigvalsh(random_GUE_matrix) )

		GOE_differences.append( GOE_eigenvalues[int(matrix_size / 2) + 1] - GOE_eigenvalues[int(matrix_size / 2)] )
		GUE_differences.append( GUE_eigenvalues[int(matrix_size / 2) + 1] - GUE_eigenvalues[int(matrix_size / 2)] )

	GOE_normed_differences = GOE_differences / np.mean(GOE_differences)
	GUE_normed_differences = GUE_differences / np.mean(GUE_differences)

	return [GOE_normed_differences, GUE_normed_differences]	





def plot_eigenvalues_differences(matrix_size, number_of_matrices):

	
	GOE_differences, GUE_differences = get_eigenvalues_differences(matrix_size, number_of_matrices)

	plt.hist(GOE_differences, color="red", lw=5, histtype='step', edgecolor="red", bins=50, normed=1, label="Gaussian Orthogonal Ensemble")
	plt.hist(GUE_differences, color="blue", lw=5, histtype='step', edgecolor="blue", bins=50, normed=1, label="Gaussian Unitary Ensemble")


	plt.legend()


	plt.autoscale()

	plt.xlim(0, 3)
	
	plt.gcf().set_size_inches(30, 24, forward=1)


	plt.savefig("plots/eigenvalues_differences.pdf")
	plt.clf()
	



if __name__ == '__main__':
	# plot_normalized_differences(N=50, number_of_matrices=10000, ensemble_type="GOE")
	plot_eigenvalues_differences(matrix_size=20, number_of_matrices=100000)
