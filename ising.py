#!/usr/bin/env python

"""
Python implementation of common model fitting operations to
analyse protein folding data. Simply automates some fitting 
and value calculation. Will be extended to include phi-value 
analysis and other common calculations.

Allows for quick model evaluation and plotting.

Also tried to make this somewhat abstract and modular to
enable more interesting calculations, such as Ising models
and such.

Requirements (recommended python 2.7+):
	- numpy
	- scipy
	- matplotlib

Lowe, A.R. 2015-2016

"""


import inspect
import numpy as np 
import scipy as sp

import matplotlib.pyplot as plt

# use a global optimisation algorithm for the Ising model fitting
from scipy.optimize import differential_evolution

# import PyFolding specific libraries
import constants
import models
import folding

__author__ = "Alan R. Lowe"
__email__ = "a.lowe@ucl.ac.uk"




def kappa(x, DeltaG, m_value):
	""" kappa: energy term for intrinsic stability of a domain """
	return np.exp(-(DeltaG - m_value*x) / constants.RT)

def tau(x, DeltaG, m_value):
	""" tau: energy term for interface stability between domains """
	return np.exp(-(DeltaG - m_value*x) / constants.RT)










class IsingDomain(object):
	""" Template class for derivative Ising domains to describe a
	protein topology.

	Protein topologies can be described using, for example, a list
	of instantiated classes to represent each of the domains.

	For example:

	p = [RepeatDomain, RepeatDomain, LoopDomain, RepeatDomain]

	Variable sharing can be performed by class re-use.

	Args:
		k_func
		t_func

	Member functions:


	Notes:

	"""
	def __init__(self, k_func=kappa, t_func=tau):
		if hasattr(k_func, '__call__'):
			self.__kappa = k_func
		else:
			raise TypeError("Ising domain kappa function must be a callable function")

		if hasattr(t_func, '__call__'):
			self.__tau = t_func
		else:
			raise TypeError("Ising domain tau function must be a callable function")

		self.DG_intrinsic = None
		self.DG_interface = None
		self.m_intrinsic = None
		self.m_interface = None

	def kappa(self, x):
		""" Return the DG_intrinsic """
		return self.__kappa(x)
 
	def tau(self, x):
		""" Return the DG_intrinsic """
		return self.__tau(x)

	def q_i(self, x, folded=1):
		""" Return the q_i matrix """
		return np.matrix([[self.kappa(x)*self.tau(x), folded],[self.kappa(x), folded]])




class RepeatDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self, k_func=kappa, t_func=tau)
		
		# bounds (DG_i, DG_ij, m_i, m_ij)
		self.bounds = ((2., 3.),(-5.,-4.),(-.8,-.5),(-.8,-.5))


class LoopDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self, k_func=kappa, t_func=tau)

		# bounds (DG_i, DG_ij, m_i, m_ij)
		self.bounds = ((2., 3.),(-5.,-4.),(-.8,-.5),(-.8,-.5))




# class GlobalFitIsing(object):
# 	"""
# 	http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html#scipy.optimize.differential_evolution
# 	"""
# 	pass

# 	def __init__(self):
# 		self.proteins = []


# 	def __call__(self, x, *args):
# 		"""This is the call used by the fitting function to evalutate the parameters passed in """

# 		if len(self.proteins) < 2:
# 			raise Error('GlobalFitIsing must have at least two proteins to fit.')

# 		res = []
# 		for protein in self.proteins:
# 			den = protein['curve'].x
# 			y = protein['curve'].y


# 			# NEED TO SET THE PARTITION FUNCTION PARAMS HERE
# 			protein['partition'].DG_intrinsic = x[0]
# 			protein['partition'].DG_interface = x[1]
# 			protein['partition'].m_intrinsic = x[2]
# 			protein['partition'].m_interface = x[3]

# 			# get the fit for this particular fit function
# 			fit = protein['partition'].theta(den)
# 			res.append( np.sum( np.abs(fit-y)**2 ) )

# 		return np.sum(res)


# 	def append(self, curve, n):
# 		""" Append an equilibrium curve, with associated data. """
# 		if isinstance(curve, folding.EquilibriumDenaturationCurve):
# 			# generate a new partition function for the data
# 			q = IsingPartitionFunction(n)
# 			self.proteins.append({'n':n, 'curve':curve, 'partition':q})
# 		else:
# 			raise TypeError('GlobalFitIsing.append must receive an EquilibriumDenaturationCurve as input')



class GlobalFitIsing(object):
	"""
	http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html#scipy.optimize.differential_evolution
	"""
	pass

	def __init__(self):
		self.proteins = []
		self.domains = []


	def __call__(self, x, *args):
		"""This is the call used by the fitting function to evalutate the parameters passed in """

		if len(self.proteins) < 2:
			raise Error('GlobalFitIsing must have at least two proteins to fit.')


		# here we can set all of the parameters for the fit
		for domain in self.domains:
			i = self.domains.index(domain)*4
			domain.DG_intrinsic = x[i]
			domain.DG_interface = x[i+1]
			domain.m_interface = x[i+2]
			domain.m_interface = x[i+3]

		# calculate the sum deviation of the fits from the raw data
		res = []
		for protein in self.proteins:
			den = protein['curve'].x
			y = protein['curve'].y

			# get the fit for this particular fit function
			fit = protein['partition'].theta(den)
			res.append( np.sum( np.abs(fit-y)**2 ) )

		return np.sum(res)


	@property 
	def bounds(self):
		""" Return the expected bounds for fitting based on the domain topologies
		used.
		"""

		for domain in domains:
			bounds = bounds + domain.bounds

		return bounds

	
	def append(self, curve, topology):
		""" Append an equilibrium curve, with associated data. Generate a topology
		and the associated partition function 
		"""
		if isinstance(curve, folding.EquilibriumDenaturationCurve):
			
			# now go through the domain topology and set up domains as necessary
			q_topology = []
			for domain in topology:
				
				if not isinstance(topology, IsingDomain):
					raise TypeError('GlobalFitIsing.append protein topologies must be specified using IsingModel classes')

				if domain not in self.domains:
					new_domain = domain() # instantiate this and append it to the list of domains
					self.domains.append( new_domain )
				else:
					new_domain = self.domains[self.domains.index(domain)]

				q_topology.append(new_domain)


			# generate a new partition function for the data
			q = IsingPartitionFunction( topology = q_topology )

			# append the new data and partition function to the list
			self.proteins.append({'n':n, 'curve':curve, 'partition':q})

		else:
			raise TypeError('GlobalFitIsing.append must receive an EquilibriumDenaturationCurve as input')


# class IsingPartitionFunction(object):
# 	""" General partition function object for Ising models.
# 	"""
# 	def __init__(self, n):
# 		if isinstance(n, int) and n>0:
# 			self.n = n
# 		else:
# 			raise TypeError('IsingPartitionFunction must be instantiated with an INTEGER (n>0) number of vertices')
# 		self.m_intrinsic = 0.
# 		self.m_interface = 0.
# 		self.__DG_intrinsic = np.ones((n,))
# 		self.__DG_interface = np.ones((n,))

# 	def __len__(self):
# 		return self.n

# 	@property 
# 	def DG_intrinsic(self):
# 		return self.__DG_intrinsic
# 	@DG_intrinsic.setter
# 	def DG_intrinsic(self, DeltaG):
# 		if isinstance(DeltaG, list):
# 			if len(DeltaG) == self.n:
# 				self.__DG_intrinsic = np.array(DeltaG)
# 			else:
# 				raise Error("DG_intrinsic array must contain the same number of elements as n")
# 		else:
# 			self.__DG_intrinsic = np.ones((self.n,)) * DeltaG

# 	@property 
# 	def DG_interface(self):
# 		return self.__DG_interface
# 	@DG_interface.setter
# 	def DG_interface(self, DeltaG):
# 		if isinstance(DeltaG, list):
# 			if len(DeltaG) == self.n:
# 				self.__DG_interface = np.array(DeltaG)
# 			else:
# 				raise Error("DG_interface array must contain the same number of elements as n")
# 		else:
# 			self.__DG_interface = np.ones((self.n,)) * DeltaG


# 	def partition(self, x, folded=np.array([])):
# 		""" Populate the partition function given the weights.
# 		"""

# 		if not folded.size:
# 			folded = np.ones((self.n,))

# 		# calculate the full partition function for the fully folded protein
# 		q_i = np.matrix([0,1])
# 		for i in xrange(self.n):
# 			k = kappa(x, self.__DG_intrinsic[i], self.m_intrinsic) 
# 			t = tau(x, self.__DG_interface[i], self.m_interface)
# 			q_i = q_i * np.matrix([[k*t, folded[i]],[k, folded[i]]])
# 		q_i = q_i * np.matrix([1,1]).T
# 		return  np.asarray( q_i.flatten() )[0][0]

# 	def subpartition(self, x, i):
# 		""" Make a subpartition function, setting one to be zero """
# 		folded = np.ones((self.n,))
# 		folded[i] = 0.
# 		return self.partition(x, folded)

# 	def theta(self, x):
# 		""" Return the fraction folded for some x value """
# 		q_n = self.partition(x)
# 		sum_q_i = np.sum( [ self.subpartition(x, i) for i in xrange(self.n)], axis=0 )
# 		theta = sum_q_i / (self.n * q_n)
# 		return 1.-theta

# 	def subpopulation(self, x, i):
# 		""" Return the fraction folded for a sub population of the states """
# 		pass





class IsingPartitionFunction(object):
	""" General partition function object for Ising models.
	"""
	def __init__(self, topology=None):
		if isinstance(topology, list):

			for domain in topology:
				if not isinstance(domain, IsingDomain):
					raise TypeError('IsingPartitionFunction: Model topologies must be specified using IsingDomain classes')

		else:
			raise TypeError('IsingPartitionFunction: model topology must be specified')


		# if everything is OK, set up the topology for this partition function		
		self.init_topology( topology )

	def init_topology(self, topology=None):
		self.topology = topology
		self.n = len(self.topology)

	def __len__(self):
		return self.n

	def partition(self, x, folded=np.array([])):
		""" Populate the partition function given the weights.
		"""

		if not folded.size:
			folded = np.ones((self.n,))

		# calculate the full partition function for the fully folded protein
		q_i = np.matrix([0,1])
		for domain in self.topology:
			q_i = q_i * domain.q_i(x, folded=folded[self.topology.index(domain)])
		q_i = q_i * np.matrix([1,1]).T
		return  np.asarray( q_i.flatten() )[0][0]

	def subpartition(self, x, i):
		""" Make a subpartition function, setting one to be zero """
		folded = np.ones((self.n,))
		folded[i] = 0.
		return self.partition(x, folded)

	def theta(self, x):
		""" Return the fraction folded for some x value """
		q_n = self.partition(x)
		sum_q_i = np.sum( [ self.subpartition(x, i) for i in xrange(self.n)], axis=0 )
		theta = sum_q_i / (self.n * q_n)
		return 1.-theta

	def subpopulation(self, x, i):
		""" Return the fraction folded for a sub population of the states """
		pass





class HomozipperIsingEquilibrium(models.FitModel):
	""" Homopolymer Zipper Ising model

	Notes:
		THIS IS THE NON-MATRIX, EXACT SOLUTION

		Aksel and Barrick. Analysis of repeat-protein folding using 
		nearest-neighbor statistical mechanical models. 
		Methods in enzymology (2009) vol. 455 pp. 95-125
	"""
	def __init__(self):
		models.FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([7, 0.2, -.53, -4.6, -0.6])
		self.constants = (('n',7),)

	def fit_func(self, x, n, DG_intrinsic, m_intrinsic, DG_interface, m_interface):
		
		# clamp to prevent instability
		if DG_intrinsic<0. or DG_interface>0.:
			return models.FIT_ERROR(x)

		k = kappa(x, DG_intrinsic, m_intrinsic) 
		t = tau(x, DG_interface, m_interface)
		pre_factor = k/(n*(k*t-1))
		numerator = n*(k*t)**(n+2) - (n+2)*(k*t)**(n+1) + (n+2)*k*t-n
		denominator = (k*t-1)**2 + k*((k*t)**(n+1) - (n+1)*k*t+n )
		theta = pre_factor * (numerator / denominator)
		return theta





def fit_heteropolymer_ising(equilibrium_curves=[], num_repeats=[]):
	""" An example script to fit a series of data sets to a heteropolymer ising model.
	"""

	fit_func = GlobalFitIsing()

	# set up the global fit
	print 'Appending {0:d} curves to GlobalFitIsing...'.format(len(equilibrium_curves))
	for protein, n in zip(equilibrium_curves, num_repeats):
		fit_func.append(protein, n)

	# bounds (DG_i, DG_ij, m_i, m_ij)
	bounds = ((2., 3.),(-5.,-4.),(-.8,-.5),(-.8,-.5))

	# do the fitting
	print 'Performing global optimisation of Ising model...'
	result = differential_evolution(fit_func, bounds, disp=True, popsize=10, tol=1e-1)
	print 'Done.'

	print result

	plot_Ising(fit_func)

	

	return result



def plot_Ising(fit_func):
	""" Function to plot fitted Ising model data.
	"""

	# plot the fits
	plt.figure()
	plt.subplot(1,2,1)
	for protein in fit_func.proteins:
		xv = protein['curve'].x
		plt.plot(xv, protein['curve'].y, 's:')
		
	plt.legend( [p['curve'].ID for p in fit_func.proteins] )

	for protein in fit_func.proteins:
		xv = protein['curve'].x
		plt.plot(xv, protein['partition'].theta(xv), 'k-')

	plt.xlabel(fit_func.proteins[0]['curve'].denaturant_label+" (M)")
	plt.ylabel('Fraction folded')


	# now plot the first derivative
	pk_max = []
	plt.subplot(1,2,2)
	for protein in fit_func.proteins:
		xv = np.linspace(0.,10.,1000)
		first_deriv = np.gradient(protein['partition'].theta(xv))
		pk_max.append( (protein['n'], xv[np.argmax(np.abs(first_deriv))]) )
		plt.plot(xv, np.abs(first_deriv), '-')
	plt.legend( [p['curve'].ID for p in fit_func.proteins] )
	plt.show()

	plt.figure()
	x,y = zip(*pk_max)
	plt.plot(x,y,'k-')
	plt.show()





def Test_EWAN():
	""" Test with Ewan's data 

	TODO:
		Predictions of other curves with fit results
		Subpopulations
		Error analysis? RMSD of curves?

	"""

	datadir = '/Users/ubcg83a/Dropbox/AlanLoweCollaboration (1)/IsingFits&Numbers/HCTPR_Datasets/'

	#files = ['H2', 'H3', 'H4', 'H5', 'H6', 'H8', 'H10']
	#helices = [5,7,9,11,13,17,21]

	files = ['H2', 'H3', 'H4', 'H5']
	helices = [5,7,9,11]

	curves = [ folding.read_equilibrium_data(datadir, f+'.csv') for f in files]

	fit_heteropolymer_ising(curves, helices)

	pass



def test_domains():
	protein = [RepeatDomain, RepeatDomain, LoopDomain, RepeatDomain]

	for domain in protein:
		print domain

	print set(protein)




def main():


	# x = np.linspace(-10., 10., 100)
	# plt.figure()
	# l = ['n='+str(i) for i in xrange(1,11)]
	# for n in xrange(1,11):
	# 	q = IsingPartitionFunction(n)
	# 	q.DG_intrinsic = 0.2
	# 	q.DG_interface = -4.6
	# 	q_n = q.theta(x)
	# 	plt.plot(x,q_n,'.-')
	# plt.legend(l)
	# plt.show()


	# # get the homozipper for comparison
	# h = HomozipperIsingEquilibrium()
	# p = h.default_params.tolist()
	# p[0] = n

	# q_h = h.fit_func(x, *p)

	# plt.figure()
	# plt.subplot(1,2,1)
	# plt.plot(x, q_n, 'bx-')
	# plt.plot(x, q_h, 'r+-')
	# plt.subplot(1,2,2)
	# plt.plot(x, q_n - q_h, 'ko-')
	# plt.show()

	# make some fake data
	x = np.linspace(0., 10., 100)

	equilibrium_curves = []
	num_repeats = []


	for n in xrange(3,30,3):
	 	q = IsingPartitionFunction(n)
	 	q.m_intrinsic = -0.53
		q.m_interface = -0.6
	 	q.DG_intrinsic = 0.0
	 	q.DG_interface = -4.6
	 	q_n = q.theta(x) + np.random.randn(len(x),)*0.01


	 	# make a fake curve
	 	eq = folding.EquilibriumDenaturationCurve(ID='Protein (n={0:d})'.format(n))
	 	eq.denaturant_label = 'Urea (M)'
		eq.curves = ['Fake']
		eq.denaturant = {'Fake': x}
		eq.signal = {'Fake': q_n}

	 	num_repeats.append(n)
	 	equilibrium_curves.append(eq)

	
	fit_heteropolymer_ising(equilibrium_curves, num_repeats)


if __name__ == "__main__":
	#Test_EWAN()
	test_domains()