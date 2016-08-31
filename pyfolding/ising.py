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

import time
import inspect
import numpy as np 
import scipy as sp

import matplotlib.pyplot as plt

# use a global optimisation algorithm for the Ising model fitting
from scipy.optimize import differential_evolution, minimize, leastsq

# import PyFolding specific libraries
import constants
import models
import core

__author__ = "Alan R. Lowe"
__email__ = "a.lowe@ucl.ac.uk"







def kappa(x, DeltaG, m_value):
	""" kappa: energy term for intrinsic stability of a domain """
	return np.exp(-(DeltaG - m_value*x) / constants.RT)

def tau(x, DeltaG, m_value):
	""" tau: energy term for interface stability between domains """
	m_value = 0.0
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

		self.DG_intrinsic = 2.
		self.DG_interface = -4.
		self.m_intrinsic = -.5
		self.m_interface = -.5

		self.labels = ['DG_i', 'DG_ij', 'm_i', 'm_ij']

		self.q_func = lambda x, folded: np.matrix([[self.kappa(x)*self.tau(x), folded],[self.kappa(x), folded]])

	def kappa(self, x):
		""" Return the DG_intrinsic """
		return self.__kappa(x, self.DG_intrinsic, self.m_intrinsic)
 
	def tau(self, x):
		""" Return the DG_intrinsic """
		return self.__tau(x, self.DG_interface, self.m_interface)

	def q_i(self, x, folded=1):
		""" Return the q_i matrix """
		return self.q_func(x, folded)
		
	


class RepeatDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self, k_func=kappa, t_func=tau)
		
		# bounds (DG_i, DG_ij, m_i, m_ij)
		self.bounds = ((1., 9.),(-9.,-1.),(-2.5,-.1),(-2.5,-.1))
		self.name = "Repeat"



class HelixDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self, k_func=kappa, t_func=tau)
		
		# bounds (DG_i, DG_ij, m_i, m_ij)
		#self.bounds = ((1., 6.),(-9.,-1.),(-.8,-.5),(-.8,-.5))
		self.bounds = ((1., 9.),(-9.,-1.),(-2.5,-.1),(-2.5,-.1))
		self.name = "Helix"


class LoopDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self, k_func=kappa, t_func=tau)

		# bounds (DG_i, DG_ij, m_i, m_ij)
		#self.bounds = ((1., 6.),(-9.,-1.),(-.8,-.5),(-.8,-.5))
		self.bounds = ((1., 9.),(-9.,-1.),(-2.5,-.1),(-2.5,-.1))
		self.name = "Loop"


class CapDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self, k_func=kappa, t_func=tau)

		# bounds (DG_i, DG_ij, m_i, m_ij)
		#self.bounds = ((1., 6.),(-9.,-1.),(-.8,-.5),(-.8,-.5))
		self.bounds = ((1., 9.),(-9.,-1.),(-2.5,-.1),(-2.5,-.1))
		self.name = "Cap"
		self.q_func = lambda x, folded: np.matrix([[self.kappa(x), folded],[self.kappa(x), folded]])



class GlobalFitWrapper(object):
	""" GlobalFitWrapper

	A wrapper for global fitting of the Ising model. Collects together all of the protein 
	topologies, shares fitting paramers and calls the respective objective functions.

	Args:

	Members:

	Member functions:
		

	Notes:
		2016/07/30 (ARL) - Need to sort out proper variable sharing

	"""

	def __init__(self):
		self.proteins = []
		self.domains = []
		self.domain_types = []


	def __call__(self, x, *args):
		"""This is the call used by the fitting function to evalutate the parameters passed in """


		if len(self.proteins) < 1:
			raise ValueError('GlobalFitIsing must have at least one curve to fit.')

		# here we can set all of the parameters for the fit
		for idx, domain in enumerate(self.domains):
			domain.DG_intrinsic = x[0]	# TODO - proper sharing of params across domain types
			domain.DG_interface = x[4*idx+1]
			domain.m_intrinsic = x[4*idx+2]
			domain.m_interface = x[4*idx+3]


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

		bounds = ()
		for domain in self.domains:
			for b in domain.bounds:
				bounds = bounds + (b,)

		return bounds

	
	def append(self, curve, topology):
		""" Append an equilibrium curve, with associated data. Generate a topology
		and the associated partition function 
		"""
		if not isinstance(curve, core.EquilibriumDenaturationCurve):
			raise TypeError('GlobalFitIsing.append must receive an EquilibriumDenaturationCurve as input')


		# now go through the domain topology and set up domains as necessary
		q_topology = []
		for domain in topology:
			
			#if not issubclass(topology, IsingDomain):
			#	raise TypeError('GlobalFitIsing.append protein topologies must be specified using IsingModel classes')

			if domain not in self.domain_types:
				new_domain = domain() # instantiate this and append it to the list of domains
				self.domains.append( new_domain )
				self.domain_types.append( domain )
			else:
				new_domain = self.domains[self.domain_types.index(domain)]

			q_topology.append(new_domain)


		# generate a new partition function for the data
		q = IsingPartitionFunction( topology = q_topology )

		# append the new data and partition function to the list
		self.proteins.append({'n':len(q_topology), 'curve':curve, 'partition':q})


	@property 
	def domain_params(self):
		""" Return the ordered list of domain parameters for printing the results of the fit """
		domain_params = []
		for domain in self.domains:
			for l in domain.labels:
				domain_params.append(domain.name+" "+l)
		return domain_params



			




class IsingPartitionFunction(object):
	""" General partition function object for Ising models.
	"""
	def __init__(self, topology=None):
		if not isinstance(topology, list):
			raise TypeError('IsingPartitionFunction: model topology must be specified')

		for domain in topology:
			if not isinstance(domain, IsingDomain):
				raise TypeError('IsingPartitionFunction: Model topologies must be specified using IsingDomain classes')


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
		for i in xrange(self.n):
			domain = self.topology[i]
			q_i = q_i * domain.q_i(x, folded=folded[i])
		q_i = q_i * np.matrix([1,1]).T
		return  np.asarray( q_i.flatten() )[0][0]

	def subpartition(self, x, i, rev=False):
		""" Make a subpartition function, setting one to be zero """
		folded = np.ones((self.n,))
		folded[i] = 0.
		if rev: folded = 1.-folded
		return self.partition(x, folded)

	def theta(self, x):
		""" Return the fraction folded for some x value """
		q_n = self.partition(x)
		sum_q_i = np.sum( [ self.subpartition(x, i) for i in xrange(self.n)], axis=0 )
		theta = sum_q_i / (self.n * q_n)
		return 1.-theta

	def subpopulation(self, x, i):
		""" Return the fraction folded for a sub population of the states """
		q_n = self.partition(x)
		q_i = self.subpartition(x, i, rev=True) 
		#sum_q_i = np.sum( [ self.subpartition(x, i, rev=True) for i in xrange(self.n)], axis=0 )
		theta = q_i / q_n
		return 1.-theta





def calculate_error_from_jacobian(jac, x=np.linspace(0,10,100)):
	"""
	Calculate Hessian from jacobian, and covariance from Hessian: 
	covar = (J^T . J)^{-1}.

	Then calculate the error based on the SE of the variance.

	This is definitely a bit wonky at the moment!

	Notes:
		This is **NOT** tested at all

	"""

	num_params = len(np.ravel(jac))

	if np.linalg.det( np.dot(np.matrix(jac).T, np.matrix(jac)) ) == 0.:
		print "Warning: Determinant of zero indicates errors are unlikely to represent true error"
		return [np.inf for p in xrange(num_params)]

	covar = np.linalg.pinv( np.dot(np.matrix(jac).T, np.matrix(jac)) )
	errors = [np.sqrt(covar[p,p]) / np.sqrt(1.*len(x)) for p in xrange(num_params)]
	return errors




class FitProgress(object):
	""" Class to take care of updating the user as to the fitting progress
	"""

	def __init__(self, update_freq=10):
		self.__iter = 0
		self.__val = None
		self.__update_freq = update_freq
		self.__timer = time.time()
		self.__iteration_times = []
		self.__convergence_tzero = None

	def __call__(self, xk, convergence=None):
		
		# keep track of the loop timing
		this_iter = time.time() - self.__timer
		self.__timer = time.time()
		self.__iteration_times.append(this_iter)

		# keep track of the convergence
		if self.__iter == 0:
			self.__convergence_tzero = convergence

		# iterate
		self.__iter+=1

		# are we finished?
		if isinstance(xk, bool):
			if xk: return

		if self.__iter % self.__update_freq == 0:

			# remove some loop times and calculate an average:
			sum_time, loop_iter = 0, 0
			while self.__iteration_times and loop_iter < self.__update_freq:
				sum_time += self.__iteration_times.pop(0)
				loop_iter += 1

			# average time per iteration
			avg_time = sum_time / loop_iter

			# NOTE: this doesn't make a lot of sense since there is a non-linear approach to the result!!
			# rate of convergence
			convergence_rate = np.abs((convergence - self.__convergence_tzero) / (avg_time * self.__iter))

			# make an assumption about time remaining
			time_remaining = (1.-convergence) / convergence_rate

			#print " - Fitting in progress (Iteration: {0:d}, Convergence: {1:.5E}, Timing: {2:2.2f}s, Remaining: {3:2.2f}s) ".format(self.__iter, convergence, avg_time, time_remaining)
			print " - Fitting in progress (Iteration: {0:d}, Convergence: {1:.5E}, Timing: {2:2.2f}s) ".format(self.__iter, convergence, avg_time)



def fit_heteropolymer(equilibrium_curves=[], topologies=[], popsize=10, tol=1e-8):
	""" An example script to fit a series of data sets to a heteropolymer ising model.

	
	Notes:
		Optimisation is performed using differential evolution (a GA)
		http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html#scipy.optimize.differential_evolution
	
	"""

	# do some parsing of the input
	if not isinstance(equilibrium_curves, list):
		raise TypeError('equilibrium_curves must be a list of EquilibriumDenaturationCurve type')

	if not isinstance(topologies, list):
		raise TypeError('topologies must be a list of IsingDomain type')	

	fit_func = GlobalFitWrapper()

	# set up the global fit
	print 'Appending {0:d} curves to GlobalFitIsing...'.format(len(equilibrium_curves))
	for protein, topology in zip(equilibrium_curves, topologies):
		fit_func.append(protein, topology)
		print ' + added {0:s} with topology {1:s}'.format(protein.ID, [d().name for d in topology])

	# do the fitting
	print '\nPerforming global optimisation of Ising model ({0:d} curves, Population size: {1:d}, Tolerance: {2:.2E})...'.format(len(equilibrium_curves), popsize, tol)

	# give the user some feedback if this is going to take some time
	if popsize > 10 or len(equilibrium_curves)>1:
		callback = FitProgress()
	else:
		callback = None

	r = differential_evolution(fit_func, fit_func.bounds, disp=False, popsize=popsize, tol=tol, callback=callback)
	

	if not r.success:
		print "Could not find a solution..."
		return None
	#TODO: calculate the errors on each of the parameters

	if hasattr(r, 'jac'):
		# the fit has the jacobian
		jac = r.jac
	else:
		r_min = minimize(fit_func, np.copy(r.x), method='L-BFGS-B', bounds=fit_func.bounds)
		jac = r_min.jac

	# calculate the errors
	r_err = calculate_error_from_jacobian(jac)

	# make a zipped list of params, fit values and estimated errors
	result = zip(fit_func.domain_params, r.x.tolist(), r_err) 

	print '\nFitting results (NOTE: Careful with the errors here): '
	for res in result:
		print u"{0:s}: {1:2.5f} \u00B1 {2:2.5f} ".format(res[0], res[1], res[2])


	plot_Ising(fit_func)
	plot_folded(fit_func.proteins[-1]['partition'])

	return result


def plot_Ising(fit_func):
	""" Function to plot fitted Ising model data.
	"""


	cmap = ['ro', 'mo', 'go', 'co', 'bo', 'ko', 'rt', 'mv', 'gv', 'cv', 'bv', 'kv']

	# plot the fits
	plt.figure(figsize=(14,8))
	#plt.subplot(1,2,1)
	ax1 = plt.subplot2grid((2,2), (0,0), rowspan=2)

	for protein in fit_func.proteins:
		idx = fit_func.proteins.index(protein)
		xv = protein['curve'].x
		plt.plot(xv, protein['curve'].y, cmap[idx], ms=10)
		
	plt.legend( [p['curve'].ID for p in fit_func.proteins] )
	
	for protein in fit_func.proteins:
		idx = fit_func.proteins.index(protein)
		xv = np.linspace(0., np.max(protein['curve'].x), 100)
		plt.plot(xv, protein['partition'].theta(xv), cmap[idx][0]+'-', lw=2)




	plt.xlabel(fit_func.proteins[0]['curve'].denaturant_label)
	plt.ylabel('Fraction unfolded')


	# now plot the first derivative
	pk_max = []
	ax2 = plt.subplot2grid((2,2), (0,1), rowspan=1)
	for protein in fit_func.proteins:
		idx = fit_func.proteins.index(protein)
		xv = np.linspace(0.,10.,1000)
		first_deriv = np.gradient(protein['partition'].theta(xv))
		pk_max.append( (protein['n'], xv[np.argmax(np.abs(first_deriv))]) )
		plt.plot(xv, np.abs(first_deriv), cmap[idx][0]+'-', lw=2)
	plt.legend( [p['curve'].ID for p in fit_func.proteins] )
	plt.xlabel(fit_func.proteins[0]['curve'].denaturant_label)
	plt.ylabel('First derivative of fit function')


	ax2 = plt.subplot2grid((2,2), (1,1), rowspan=1)
	h,x = plot_folded(fit_func.proteins[-1]['partition'])
	for i in xrange(h.shape[1]):
		plt.plot(x, h[:,i])
	dn = iter(['_{0:d}'.format(i) for i in xrange(len(fit_func.proteins[-1]['partition'].topology))])
	plt.legend([d.name+dn.next() for d in fit_func.proteins[-1]['partition'].topology ])
	plt.xlabel(fit_func.proteins[0]['curve'].denaturant_label)
	plt.ylabel('Fraction unfolded (subpopulation)')
	plt.show()




def plot_folded(partition):
	h = np.zeros((100,partition.n))
	x = np.linspace(0.,10.,100)
	for i in xrange(partition.n):
		p = partition.subpopulation(x,i)
		h[:,i] = p 
	return h,x









if __name__ == "__main__":
	Test_EWAN()
	#test_domains()