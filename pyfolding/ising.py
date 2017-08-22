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

import sys
import time
import inspect
import numpy as np
import scipy as sp

import matplotlib.pyplot as plt

# use a global optimisation algorithm for the Ising model fitting
from scipy.optimize import differential_evolution, minimize, leastsq, curve_fit

# import PyFolding specific libraries
import constants
import models
import core

__author__ = "Alan R. Lowe"
__email__ = "a.lowe@ucl.ac.uk"


# bounds (DG_i, DG_ij, m_i, m_ij)
PARAMETER_LABELS = ('DG_i', 'DG_ij', 'm_i', 'm_ij')
#ACCEPTABLE_BOUNDS = ((-2., 15.),(-15.,2.),(-3.,-.01),(-3.,-.01))
ACCEPTABLE_BOUNDS = ((-10., 10.),(-10.,10.),(-3.,0.1),(-3.,0.1))





def free_energy_m(x, DeltaG, m_value):
	return np.exp(-(DeltaG - m_value*x) / core.temperature.RT)

def free_energy(x, DeltaG, m_value):
	return np.exp( -DeltaG / core.temperature.RT )

def dummy_free_energy(x):
	return np.exp( -0. / core.temperature.RT )









class IsingDomain(object):
	"""
	IsingDomain

	Template class for derivative Ising domains to describe a
	protein topology.

	Protein topologies can be described using, for example, a list
	of instantiated classes to represent each of the domains.

	For example:

	p = [RepeatDomain, RepeatDomain, LoopDomain, RepeatDomain]

	Variable sharing can be performed by class re-use, i.e. pointers
	to instantiated domains.

	Args:
		k_func
		t_func

	Member functions:
		tau - the tau function for this domain (coupling)
		kappa - the kappa function for this domain (intrinsic)
		q_i - calculate the weighted matrix for this domain, using the
			intrinsic kappa function and a tau function specified by the
			topology.

	Properties:
		bounds - the bounds for the optimisation algorithm
		labels - the parameter labels of the domain

	Notes:

		t_func and k_func can be intialised with m-value dependencies,
		or not, by specifying the free_energy or free_energy_m funcs.

	"""
	def __init__(self, k_func=free_energy_m, t_func=free_energy):
		if hasattr(k_func, '__call__'):
			self.__kappa = k_func
		else:
			raise TypeError("Ising domain kappa function must be a callable function")

		if hasattr(t_func, '__call__'):
			self.__tau = t_func
		else:
			raise TypeError("Ising domain tau function must be a callable function")

		self.DG_i = 2.
		self.DG_ij = -4.
		self.m_i = -.5
		self.m_ij = -.5

		self.used_variables = (0, 1, 2) # add 3 for m_ij term

		# restrict the default bounds to those used variables
		self.bounds = ACCEPTABLE_BOUNDS

		self.q_func = lambda x, tau, folded: np.matrix([[self.kappa(x)*tau(x), folded],[self.kappa(x), folded]])


	def kappa(self, x):
		""" Return the DG_intrinsic """
		return self.__kappa(x, self.DG_i, self.m_i)

	def tau(self, x):
		""" Return the DG_interface """
		return self.__tau(x, self.DG_ij, self.m_ij)

	def q_i(self, x, tau, folded=1):
		""" Return the q_i matrix """
		return self.q_func(x, tau, folded)


	@property
	def bounds(self):
		return tuple([ self.__bounds[i] for i in self.used_variables ])
	@bounds.setter
	def bounds(self, bounds):
		if not isinstance(bounds, tuple):
			raise TypeError("Fit bounds must be specified as a tuple")
		self.__bounds = bounds

	@property
	def labels(self):
		return [ PARAMETER_LABELS[i] for i in self.used_variables ]



""" Pre-defined domain topologies for GlobalFitIsing """

def list_models():
	""" Print out which models are pre-defined here """
	clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)
	fit_models = [ cls[0] for cls in clsmembers if cls[1].__bases__[0] == IsingDomain ]
	return fit_models




class HelixDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self)
		self.name = "Helix"


class RepeatDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self)
		self.name = "Repeat"


class MutantRepeatDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self)
		self.name = "MutantRepeat"


class LoopDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self)
		self.name = "Loop"


class MutantLoopDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self)
		self.name = "MutantLoop"


class CapDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self)
		self.name = "Cap"
		self.q_func = lambda x, tau, folded: np.matrix([[self.kappa(x), folded],[self.kappa(x), folded]])
		self.used_variables = (0, 2)


class MutantCapDomain(IsingDomain):
	def __init__(self):
		IsingDomain.__init__(self)
		self.name = "MutantCap"
		self.q_func = lambda x, tau, folded: np.matrix([[self.kappa(x), folded],[self.kappa(x), folded]])
		self.used_variables = (0, 2)








class GlobalFitWrapper(object):
	"""
	GlobalFitWrapper

	A wrapper for global fitting of the Ising model. Collects together all of the protein
	topologies, shares fitting paramers and calls the respective objective functions.

	Args:

	Member functions:
		append - take a protein topology and the equilibrium denaturation curve and append
				the data to the set. Automatically calculate domain reuse for variable
				sharing. Also generates the partition function (IsingPartitionFunction) for
				this new protein.

	Properties:
		bounds - the parameter bounds for this global optimisation
		domain_params - return the list of parameters required for the optimisation

	Notes:
		2016/07/30 (ARL) - Need to sort out proper variable sharing. DONE!

	"""

	def __init__(self):
		""" Initialise variables for the global optimisation """
		self.proteins = []
		self.domains = []
		self.domain_types = []


	def __call__(self, x, *args):
		""" This is the call used by the fitting function to evalutate the parameters passed in """


		if len(self.proteins) < 1:
			raise ValueError('GlobalFitIsing must have at least one curve to fit.')

		p_val = iter(x)

		# here we can set all of the parameters for the fit
		for idx, domain in enumerate(self.domains):

			# set the shared parameters for this fit
			for p in domain.labels:
				setattr(domain, p, p_val.next())


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
	"""
	General partition function object for Ising models. Based on the
	matrix formulation for a heteropolymer ising model. Briefly:

	\begin{equation}
	\begin{aligned}
	\kappa(x) &= \exp(-(\Delta G_{intrinsic} - m_{intrinsic}x) / RT) \\
	\tau(x) &= \exp(-\Delta G_{interface}) / RT) \\
	q(i) &=
	\begin{bmatrix} 0 & 1\end{bmatrix}
	\begin{bmatrix} \kappa_1\tau_{-1} & 1\\ \kappa & 1 \end{bmatrix}
	...
	\begin{bmatrix} \kappa_n\tau_{n-1} & 1\\ \kappa & 1 \end{bmatrix}
	\begin{bmatrix} 1 \\ 1\end{bmatrix} \\
	\theta &=  \frac{1}{nq(n)} \sum_{i=0}^{n}{q(i)}
	\end{aligned}
	\end{equation}

	Each protein has it's own partition function. The partition function
	relies on a fully defined topology, i.e. a 1D (Python) list of
	instantiated IsingDomain objects. The full partition function can
	be calculated from this assuming i,i-1 coupling.

	Subpartition functions can also be calculated as well as the
	fraction of folded molecules (theta).

	In the futute, it would be nice to consider other topologies,
	perhaps as graphs.

	Args:
		topology - essentially a python list of IsingDomain objects


	Member functions:
		partition - calculate the full (or sub) partition function given
			the topology
		subpartition - calculate a subpartition function (calls partition)
		theta - calculate the fraction of folded molecules at some value of x
		subpopulation - return the fraction folded for a subpopulation

	Properties:
		__len__ - the number of domains in this topology.

	Notes:
		Matrix formulation is from:
		Aksel and Barrick. Analysis of repeat-protein folding using
		nearest-neighbor statistical mechanical models.
		Methods in enzymology (2009) vol. 455 pp. 95-125

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


	def tau_func(self):
		""" Return the tau functions
		This is the sequence of tau functions, as a list of pointers to functions
		"""
		return [dummy_free_energy] + [self.topology[i].tau for i in xrange(self.n-1)]

	def partition(self, x, folded=np.array([])):
		""" Populate the partition function given the weights.
		"""

		if not folded.size:
			folded = np.ones((self.n,))

		# this is the sequence of tau functions, as a list of pointers to functions
		# tau_func = [dummy_free_energy] + [self.topology[i].tau for i in xrange(self.n-1)]
		tau_func = self.tau_func()

		# calculate the full partition function for the fully folded protein
		q_i = np.matrix([0,1])
		for i in xrange(self.n):
			domain = self.topology[i]
			q_i = q_i * domain.q_i(x, tau_func[i], folded=folded[i])
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
		q_i = self.subpartition(x, i, rev=False)
		theta = q_i / q_n
		return 1.-theta





def calculate_fit_residuals(fit_func):
	"""
	Function to calculate the residuals of a fit.
	"""
	r_res = np.array([])
	r_sq = []
	for protein in fit_func.proteins:
		y_data = protein['curve'].y
		y_fit = protein['partition'].theta( protein['curve'].x )

		r_res = np.concatenate((r_res, y_data-y_fit))
		r_sq.append( core.r_squared(y_data=y_data, y_fit=y_fit) )

	return r_res, r_sq



def calculate_error_from_jacobian(jac):
	"""
	Calculate Hessian from jacobian, and covariance from Hessian:
	covar = (J^T . J)^{-1}.

	Then calculate the error based on the SE of the variance.
	This is definitely a bit wonky at the moment!

	Args:
		jac - the jacobian matrix from the minimisation.

	Returns:
		covar - the covariance matrix for the parameters of the fit.

	Notes:
		This is **NOT** tested very well

		http://stats.stackexchange.com/questions/71154/when-an-analytical-jacobian-is-available-is-it-better-to-approximate-the-hessia
	"""

	num_params = len(np.ravel(jac))

	if np.linalg.det( np.dot(np.matrix(jac).T, np.matrix(jac)) ) == 0.:
		print "Warning: Determinant of zero indicates that this is a non-unique, poor solution!"
		return np.zeros((num_params,num_params))+np.inf

	covar = np.linalg.pinv( np.dot(np.matrix(jac).T, np.matrix(jac)) )
	return covar



class FitProgress(object):
	"""
	FitProgress

	A callback to take care of updating the user as to the fitting progress.
	Estimates the time taken per iteration and the convergence of the fit,
	to give an idea of what remains and how long it may take.
	"""

	def __init__(self, update_freq=100):
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
			print " - Fitting in progress (Iteration: {0:d}, Convergence: {1:.5E}, Timing: {2:2.2f}s per iteration) ".format(self.__iter, convergence, avg_time)






def fit_homopolymer(equilibrium_curves=[], topologies=[], p0=[5, 3.3,.1,-5.], **kwargs):
	"""
	Fit a homopolymer model to a dataset.

	TODO: Improve this.

	"""

	global_fit = core.GlobalFit()
	global_fit.fit_funcs = [models.HomozipperIsingEquilibrium for i in xrange(len(equilibrium_curves))]
	global_fit.constants = [(('n',n),) for n in topologies]

	global_fit.x = [p.x for p in equilibrium_curves]

	x = np.concatenate([p.x for p in equilibrium_curves])
	y = np.concatenate([p.y for p in equilibrium_curves])

	# fit the curve
	out, covar = curve_fit(global_fit, x, y, p0=p0, bounds=((0,0,-1.,-10.),(21,10.,1.,0)) )

	# somewhere to store the results
	results = []
	covar = covar[1:,1:]

	x = np.linspace(0.,10.,100)

	# calculate errors
	for i, protein in enumerate(equilibrium_curves):

		# calculate the fit
		n = topologies[i]
		fit_args = [n] + out[1:].tolist()
		fit_y = global_fit.fit_funcs[i](global_fit.x[i], *fit_args)

		result = core.FitResult(fit_name="Homopolymer Ising Model", fit_args=global_fit.fit_funcs[0].fit_func_args[1:])
		result.ID = protein.ID
		result.fit_params = out[1:].tolist()
		result.method = "scipy.optimize.curve_fit"
		#result.y = protein['partition'].theta( result.x )
		result.y = global_fit.fit_funcs[i](x, *fit_args)
		result.covar = covar
		result.residuals = protein.y - fit_y
		#result.r_squared = r_squared[i]

		results.append( result )

	print '\nFitting results: '
	for r_arg, r_val, r_err in results[0].details:
	 	print u"{0:s}: {1:2.5f} \u00B1 {2:2.5f} ".format(r_arg, r_val, r_err)


	# plot some of the results
	plt.figure(figsize=(14,8))
	for i, p in enumerate(equilibrium_curves):
		plt.plot(p.x, p.y, 'o', x, results[i].y, '-')
	plt.xlabel(p.denaturant_label, fontsize=constants.FONT_SIZE)
	plt.ylabel('Fraction unfolded', fontsize=constants.FONT_SIZE)
	plt.show()


	if 'save' in kwargs:

		# import the CSV writer
		import csv

		save_filename = kwargs['save']
		if not isinstance(save_filename, basestring):
			raise TypeError("Save path must be a string")

		if not save_filename.endswith(('.csv', '.CSV')):
			save_filename = save_filename+".csv"

		# make a list of the field names for the csv file
		fnames = [equilibrium_curves[0].denaturant_label]
		for result in results:
			fnames.append(result.ID)

		# now put iterators into a dictionary

		csv_results = {fnames[0]: iter(x) }
		for result in results:
			# NOTE: we should update this to cover the entire x range
			csv_results[result.ID] = iter(result.y)

		# note to interested reader - this is a really weird way of doing this!
		with open(save_filename, 'wb') as csvfile:
			r = csv.DictWriter(csvfile, fieldnames=fnames, dialect=csv.excel_tab, delimiter=',')
			r.writeheader()

			for i in xrange(len(x)):
				h = {k: csv_results[k].next() for k in csv_results.keys()}
				r.writerow(h)

		print "Written out .csv file of fits..."

	return results
















def fit_heteropolymer(equilibrium_curves=[], topologies=[], popsize=10, tol=1e-8, maxiter=None, **kwargs):
	"""

	An example script to fit a series of data sets to a heteropolymer ising model.

	Args:
		equilibrium_curves
		topologies
		popsize
		tol
		save


	Notes:
		Optimisation is performed using differential evolution (a GA)
		http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html#scipy.optimize.differential_evolution

		TODO (2017/04/10) - Implement the saving using the core save function

	"""

	# do some parsing of the input
	if not isinstance(equilibrium_curves, list):
		raise TypeError('equilibrium_curves must be a list of EquilibriumDenaturationCurve type')

	if not isinstance(topologies, list):
		raise TypeError('topologies must be a list of IsingDomain type')

	results = []
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

	# perform the actual fitting operation
	r = differential_evolution(fit_func, fit_func.bounds, disp=False, popsize=popsize, tol=tol, callback=callback, maxiter=maxiter)

	if not r.success:
		print "Could not find a solution..."
		return None

	# calculate the errors on each of the parameters
	if hasattr(r, 'jac'):
		# the fit has the jacobian
		jac = r.jac
	else:
		r_min = minimize(fit_func, np.copy(r.x), method='L-BFGS-B', bounds=fit_func.bounds)
		jac = r_min.jac

	# calculate the errors
	r_res, r_squared = calculate_fit_residuals(fit_func)
	r_cov = calculate_error_from_jacobian(jac)

	for i, protein in enumerate(fit_func.proteins):
		result = core.FitResult(fit_name="Heteropolymer Ising Model", fit_args=fit_func.domain_params)
		result.ID = protein['curve'].ID
		result.fit_params = r.x.tolist()
		result.method = "scipy.optimize.differential_evolution"
		result.y = protein['partition'].theta( result.x )
		result.covar = r_cov
		result.residuals = r_res		# NOTE: these are the total residuals as opposed to the per curve residuals
		result.r_squared = r_squared[i]

		results.append( result )

	print '\nFitting results (NOTE: Careful with the errors here): '
	for r_arg, r_val, r_err in results[0].details:
	 	print u"{0:s}: {1:2.5f} \u00B1 {2:2.5f} ".format(r_arg, r_val, r_err)
	for result in results:
		print u"{0:s} R^2: {1:2.5f}".format(result.ID, result.r_squared)

	# now plot the output
	plot_domains(topologies, labels=[protein.ID for protein in equilibrium_curves])
	plot_Ising(fit_func)
	plot_folded(fit_func.proteins[-1]['partition'])


	if "save" in kwargs:

		# import the CSV writer
		import csv

		save_filename = kwargs['save']
		if not isinstance(save_filename, basestring):
			raise TypeError("Save path must be a string")

		if not save_filename.endswith(('.csv', '.CSV')):
			save_filename = save_filename+".csv"

		# make a list of the field names for the csv file
		fnames = [fit_func.proteins[0]['curve'].denaturant_label]
		for p in equilibrium_curves:
			fnames.append(p.ID)

		# now put iterators into a dictionary
		x = np.linspace(0.,10.,100)
		csv_results = {fnames[0]: iter(x) }
		for protein in fit_func.proteins:
			csv_results[protein['curve'].ID] = iter( protein['partition'].theta(x) )

		# note to interested reader - this is a really weird way of doing this!
		with open(save_filename, 'wb') as csvfile:
			r = csv.DictWriter(csvfile, fieldnames=fnames, dialect=csv.excel_tab, delimiter=',')
			r.writeheader()

			for i in xrange(len(x)):
				h = {k: csv_results[k].next() for k in csv_results.keys()}
				r.writerow(h)

		print "Written out .csv file of fits..."


	return results


def plot_Ising(fit_func):
	"""
	Function to plot fitted Ising model data.

	"""

	cmap = ['ro', 'mo', 'go', 'co', 'bo', 'ko', 'rv', 'mv', 'gv', 'cv', 'bv',
			'kv', 'rs', 'ms', 'gs', 'cs', 'bs', 'ks', 'r.', 'm.', 'g.', 'c.',
			'b.', 'k.']

	# plot the fits
	plt.figure(figsize=(14,8))
	ax1 = plt.subplot2grid((2,2), (0,0), rowspan=2)

	for protein in fit_func.proteins:
		idx = fit_func.proteins.index(protein)
		xv = protein['curve'].x
		ax1.plot(xv, protein['curve'].y, cmap[idx % len(cmap)], ms=4)

	#plt.legend( [p['curve'].ID for p in fit_func.proteins], loc='upper left')

	for protein in fit_func.proteins:
		idx = fit_func.proteins.index(protein)
		xv = np.linspace(0., np.max(protein['curve'].x), 100)
		ax1.plot(xv, protein['partition'].theta(xv), cmap[idx % len(cmap)][0]+'-',
			lw=2, label=protein['curve'].ID)

	ax1.set_xlabel(fit_func.proteins[0]['curve'].denaturant_label)
	ax1.set_ylabel('Fraction unfolded')


	# now plot the first derivative
	pk_max = []
	ax2 = plt.subplot2grid((2,2), (0,1), rowspan=1)
	for protein in fit_func.proteins:
		idx = fit_func.proteins.index(protein)
		xv = np.linspace(0.,10.,1000)
		first_deriv = np.gradient(protein['partition'].theta(xv))
		pk_max.append( (protein['n'], xv[np.argmax(np.abs(first_deriv))]) )
		ax2.plot(xv, np.abs(first_deriv), cmap[idx % len(cmap)][0]+'-', lw=2,
			label=protein['curve'].ID)

	ax2.set_xlabel(fit_func.proteins[0]['curve'].denaturant_label)
	ax2.set_ylabel('First derivative of fit function')


	ax3 = plt.subplot2grid((2,2), (1,1), rowspan=1)
	h,x = plot_folded(fit_func.proteins[-1]['partition'])
	dn = iter(['{0:s}_{1:d}'.format(d.name,i) for i,d in enumerate(fit_func.proteins[-1]['partition'].topology)])
	for i in xrange(h.shape[1]):
		plt.plot(x, h[:,i], cmap[i % len(cmap)]+'-', lw=2, markersize=4,
			label=dn.next())
	ax3.set_xlabel(fit_func.proteins[0]['curve'].denaturant_label)
	ax3.set_ylabel('Fraction unfolded (subpopulation)')

	# do some formatting with the smaller plots to fit legends beside them
	# http://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
	box2 = ax2.get_position()
	ax2.set_position([box2.x0, box2.y0, box2.width*0.8, box2.height])
	box3 = ax3.get_position()
	ax3.set_position([box3.x0, box3.y0, box3.width*0.8, box3.height])

	# Put a legend below current axis
	ax2.legend(loc='center left', bbox_to_anchor=(1., 0.5))
	ax3.legend(loc='center left', bbox_to_anchor=(1., 0.5))

	plt.show()




def plot_folded(partition):
	h = np.zeros((100,partition.n))
	x = np.linspace(0.,10.,100)
	for i in xrange(partition.n):
		p = partition.subpopulation(x,i)
		h[:,i] = p
	return h,x


def collapse_topology(topology):
	"""
	Collapse a topology into a compact representation.
	"""

	# make a copy so that we don't do this in place
	topology_full = list(topology)
	compact = []
	counter = 0

	domain = topology_full.pop(0)

	while topology_full:
		next_domain = topology_full.pop(0)
		if next_domain is domain:
			counter += 1
		else:
			compact.append((domain, counter+1))
			counter = 0
			domain = next_domain

	# clean up
	compact.append((domain, counter+1))

	return compact




def brace(x,y, rev=1, length=0.95):
	xb = [x+(rev*.4), x+(rev*.5), x+(rev*.5), x+(rev*.4)]
	yb = [y+length/2, y+length/2, y-length/2, y-length/2]
	return xb,yb




def plot_domains(topologies, labels=None, collapse=False, **kwargs):
	"""
	Function to generate a pretty plot of the domain architecture of
	the various protein topologies presented.


	Args:
		topologies: the protein topologies
		labels: protein names
		fold: a boolean that collapses sequences of similar domains for smaller
			representation.

	"""

	from matplotlib.patches import Patch

	if not labels:
		labels = ['Protein {0:d}'.format(i) for i in xrange(len(topologies))]

	if 'fold' in kwargs:
		raise DeprecationWarning('Fold keyword is being replaced with collapse.')
		collapse = kwargs['fold']


	# collapse the topologies
	compact = []
	domain_types = set()

	for l,t in zip(labels, topologies):

		domain_types.update([d().name for d in t])

		c = collapse_topology(t)
		if not collapse:
			tc = []
			for d,n in c: tc += [(d,1)]*n
			c = tc

		compact.append((l, c))


	# set up the plotting
	domain_types = list(domain_types)

	d_color = lambda name: cmap[ domain_types.index(name) ]
	cmap = ['r', 'k', 'g', 'b', 'c', 'm', 'y']

	# now we want to plot these
	plt.figure(figsize=(14, len(topologies)*1.5))
	ax = plt.subplot(111)

	for y, (protein, topology) in enumerate(compact):
		for x, (domain, d_cnt) in enumerate(topology):

			name = domain().name
			c = plt.Circle((x*1.5, y), 0.45, edgecolor=d_color(name),
						facecolor='w', label=name)
	 		ax.add_artist(c)


	 		# if we're folding, then plot braces to show that
	 		if d_cnt > 1:
	 			lb_x, lb_y = brace(x*1.5,y, rev=-1)
	 			rb_x, rb_y = brace(x*1.5,y, rev=1)
	 			b = plt.Line2D(lb_x, lb_y, color='k')
	 			ax.add_artist(b)
	 			b = plt.Line2D(rb_x, rb_y, color='k')
	 			ax.add_artist(b)
	 			ax.text(x*1.5+.55, y-.4, str(d_cnt), color='k')


	 	 	i_str = str(x)
			ij_str = str(x)

			# add on any labels:
			if hasattr(domain(), 'm_i'):
				ax.text(x*1.5,y+.1,'$m^{'+i_str+'}_{i}$',
					horizontalalignment='center', fontsize=18, color=d_color(name))
			if hasattr(domain(), 'DG_i'):
				ax.text(x*1.5,y-.3,'$\Delta G^{'+i_str+'}_{i}$',
					horizontalalignment='center', fontsize=18, color=d_color(name))

			# draw arrows to represent interactions
			if hasattr(domain(), 'DG_ij') and x!=len(topology)-1:
				ax.text(x*1.5+.9,y,'$\Delta G^{'+ij_str+'}_{ij}$',
					horizontalalignment='center', fontsize=12, rotation=90,
					color=d_color(name))



	ax.set_yticks(np.arange(len(topologies)), minor=False)
	ax.set_xticklabels([], minor=False, rotation='vertical')
	ax.set_yticklabels(labels, minor=False)

	# make the legend entries
	l = [ Patch(edgecolor=d_color(d), facecolor='w', label=d) for d in domain_types ]


	ax.set_xlim([-1.,11.])
	ax.set_ylim([-1., len(topologies)])
	ax.set_aspect('equal', adjustable='box')
	plt.legend(handles=l, loc='lower right')
	plt.title('Ising model domain topologies')

	if 'save' in kwargs:
		save_filename = kwargs['save']
		if not isinstance(save_filename, basestring):
			raise TypeError("Save path must be a string")
		if not save_filename.endswith(('.pdf', '.PDF')):
			save_filename = save_filename+".pdf"

		plt.savefig(save_filename, dpi=144., pad_inches=0)

	plt.show()
	return






if __name__ == "__main__":
	pass
