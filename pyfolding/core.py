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

Lowe, A.R. 2015

"""



import os
import csv
import inspect

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from scipy import optimize
from collections import OrderedDict

import constants

__author__ = "Alan R. Lowe"
__email__ = "a.lowe@ucl.ac.uk"
__version__ = 0.03



# set some global plotting params
plt.rc('xtick', labelsize=constants.LABEL_SIZE)
plt.rc('ytick', labelsize=constants.LABEL_SIZE)

# by default turn off autoscrolling if it exists
import utils
utils.disable_autoscroll()






class __Temperature(object):
	""" Maintain temperature information across all functions. Meant to be a
	wrapper for a global variable, temperature which is used across different
	models.

	Args:
		temperature: set the temperature in celsius

	Properties:
		temperature: return the temperature in celsius
		RT:	return the product of the ideal gas constant and the temperature

	Notes:
		This changes the temperature globally, so be careful when using it,
		since all subsequent calculations will use this temperature.

		This is also not to be used by USERS!
	"""
	def __init__(self):
		self.__temperature = constants.TEMPERATURE_CELSIUS

	@property
	def temperature(self): return self.__temperature
	@temperature.setter
	def temperature(self, value=constants.TEMPERATURE_CELSIUS):
		import numbers
		if not isinstance(value, numbers.Real):
			return TypeError("Temperature must be specified as a number")
		if value < 0 or value > 100:
			raise ValueError("Temperature ({0:2.2f}) is not in valid range (0-100 degrees C)".format(value))
		self.__temperature = value

	@property
	def RT(self):
		return constants.IDEAL_GAS_CONSTANT_KCAL * (constants.ZERO_KELVIN + self.temperature)

temperature = __Temperature()


def set_temperature(value=constants.TEMPERATURE_CELSIUS):
	""" Set the temperature.

	Args:
		temperature: set the temperature in celsius

	Returns:
		None

	Usage:
		>> pyfolding.set_temperature( 10.2 )
	"""
	temperature.temperature = value
	print u"Set temperature to {0:2.2f}\u00B0C".format(value)
	print "(NOTE: Careful, this sets the temperature for all subsequent calculations)"



"""
===========================================================
TEST FUNCTION
===========================================================
"""


def test(protein_ID='Simulated protein'):
	"""
	Test function to make sure that PyFolding is installed correctly
	and functioning as it should. Generates a simulated data set
	using known parameters and noise, and then fits and plots the
	data comparing these to the ground truth.
	"""

	import models

	# initialise the data structures
	chevron = Chevron(ID=protein_ID)
	equilibrium = EquilibriumDenaturationCurve(ID=protein_ID)

	acceptible_error = 1e-2
	truth = {'eq':[1.5, 5.], 'kin': [100., 1., 0.005, 1.]}


	# denaturant concentrations
	den = np.linspace(0.,10.,100)

	# generate a two-state equilibrium curve, with Gaussian noise
	# alpha_f, beta_f, alpha_u, beta_u, m, d50
	eq_model = models.TwoStateEquilibrium()
	eq_raw = eq_model.fit_func(den, *truth['eq'])
	eq_sim = eq_raw + np.random.randn(100,)*0.0001

	equilibrium.denaturant_label = "[Denaturant] (M)"
	equilibrium.curves = ['e1']
	equilibrium.denaturant = {'e1': den}
	equilibrium.signal = {'e1': eq_sim}
	equilibrium.fit_func = models.TwoStateEquilibrium


	# generate a two-state chevron curve, with Gaussian noise
	# kf, mf, ku, mu
	kin_model = models.TwoStateChevron()
	kin_raw = kin_model.fit_func(den, *truth['kin'])
	kin_sim = np.exp( np.log(kin_raw) + np.random.randn(100,)*0.001 )


	chevron.denaturant_label = "[denaturant] (M)"
	chevron.phases = ['k1']
	chevron.denaturant = {'k1': den}
	chevron.rates = {'k1': kin_sim}
	chevron.fit_func = models.TwoStateChevron

	# fit the equilibrium data to a two-state model
	equilibrium.fit()

	# use the midpoint (D_50) of the equilibrium curve as the kinetic midpoint
	chevron.midpoint = equilibrium.midpoint

	# now fit the chevron to a two-state model
	chevron.fit()

	# get the parameters and check that they are the same as the
	# ground truth set

	for p_truth, p_fit in zip(truth['eq'], equilibrium.fit_params):
		if (p_truth - p_fit)**2 > acceptible_error:
			raise ValueError("PyFolding self-test failed. Fitting error ({0:f}) exceeds \
				bounds ({1:f}) \n".format((p_truth - p_fit)**2, acceptible_error))

	for p_truth, p_fit in zip(truth['kin'], chevron.fit_params):
		if (p_truth - p_fit)**2 > acceptible_error:
			raise ValueError("PyFolding self-test failed. Fitting error ({0:f}) exceeds \
				bounds ({1:f}) \n".format((p_truth - p_fit)**2, acceptible_error))


	print 'Test completed!'

	# plot the output
	plot_figure(equilibrium, chevron, display=True)




"""
===========================================================
FILE I/O OPERATIONS
===========================================================
"""





def check_filename(directory, filename):
	""" Check the filename for consistency.
	"""

	if not isinstance(directory, basestring) or not isinstance(filename, basestring):
		raise TypeError('Pyfolding expects a filename as a string')

	if not filename.lower().endswith(('.csv', '.CSV')):
		raise IOError('PyFolding expects a .CSV file as input: {0:s}'.format(filename))

	if not os.path.exists(os.path.join(directory, filename)):
		raise IOError('PyFolding could not find the file: {0:s}'.format(os.path.join(directory, filename)))




def read_kinetic_data(directory=None, filename=None):
	""" Read in kinetic data in the form of an .csv worksheet. It
	should be arranged such that each file is a different protein,
	and columns represent the following:

	[den] k1 k2 ...

	This function then returns a chevron object with the data
	"""

	# check the filename
	check_filename(directory, filename)

	protein_ID, ext = os.path.splitext(filename)
	chevron = Chevron(ID=protein_ID)

	with open(os.path.join(directory, filename), 'rU') as chevron_file:
		chevron_reader = csv.reader(chevron_file, delimiter=',', quotechar='|')
		header = chevron_reader.next()

		# set up the chevron with the correct number of phases
		chevron.denaturant_label = header[0]
		chevron.phases = header[1:]
		chevron.denaturant = {phase:[] for phase in chevron.phases}
		chevron.rates = {phase:[] for phase in chevron.phases}

		# now group all of the data and insert into object
		for row in chevron_reader:
			den_conc = float( row.pop(0) )
			for phase, rate in zip(chevron.phases, row):
				if rate:
					chevron.denaturant[phase].append(den_conc)
					chevron.rates[phase].append(float(rate))

	return chevron


def read_equilibrium_data(directory=None, filename=None):
	""" Read in an equilbrium denaturation curve from a .csv
	worksheet. It should be arranged such that each file is a
	different protein, and columns represent the following:

	[den] unfolding

	This function then returns an equilbrium curve object.
	"""

	# check the filename
	check_filename(directory, filename)

	protein_ID, ext = os.path.splitext(filename)
	equilibrium = EquilibriumDenaturationCurve(ID=protein_ID)

	with open(os.path.join(directory, filename), 'rU') as equilbrium_file:
		equilibrium_reader = csv.reader(equilbrium_file, delimiter=',', quotechar='|')
		header = equilibrium_reader.next()

		# set up the chevron with the correct number of phases
		equilibrium.denaturant_label = header[0]
		equilibrium.curves = header[1:]
		equilibrium.denaturant = {curve:[] for curve in equilibrium.curves}
		equilibrium.signal = {curve:[] for curve in equilibrium.curves}

		# now group all of the data and insert into object
		for row in equilibrium_reader:
			den_conc = float( row.pop(0) )
			for curve, signal in zip(equilibrium.curves, row):
				if signal:
					equilibrium.denaturant[curve].append(den_conc)
					equilibrium.signal[curve].append(float(signal))

	return equilibrium









"""
===========================================================
BASE CLASSES
===========================================================
"""



class Protein(object):
	""" Protein wrapper object
	"""
	def __init__(self, ID=None):
		self.ID = ID
		self.chevron = None
		self.equilibrium = None

	@property
	def deltaG(self): return self.equilibrium.deltaG

	@property
	def kf_H20(self): return self.chevron.results.y[0]






class FitResult(object):
	""" Fitting result object.

	This is an internal class that collates fit results and enables calculation
	of errors, residuals and other fun things.

	Args:
		name: a name for the fit result (e.g. TwoStateChevron)
		fit_args: the fit function arguments
		ID: The identifier of the protein

	Properties:
		method: the name of the optimisation algorithm used
		errors: the calculated errors (SEM) for the fit arguments
		details: return a zipped list of argument, value, error tuples
		standard_error: the standard_error of the overall fit
		covar: covariance matrix following optimisation
		residuals: residuals of fit to data
		r_squared: r^2 value for the fit

	Members:
		display: return a formatted output of the fitting statistics

	Notes:
		TODO(arl): implement an export function

	"""

	def __init__(self, fit_name=None, fit_args=None):
		self.ID = None
		self.fit_args = fit_args
		self.name = fit_name
		self.x = np.linspace(0., 10., 100)
		self.y = None
		self.covar = None
		self.residuals = None
		self.r_squared = None

		self.fit_errors = None
		self.fit_params = None

		self.__method = "scipy.optimize.curve_fit"

	@property
	def method(self): return self.__method
	@method.setter
	def method(self, method=None):
		if not isinstance(method, basestring):
			raise TypeError("FitResult: Method must be a string")
		self.__method = method

	def display(self):
		""" Print the errors and fit values """

		table_width = max([len("Model: "+self.name), len(" Fitting results "), 50])

		print u"="*table_width
		print u" Fitting results "
		print u"="*table_width
		if self.ID: print u"ID: {0:s}".format(self.ID)
		print u"Model: {0:s}".format(self.name)
		print u"Method: {0:s} \n".format(self.method)
		for e in self.details:
			fit_arg, fit_val, fit_err = e
			print u"{0:s}: {1:2.5f} \u00B1 {2:2.5f}".format(fit_arg, fit_val, fit_err)

		print u"-"*table_width
		print u"R^2: {0:2.5f}".format(self.r_squared)
		print u"="*table_width

	@property
	def errors(self):
		""" Calculates the SEM of the fitted parameters, using the (hopefully
			well formed) covariance matrix from the curve fitting algorithm.
		"""

		if not isinstance(self.covar, np.ndarray):
			raise ValueError("FitResult: Covariance matrix is not defined")

		if not isinstance(self.residuals, np.ndarray):
			raise ValueError("FitResult: Residuals are not defined")

		num_params = len(self.fit_args)
		errors = [np.sqrt(float(self.covar[p,p]) * np.var(self.residuals)) / np.sqrt(1.*len(self.residuals)) for p in xrange(num_params)]
		return errors

	@property
	def details(self):
		""" Return a zipped list of the fit arguments, values and errors """
		return zip(self.fit_args, self.fit_params, self.errors)

	@property
	def standard_error(self):
		""" Return the standard error of the fit """
		return np.std(self.residuals) / np.sqrt(1.*len(self.residuals))

	def export(self, filename):
		raise NotImplementedError










class FoldingData(object):
	""" Base class fo chevrons and equilibrium denaturation
	curves. Takes care of common functions such as fitting
	of models.
	"""
	def __init__(self):
		self.__fit_func = None
		self.fit_params = None
		self.fit_covar = None
		self.__fit = None
		self.__fit_residuals = None
		self.components = None

	@property
	def fit_func(self): return self.__fit_func.name
	@fit_func.setter
	def fit_func(self, fit_func=None):
		if hasattr(fit_func, "__call__"):
			self.__fit_func = fit_func()
		else:
			raise AttributeError("Fit function must be callable")

	@property
	def fit_func_args(self):
		if self.__fit_func:
			return self.__fit_func.fit_func_args

	@property
	def results(self):
		return self.__fit
	@results.setter
	def results(self, result):
		if not isinstance(result, FitResult):
			raise TypeError("Results must be of type FitResult")

		print "Warning: overwriting fit result for {0:s}".format(self)
		self.__fit = result

	def fit(self, p0=None, constants=None):
		""" Fit the data to the defined model. Use p0 to introduce the estimated
		 start values.

		 TODO(arl): Should sort out how constants are passed to curve fit,
		 since they will affect errors in current form.
		"""
		if self.__fit_func:

			# reset components
			self.components = None

			# set the default fitting parameters
			if not p0: p0 = self.__fit_func.default_params

			# set any constants
			if constants:
				# TODO: error checking/parsing here
				self.__fit_func.constants = constants

			# perform the actual fitting
			try:
				out = optimize.curve_fit(self.__fit_func, self.x, self.y, p0=p0,
					maxfev=20000)
			except RuntimeWarning:
				raise Exception("Optimisation could not complete. Try adjusting"
				 		" your fitting parameters (p0)")

			# create a results structure
			self.__fit = FitResult(fit_name=self.fit_func, fit_args=self.fit_func_args)
			self.__fit.ID = self.ID
			self.__fit.method = "scipy.optimize.curve_fit"

			# use get_fit_params in case we have any constants defined
			self.__fit.fit_params = np.array( self.__fit_func.get_fit_params(self.x, *list(out[0])) )
			self.__fit.y = self.__fit_func.fit_func(self.__fit.x, *list(self.__fit.fit_params))
			self.__fit.covar = out[1]

			fit_y_data_x = self.__fit_func.fit_func(self.x, *list(self.__fit.fit_params))
			self.__fit.residuals = self.y_raw - fit_y_data_x
			self.__fit.r_squared = r_squared(y_data=self.y_raw, y_fit=fit_y_data_x)

			# TODO(arl): this is obsolete now - we should be storing the FitResult
			# store locally ?
			self.fit_params = np.array( self.__fit_func.get_fit_params(self.x, *list(out[0])) )
			self.fit_covar = out[1]

			if hasattr(self.__fit_func, "components"):
				self.components = self.__fit_func.components(np.linspace(0.,10.,100), *list(self.fit_params))

		else:
			raise AttributeError("Fit function must be defined")


		self.__fit.display()


	@property
	def fitted_x(self):
		raise DeprecationWarning("Warning, this feature will be deprecated soon. Use .results for the fit results")
		return np.linspace(0., 10., 100)

	@property
	def fitted(self):
		raise DeprecationWarning("Warning, this feature will be deprecated soon. Use .results for the fit results")
		if isinstance(self.fit_params, np.ndarray):
			x = self.fitted_x
			return np.copy( self.__fit_func.fit_func(x, *list(self.fit_params)) )

	def print_fit_params(self):
		if isinstance(self.fit_params, np.ndarray):
			print self.fit_params

	def plot(self, **kwargs):
		""" Plot a simple figure of the data, this is context dependent
		title='', marker='wo', display_fit=True
		"""

		# make this cleaner by calling an independent function. User can also
		# call these functions
		if isinstance(self, Chevron):
			plot_chevron(self, **kwargs)
		else:
			plot_equilibrium(self, **kwargs)


	def save_fit(self, filename):
		d = [('x',self.results.x), ('y',self.results.y)]
		# order this dictionary so that we have x first
		data = OrderedDict( d )
		utils.write_CSV(filename, data)







class Chevron(FoldingData):
	""" Chevron plot for protein folding kinetics.

	Args:

	Methods:

	Notes:

	"""
	def __init__(self, ID=None):
		FoldingData.__init__(self)
		self.ID = ID
		self.denaturant_label = None
		self.denaturant = None
		self.rates = None
		self.phases = None
		self.__midpoint = None


	@property
	def x(self): return np.array(self.denaturant[self.phases[0]])
	@property
	def y(self): return np.array(np.log(self.rates[self.phases[0]]))

	@property
	def y_raw(self): return np.array(self.rates[self.phases[0]])


	@property
	def midpoint(self):
		""" Return a calculated midpoint for the chevron. Unless
		we have set one using equilibrium data.
		"""
		if not self.__midpoint and self.denaturant:
			return self.denaturant['k1'][ np.argmin(self.rates['k1']) ]
		else:
			return self.__midpoint

	@midpoint.setter
	def midpoint(self, midpoint=0.0):
		if isinstance(midpoint, float):
			if midpoint>0. and midpoint<10.: self.__midpoint = midpoint
		else:
			raise Exception("Midpoint must be a float and 0<x<10")

	def unfolding_limb(self, phase=None):
		""" Return only the unfolding limb data
		"""
		if not phase:
			phase = self.phases[0]
		elif phase not in self.phases:
			return None
		denaturant, rates = [], []
		for d,r in zip(self.denaturant[phase], self.rate(phase)):
			if d > self.midpoint:
				denaturant.append(d)
				rates.append(r)

		return denaturant, rates

	def refolding_limb(self, phase=None):
		""" Return only the refolding limb data
		"""
		if not phase:
			phase = self.phases[0]
		elif phase not in self.phases:
			return None
		denaturant, rates = [], []
		for d,r in zip(self.denaturant[phase], self.rate(phase)):
			if d <= self.midpoint:
				denaturant.append(d)
				rates.append(r)

		return denaturant, rates

	def chevron(self, phase=None):
		""" Return the entire phase of a chevron
		"""
		if not phase:
			phase = self.phases[0]
		elif phase not in self.phases:
			return None
		return self.denaturant[phase], self.rate(phase)

	def rate(self, phase=None):
		return np.log(self.rates[phase])







class EquilibriumDenaturationCurve(FoldingData):
	""" Equilibrium Denaturation curve
	Args:

	Methods:

	Notes:

	"""
	def __init__(self, ID=None):
		FoldingData.__init__(self)
		self.ID = ID
		self.curves = None
		self.denaturant_label = None
		self.denaturant = None
		self.signal = None


	@property
	def x(self): return np.array(self.denaturant[self.curves[0]])
	@property
	def y(self): return np.array(self.signal[self.curves[0]])

	@property
	def y_raw(self): return self.y

	@property
	def normalised(self):
		""" TODO(arl): Return a normalised equilbrium curve.
		"""
		raise NotImplementedError

	@property
	def m_value(self):
		if isinstance(self.fit_params, np.ndarray):
			return self.fit_params[ self.fit_func_args.index('m') ]
		return None

	@property
	def midpoint(self):
		if isinstance(self.fit_params, np.ndarray):
			return self.fit_params[ self.fit_func_args.index('d50') ]
		else:
			return None

	@property
	def two_state(self):
		""" Return whether this is a two state model or not """
		return 'd50' in self.fit_func_args

	def point(self, fraction_folded=0.5):
		""" Return the denaturant concentration for a particular
		fraction folded. Assumes a two-state transition since I
		had to derive this equation by hand.
		"""
		if self.m_value and self.midpoint:
			if fraction_folded<0. or fraction_folded>1.:
				raise ValueError("Fraction folded must be in the range 0.<x<1.")
			return (np.log((1.-fraction_folded)/fraction_folded) / self.m_value) + self.midpoint
		else:
			return None

	@property
	def deltaG(self):
		""" Return the deltaG value based on the fit of the data """
		if self.m_value and self.midpoint:
			return self.m_value * self.midpoint
		else:
			return None







def FIT_ERROR(x):
	""" Return a generic fit error """
	if isinstance(x, np.ndarray):
		return np.ones(x.shape)*constants.FITTING_PENALTY
	else:
		return None








class GlobalFit(object):
	""" Wrapper function to perform global fitting.  This acts as a wrapper for
	multiple FitModels, enabling the user to pair datasets and models and share
	data or arguments.

	Args:
		x: concatenated x data
		y: concatenated y data

	Properties:
		fit_funcs: the fit functions
		constants: constants for the fitting

	Members:
		__call__: evaluates the fit functions

	Notes:

	"""
	def __init__(self):
		self.x = []
		self.y = []
		self.__fit_funcs = []

	@property
	def fit_funcs(self): return self.__fit_funcs
	@fit_funcs.setter
	def fit_funcs(self, fit_funcs):
		for fit_func in fit_funcs:
			if not hasattr(fit_func, "__call__"): continue
			# append it and instantiate it
			self.__fit_funcs.append( fit_func() )

	@property
	def constants(self):
		return [f.constants for f in self.__fit_funcs]
	@constants.setter
	def constants(self, constants=None):
		if len(constants) != len(self.__fit_funcs):
			raise ValueError("Number of constants should be the same as number of fit functions")

		for constant, fit_func in zip(constants, self.__fit_funcs):
			fit_func.constants = constant

	def __call__(self, *args):
		""" Dummy call for all fit functions """
		x = args[0]
		fit_args = args[1:]
		ret = np.array(())
		for i, fit_func in enumerate(self.fit_funcs):
			x_this = np.array( self.x[i] )
			ret = np.append( ret, fit_func(x_this, *fit_args) )
		return ret









class FitModel(object):
	""" FitModel class

	A generic fit model to enable a common core set of functions
	but specific new member functions to be enabled in derived
	classes.

	Can define parameters in this manner:
	('kf',0), ('mf',1), ... in order to enable paramater sharing
	in global fitting. By default the model just gets the params
	in the order they are defined in the function defininition

	Note: this must be subclassed to work.

	Args:
		constants: constants for fitting

	Properties:
		params: the parameters of the fit model
		name: the name of the fit model
		default_params: the default starting parameters for the fit
		fit_func_args: the names of the fit function arguments
		equation: a LaTeX formatted string of the model

	Methods:
		__call__: evaluates the fit function with the given args
		print_equation: (static) prints the equation to the stdout
		fit_func: (not defined) the actual fit function
		error_func: (not defined) the error function

	Notes:

	"""
	def __init__(self):
		self.__params = None
		self.__param_names = None
		self.__default_params = None

		self.fit_params = None
		self.fit_covar = None
		self.constants = None

		# has this model been verified
		self.verified = False

	@property
	def params(self): return self.__params
	@params.setter
	def params(self, params=None):
		if isinstance(params, tuple):
			self.__params, self.__param_names = [], []
			for key,value in params:
				self.__param_names.append(key)
				self.__params.append(value)
		else:
			raise Warning("Fit parameters must be a tuple")


	@property
	def name(self): return self.__class__.__name__

	def __call__(self, x, *args):
		""" Parse the fit arguments and pass onto the
		fitting function
		"""
		fit_args = self.get_fit_params(x, *args)
		return self.error_func( self.fit_func(x, *fit_args) )


	def fit_func(self, x, *args):
		""" The fit function should be defined here """
		raise Exception("Fit function must be defined")

	def error_func(self, y):
		""" The error function should be defined here """
		return y

	def get_fit_params(self, x, *args):
		fit_args = [args[v] for v in self.__params]
		# if we have constants replace the arguments
		# with the constants
		if self.constants:
			for arg, constant in self.constants:
				if arg in self.__param_names:
					idx = self.__params[ self.__param_names.index(arg) ]
					fit_args[idx] = constant
		return fit_args

	@property
	def default_params(self):
		""" Give back either the set starting parameters,
		or set all to 1.0
		"""
		if isinstance(self.__default_params, np.ndarray):
			return self.__default_params
		else:
			return np.ones((len(self.params),1))
	@default_params.setter
	def default_params(self, params):
		if isinstance(params, np.ndarray):
			self.__default_params = params

	@property
	def fit_func_args(self):
		return inspect.getargspec(self.fit_func).args[2:]

	@property
	def equation(self):
		raise NotImplementedError

	@staticmethod
	def print_equation(self):
		from IPython.display import display, Math, Latex
		display(Math(self.equation))
		return None









def r_squared(y_data=None, y_fit=None):
	return 1. - np.sum((y_data - y_fit)**2) / np.sum((y_data - np.mean(y_data))**2)


def phi(ref_protein, mut_protein):
	""" Makes this easier to use! """
	from phi import phi
	return phi(ref_protein, cmp_protein)






def plot_figure(equilibrium, chevron, pth=None, display=False, save=False):
	""" Plot a figure with the data and fits combined.

	TODO: Clean this up. It is terrible.
	"""
	if equilibrium.two_state:
		dfifty = equilibrium.midpoint
		dfive = equilibrium.point(0.05)
		dninetyfive = equilibrium.point(0.95)

	fity = chevron.results.y
	res = chevron.results.residuals
	fiteq = equilibrium.results.y

	#
	eq_y_range = (np.min(equilibrium.y)-.1, np.max(equilibrium.y)+.1)
	kin_x_range = (-.1, np.max(chevron.x)+.1)
	# NOTE (ergm) changed limits of the dfifity to make a prettier picture 1/9/2017
	kin_y_range = (np.min(chevron.y_raw), np.max(chevron.y_raw))


	# PLOT THE EQUILIBRIUM DATA
	plt.figure(figsize=(14,8))
	plt.subplot2grid((4,2),(0,0))
	if equilibrium.two_state:
		plt.plot([dfifty,dfifty],eq_y_range,'b-',[dfive,dfive],eq_y_range,'b:',
				[dninetyfive,dninetyfive],eq_y_range,'b:')
	plt.plot(equilibrium.x, equilibrium.y,'ko' , markersize=constants.MARKER_SIZE)
	plt.plot(equilibrium.results.x, fiteq, 'r-', linewidth=constants.LINE_WIDTH)
	plt.ylabel('Signal (A.U.)', fontsize=constants.FONT_SIZE)
	plt.title(chevron.ID, fontsize=constants.FONT_SIZE)
	plt.xlim(kin_x_range)


	# PLOT THE CHEVRON
	plt.subplot2grid((4,2),(1,0), rowspan=2)
	plt.semilogy()
	if equilibrium.two_state:
		# NOTE (ergm) changed limits of the dfifity to make a prettier picture 1/9/2017
		plt.plot([dfifty,dfifty],kin_y_range,'b-',[dfive,dfive],kin_y_range,
				'b:', [dninetyfive,dninetyfive],kin_y_range,'b:')
		if 'k2' in chevron.phases:
			plt.plot([dfifty,dfifty],[1e-2, 1.5e3],'b-',[dfive,dfive],[1e-2, 1.5e3],
					'b:', [dninetyfive,dninetyfive],[1e-2, 1.5e3],'b:')

	if chevron.components:
		for c in chevron.components:
			plt.plot(chevron.results.x, chevron.components[c], 'r--', linewidth=constants.LINE_WIDTH)

	if 'k2' in chevron.phases:
		plt.plot(chevron.denaturant['k2'], chevron.rates['k2'], 'b^', markersize=constants.MARKER_SIZE)

	plt.plot(chevron.x, chevron.y_raw, 'ko' , markersize=constants.MARKER_SIZE)
	plt.plot(chevron.results.x, fity, 'r-', linewidth=constants.LINE_WIDTH)
	plt.ylabel(r'$k_{obs} (s^{-1})$', fontsize=constants.FONT_SIZE)
	plt.xlim(kin_x_range)
	# NOTE (ergm) changed to make a prettier picture 1/9/2017
	plt.ylim(np.min(chevron.y_raw-0.5e-2), np.max(chevron.y_raw+5))

	# NOTE (ergm) changed to make a prettier picture 1/9/2017
	if 'k2' in chevron.phases:
		plt.ylim(np.min(chevron.y_raw-10), np.max(chevron.y_raw+1.5e3))

	# PLOT THE RESIDUALS
	plt.subplot2grid((4,2),(3,0))
	markerline, stemlines, baseline = plt.stem(chevron.x, res, 'k:')
	plt.setp(markerline, 'markerfacecolor', 'k')
	plt.setp(stemlines, 'color', 'k')
	plt.setp(baseline, 'color', 'k')
	plt.plot([0.,10.],[0.,0.],'k-', linewidth=constants.LINE_WIDTH)
	plt.xlabel(chevron.denaturant_label, fontsize=constants.FONT_SIZE)
	plt.ylabel(r'$k_{fit}-k_{obs} (s^{-1})$', fontsize=constants.FONT_SIZE)
	plt.xlim(kin_x_range)
	# NOTE (ergm) editted out to make a prettier picture 4/9/2017
	# plt.ylim((-0.5,0.5))

	# now plot some output
	t = u"Data-set: {0:s} \n".format(equilibrium.ID)
	t+= u"\n"
	t+= u"Equilibrium Model: {0:s} \n".format(equilibrium.fit_func)
	for e in equilibrium.results.details:
		fit_arg, fit_val, fit_err = e
		t+= u"{0:s}: {1:2.5f} \u00B1 {2:2.5f} \n".format(fit_arg, fit_val, fit_err)
	if equilibrium.two_state:
		t+= u"Folding midpoint: {0:2.2f} M\n".format(equilibrium.midpoint)
	t+= u"$R^2$: {0:2.2f} \n".format(equilibrium.results.r_squared)
	t+= u"\n"
	try:
		t+= u"Kinetic Model: {0:s} \n".format(chevron.fit_func)
	except:
		pass
	t+= u"Fit Standard Error: {0:2.2f} \n".format(chevron.results.standard_error)
	for e in chevron.results.details:
		fit_arg, fit_val, fit_err = e
		t+= u"{0:s}: {1:.2e} \u00B1 {2:.2e} \n".format(fit_arg, fit_val, fit_err)
	t+= u"$R^2$: {0:2.2f} \n".format(chevron.results.r_squared)


	ax = plt.gcf()
	ax.text(0.65, 0.95, t, horizontalalignment='left', verticalalignment='top', fontsize=constants.FONT_SIZE)
	plt.tight_layout()
	if save:
		plt.savefig(os.path.join(pth,"Fitting"+equilibrium.ID+"_{0:s}.pdf".format(chevron.fit_func)))
	if display:
		plt.show()
	plt.close()



def plot_chevron(protein, components=False,  **kwargs):
	""" plot_chevron

	Plots a chevron.
	"""

	if not isinstance(protein, Chevron):
		raise TypeError("protein is not of type Chevron")

	plt.figure(figsize=(8,5))
	plt.semilogy(protein.x, protein.y_raw, 'ko', markersize=constants.MARKER_SIZE)

	# TODO: not sure I like this:
	for phase in protein.phases[1:]:
		k_x, k_y = protein.denaturant[phase], np.exp( protein.rate(phase) )
		plt.semilogy(k_x, k_y, 'ko', markersize=constants.MARKER_SIZE)

	if protein.results:
		plt.semilogy(protein.results.x, protein.results.y, 'r-', linewidth=constants.LINE_WIDTH)

	if protein.components and components:
		x = np.linspace(0., 10., 100)
		for c in protein.components:
			plt.plot(x, protein.components[c], 'r--', linewidth=constants.LINE_WIDTH)

	plt.title('{0:s} Chevron Plot'.format(protein.ID), fontsize=constants.FONT_SIZE)

	# only need to scale y axis if the components are plotted too
	if components:
		plt.ylim([np.min(protein.y_raw)-10., np.max(protein.y_raw)+10.])

	plt.grid(False)
	plt.xlabel(protein.denaturant_label, fontsize=constants.FONT_SIZE)
	plt.ylabel(r'$\ k_{obs}$ $(s^{-1})$', fontsize=constants.FONT_SIZE)
	plt.show()


def plot_equilibrium(protein, **kwargs):
	""" plot_equilibrium

	Plots an equilibrium curve.
	"""

	if not isinstance(protein, EquilibriumDenaturationCurve):
		raise TypeError("protein is not of type EquilibriumDenaturationCurve")

	plt.figure(figsize=(8,5))
	plt.plot(protein.x, protein.y,'ko', markersize=constants.MARKER_SIZE)
	if protein.results:
		plt.plot(protein.results.x, protein.results.y, 'r-', linewidth=constants.LINE_WIDTH)
	plt.title('{0:s} Denaturation'.format(protein.ID), fontsize=constants.FONT_SIZE)
	plt.xlim([-0.1, 9])
	plt.grid(False)
	plt.xlabel(protein.denaturant_label, fontsize=constants.FONT_SIZE)
	plt.ylabel('Signal', fontsize=constants.FONT_SIZE)
	plt.show()



if __name__ == "__main__":
	test()
