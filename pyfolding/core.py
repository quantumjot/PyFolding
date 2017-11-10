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
from collections import OrderedDict
import numpy as np
from scipy import optimize
from scipy.stats import t as t_distrb

# pyfolding imports
import utils
import constants
from plotting import *

__author__ = "Alan R. Lowe"
__email__ = "a.lowe@ucl.ac.uk"
__version__ = constants.VERSION

# by default turn off autoscrolling if it exists
utils.disable_autoscroll()

# set up a global temperature object
temperature = utils.__Temperature()




"""
===========================================================
FILE I/O OPERATIONS
===========================================================
"""

def read_kinetic_data(directory=None, filename=None):
	""" Read in kinetic data in the form of an .csv worksheet. It
	should be arranged such that each file is a different protein,
	and columns represent the following:

	[den] k1 k2 ...

	This function then returns a chevron object with the data
	"""
	reader = utils.DataImporter(datatype='Chevron')
	return reader.load(os.path.join(directory,filename))



def read_equilibrium_data(directory=None, filename=None):
	""" Read in an equilbrium denaturation curve from a .csv
	worksheet. It should be arranged such that each file is a
	different protein, and columns represent the following:

	[den] unfolding

	This function then returns an equilbrium curve object.
	"""
	reader = utils.DataImporter(datatype='EquilibriumDenaturationCurve')
	return reader.load(os.path.join(directory,filename))


def read_generic_data(directory=None, filename=None):
	""" Read in a generic dataset from a .csv
	worksheet. It should be arranged such that each file is a
	different protein, and columns represent the following:

	x y_0 y_1 ....

	This function then returns a generic data object.
	"""
	reader = utils.DataImporter()
	return reader.load(os.path.join(directory,filename))




"""
===========================================================
SETTING CALCULATION TEMPERATURE
===========================================================
"""

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
BASE CLASSES
===========================================================
"""

class DataTemplate(object):
	""" DataTemplate

	Base class fo chevrons, equilibrium denaturation curves and generic data.
	Takes care of common functions such as fitting of models.


	Data is stored internally as a dictionary and associated list, for example:

	labels = ['denaturant', 'k1', 'k2']
	data = {'x': {'k1': [], 'k2': []}, 'y': {'k1':[], 'k2':[]}

	Subclassed objects may use these data in different ways, for example as
	Chevron plots or Equilibrium denaturation curves.


	Usage:

		>> data['k1']

		returns a tuple of x['k1'] and y['k1']



	Properties:
		datasets - return a list of datasets in the model
		fit_func - return/set the fit function for the dataset
		fit_func_args - return the fit function arguments
		fit_params - return the final parameters following the fit
		results - a FitResult object following fitting


	Members:

	Notes:

	"""
	def __init__(self):
		# store the raw data in a dictionary
		self.labels = []
		self.data = {}

		# store associated fit functions
		self.__fit_func = None
		self.__fit = None
		self.__fit_residuals = None
		self.components = None

	def initialise(self):
		raise NotImplementedError

	def __getitem__(self, dataset):
		""" Return an XY pair from the dataset, based on the label """
		if not isinstance(dataset, basestring):
			raise TypeError('Dataset must be specified as as string')
		if dataset not in self.datasets:
			raise ValueError('Dataset {0:s} not found'.format(dataset))

		return ( np.array(self.data['x'][dataset], dtype='float'),
				np.array(self.data['y'][dataset], dtype='float') )

	@property
	def datasets(self): return self.labels[1:]

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
	def fit_params(self):
		return [p.value for p in self.__fit.fit_params]

	@property
	def results(self):
		return self.__fit
	@results.setter
	def results(self, result):
		if not isinstance(result, FitResult):
			raise TypeError("Results must be of type FitResult")

		print "Warning: overwriting fit result for {0:s}".format(self)
		self.__fit = result

	def fit(self, p0=None, const=None):
		""" Fit the data to the defined model. Use p0 to introduce the estimated
		 start values.
		"""
		if self.__fit_func:

			# reset components
			self.components = None

			# set the default fitting parameters
			if not p0: p0 = self.__fit_func.default_params

			# set up the fit
			f = GlobalFit()
			f.fit_funcs = [self.__fit_func]
			if const: f.constants = [const]
			f.shared = [] # no shared parameters by default
			f.x = [self.x]
			f.y = [self.y]
			f.ID = [self.ID]

			out, covar = f.fit( p0=p0 )
			self.__fit = f.results[0]

			if hasattr(self.__fit_func, "components"):
				self.components = self.__fit_func.components(constants.XSIM, *out.tolist())

		else:
			raise AttributeError("Fit function must be defined first.")


		self.__fit.display()


	@property
	def fitted_x(self):
		raise DeprecationWarning("This feature will be deprecated soon.")

	@property
	def fitted(self):
		raise DeprecationWarning("This feature will be deprecated soon.")

	def print_fit_params(self):
		raise DeprecationWarning("This feature will be deprecated soon.")
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
		""" Export the fit. """
		exporter = utils.FitExporter()
		exporter.export(filename, self.results)




class Protein(object):
	""" Protein wrapper object.

	This class wraps different types of data and acts as a container object for
	a single protein.  It can contain equilbrium, kinetic and other types of
	data. The object can be passed to higher-order functions, such as 'phi' that
	use multiple datasets for calculations.

	Properties:
		deltaG - the equilbrium deltaG value from equilbrium data
		kf_H20 - the observed folding rate in water

	Notes:
		None
	"""
	def __init__(self, ID=None):
		self.ID = ID
		self.chevron = None
		self.equilibrium = None
		self.other = None

	@property
	def deltaG(self): return self.equilibrium.deltaG

	@property
	def kf_H20(self): return self.chevron.results.y_fit[0]




class GenericData(DataTemplate):
	""" A generic data model.
	"""
	def __init__(self, ID=None):
		DataTemplate.__init__(self)
		self.ID = ID

	@property
	def x(self): return self[self.datasets[0]][0]
	@property
	def y(self): return self[self.datasets[0]][1]

	@property
	def y_raw(self): return self.y

	def initialise(self):
		pass



class Chevron(DataTemplate):
	""" Chevron plot for protein folding kinetics.

	Args:

	Methods:

	Notes:

	"""
	def __init__(self, ID=None):
		DataTemplate.__init__(self)
		self.ID = ID
		self.__midpoint = None


	@property
	def denaturant_label(self): return self.labels[0]
	@property
	def phases(self): return self.datasets
	@property
	def rates(self): return {k:self[k][1] for k in self.phases}
	@property
	def denaturant(self): return {k:self[k][0] for k in self.phases}




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







class EquilibriumDenaturationCurve(DataTemplate):
	""" Equilibrium Denaturation curve
	Args:

	Methods:

	Notes:

	"""
	def __init__(self, ID=None):
		DataTemplate.__init__(self)
		self.ID = ID

	@property
	def denaturant_label(self): return self.labels[0]
	@property
	def curves(self): return self.datasets
	@property
	def signal(self): return {k:self[k][1] for k in self.curves}
	@property
	def denaturant(self): return {k:self[k][0] for k in self.curves}


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
		if isinstance(self.fit_params, list):
			return self.fit_params[ self.fit_func_args.index('m') ]
		return None


	@property
	def midpoint(self):
		if isinstance(self.fit_params, list):
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






"""
===========================================================
MODEL FITTING FUNCTIONS
===========================================================
"""


def FIT_ERROR(x):
	""" Return a generic fit error """
	if isinstance(x, np.ndarray):
		return np.ones(x.shape)*constants.FITTING_PENALTY
	else:
		return None



class FitParameter(object):
	""" Object to store parameter error information """
	def __init__(self, name, value, param_type='free'):
		self.name = name
		self.value = value
		self.type = param_type
		self.DoF = None
		self.SE = 0
		self.CI = [-np.inf, np.inf]
		self.covar = None
		self.r_squared = None

	@property
	def name(self): return self.__name
	@name.setter
	def name(self, arg_name):
		if not isinstance(arg_name, basestring):
			raise TypeError('Arg name must be of type string')
		self.__name = arg_name

	@property
	def type(self): return self.__type
	@type.setter
	def type(self, arg_type):
		if not isinstance(arg_type, basestring):
			raise TypeError('Arg type must be of type string')
		if arg_type not in ['free', 'shared', 'constant']:
			raise ValueError('Arg type must be either free, shared or constant')
		self.__type = arg_type

	@property
	def CI_low(self): return self.CI[0]
	@property
	def CI_high(self): return self.CI[1]



class GlobalFit(object):
	""" GlobalFit

	Wrapper function to perform global fitting.  This acts as a wrapper for
	multiple FitModels, enabling the user to pair datasets and models and share
	data or arguments.

	For each fit function, a list of arguments is compiled. Those belonging to
	the shared or constant type are set respectively.

	Note that a single or individual fit is just a special case of a global fit
	where there are no shared values and only one dataset. This wrapped can be
	used for that purpose too...

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
		self.ID = []
		self.x = []
		self.y = []
		self.__fit_funcs = []
		self.__shared = []
		self.__initialised = False
		self.__params = None
		self.__results = None

		self.covar = None

	@property
	def fit_funcs(self): return self.__fit_funcs
	@fit_funcs.setter
	def fit_funcs(self, fit_funcs):
		for fit_func in fit_funcs:
			if not hasattr(fit_func, "__call__"): continue
			# append it and instantiate it
			if isinstance(fit_func, FitModel):
				self.__fit_funcs.append(fit_func)
			else:
				self.__fit_funcs.append( fit_func() )

	@property
	def constants(self):
		return [f.constants for f in self.__fit_funcs]
	@constants.setter
	def constants(self, const=None):
		if len(const) != len(self.__fit_funcs):
			raise ValueError("Number of constants should be the same as number of fit functions")

		for constant, fit_func in zip(const, self.__fit_funcs):
			fit_func.constants = constant


	@property
	def shared(self):
		return self.__shared
	@shared.setter
	def shared(self, shared_args=[]):
		""" Set the shared arguments for the global fit """
		if not isinstance(shared_args, (list, tuple)):
			raise TypeError('Shared args must be of type list or tuple')
		if not all([isinstance(a, basestring) for a in shared_args]):
			raise TypeError('Shared args must be a list of strings.')

		# TODO(arl): check that these shared params exist in the fit functions
		# and report an error if incorrect...

		self.__shared = list(set(shared_args))


	@property
	def params(self): return self.__params


	def __call__(self, *args):
		""" Dummy call for all fit functions """
		if not self.__initialised: self.initialise()

		x = args[0]
		fit_args = args[1:]

		# now set the values of the objects
		for p, p_val in zip(self.params, fit_args):
			self.__params[p].value = p_val

		ret = np.array(())
		for i, fit_func in enumerate(self.fit_funcs):
			ret = np.append( ret, self.eval_func(i) )
		return ret

	def initialise(self):
		""" Set up all of the shared, constant and free parameters """

		if len(self.ID) != len(self.x):
			self.ID = ['protein_{0:d}'.format(i) for i in xrange(len(self.x))]

		shared = {s:FitParameter(s, 0.0, param_type='shared') for s in self.shared}

		# set up an ordered dictionary of the parameter objects
		all_params = OrderedDict(shared)

		for f in self.fit_funcs:
			fit_func_params = []
			const = [c[0] for c in f.constants]
			for arg in f.fit_func_args:
				if arg in shared:
					fit_func_params.append(shared[arg])
				elif arg in const:
					c_val = f.constants[const.index(arg)][1]
					fit_func_params.append(FitParameter(arg, c_val, param_type='constant'))
				else:
					fit_func_params.append(FitParameter(arg, 0.0, param_type='free'))

			f.rick_and_morty = fit_func_params
			# print f.name, [(g.name, g.type) for g in f.rick_and_morty]

		# now make the master list of params
		for i, f in enumerate(self.fit_funcs):
			for p in f.rick_and_morty:
				if p.type=='shared' and p.name not in all_params:
					all_params[p.name] = p
				elif p.type not in ('shared','constant'):
					# all_params[p.name+'_'+str(i)] = p
					all_params[p.name+'_{'+self.ID[i]+'}'] = p

		# save this ordered dict for later
		self.__params = all_params

		# set the flag so that we don't do this again
		self.__initialised = True


	def eval_func(self, i):
		""" Evaluate the fit function """
		if i<0 or i>len(self.fit_funcs):
			raise ValueError('Cannot evaluate fit function {0:d}'.format(i))
		fit_func = self.fit_funcs[i]
		x_this = np.array( self.x[i] )
		args_this = [a.value for a in fit_func.rick_and_morty]
		return fit_func(x_this, *args_this)

	def fit(self, p0=[], bounds=None):
		""" Run the fit. """
		x = np.concatenate([x for x in self.x])
		y = np.concatenate([y for y in self.y])

		# fit the data
		if bounds:
			out, covar = optimize.curve_fit(self, x, y, p0=p0, bounds=bounds, max_nfev=20000)
		else:
			out, covar = optimize.curve_fit(self, x, y, p0=p0, maxfev=20000)


		# now finalise and set up the results
		self.all_residuals = residuals(y_data=y, y_fit=self(x, *out))
		self.finalise(out, covar)

		return out, covar


	def finalise(self, out, covar):
		""" Take the results of the fitting, set the parameter values and
		calculate errors.
		"""

		# put the parameter values in
		for i, p in enumerate(self.params):
			self.params[p].value = out[i]
			self.params[p].covar = covar[i,i]

		self.covar = covar
		self.__results = []

		# set up the fit result objects
		for i,f in enumerate(self.fit_funcs):
			result = FitResult(fit_name=f.name, fit_params=f.rick_and_morty)
			result.ID = self.ID[i]
			result.method = "pyfolding.GlobalFit and scipy.optimize.curve_fit"
			result.y = self.eval_func(i)
			result.x_fit = constants.XSIM
			result.y_fit = f(constants.XSIM, *[a.value for a in f.rick_and_morty])
			result.covar = covar
			result.residuals = residuals(y_data=self.y[i], y_fit=result.y)
			result.r_squared = r_squared(y_data=self.y[i], y_fit=result.y)
			result.all_residuals = self.all_residuals

			self.__results.append(result)


	@property
	def results(self):
		return self.__results




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
		all_residuals: all residuals for a global fit (same as residuals if an
						individual fit)
		r_squared: r^2 value for the fit

	Members:
		display: return a formatted output of the fitting statistics

	Notes:
		TODO(arl): implement an export function

	"""

	def __init__(self, fit_name=None, fit_params=None):
		self.ID = None
		self.fit_params = fit_params
		self.name = fit_name

		self.covar = None
		self.residuals = None
		self.all_residuals = None
		self.r_squared = None

		self.x_fit = None
		self.y_fit = None

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

		table_width = max([len("Model: "+self.name), len(" Fitting results "), 80])

		nl = 0
		for p in self.details:
			nl = max(nl, len(p.name))


		print u"="*table_width
		print u"Fitting results"
		print u"="*table_width
		if self.ID: print u"ID: {0:s}".format(self.ID)
		print u"Model: {0:s}".format(self.name)
		print u"Optimiser: {0:s}".format(self.__method)
		print u"Temperature: {0:2.2f}\u00B0C\n".format(temperature.temperature)

		for p in self.details:
			self.display_row(p, nl)


		print u"-"*table_width
		print u"R^2: {0:2.5f}".format(self.r_squared)
		print u"="*table_width
		print "\n"

	def display_row(self, p, max_name_len):
		""" Take a parameter and display a row of the table """

		p_name = p.name.ljust(max_name_len)

		if p.type is 'constant':
			print u"({0:s}) {1:s} {2:>10.5f}".format(p.type[0], p_name, p.value)
			return

		print u"({0:s}) {1:s} {2:>10.5f} \u00B1 {3:<10.5f}" \
			u" \t 95\u0025 CI[{4:>10.5f}, {5:>10.5f}]".format(p.type[0], p_name, p.value, p.SE, p.CI_low, p.CI_high)


	def confidence(self, i):
		""" Return the 95 per cent confidence interval for a fitted parameter
		https://stats.stackexchange.com/questions/72047/when-fitting-a-curve-how-do-i-calculate-the-95-confidence-interval-for-my-fitt
		[BestFit(Pi) +/- t(95%,DF)*SE(Pi)

		NOTES:
			TODO(arl): make this a user defined interval
		"""
		ci = constants.CONFIDENCE_INTERVAL / 100.0
		conf = t_distrb.pdf(0.95, self.DoF) * self.SE(i)
		return (self.fit_params[i].value-conf, self.fit_params[i].value+conf)


	def SE(self, i):
		""" Return the SE for parameter i
		SE(Pi) = sqrt[ (SS/DF) * Cov(i,i) ]
		"""
		SE = np.sqrt( (np.sum(self.all_residuals**2) / self.DoF) * self.fit_params[i].covar )
		# print 'SS: {0:2.5f}, DoF: {1:2.5f}, covar: {2:2.5f}, SE: {3:2.5f}'.format(np.sum(self.all_residuals**2), self.DoF, self.fit_params[i].covar, SE)
		return SE

	@property
	def DoF(self):
		""" Return the number of degrees of freedom, essentially the difference
		between the number of data points and the number of fit parameters
		"""
		return len(self.all_residuals) - len(self.fit_params)

	@property
	def details(self):
		""" Return a zipped list of the fit arguments, values and errors """
		details = []
		for i, f in enumerate(self.fit_params):
			if f.type is 'constant':
				details.append(f)
				continue
			f.DoF = self.DoF
			f.SE = self.SE(i)
			f.CI = self.confidence(i)
			# f.covar = self.covar[i,i]
			details.append(f)
		return details

	@property
	def standard_error(self):
		""" Return the standard error of the fit """
		return np.std(self.residuals) / np.sqrt(1.*len(self.residuals))

	def export(self, filename):
		raise NotImplementedError





















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
		self.constants = []

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
		# fit_args = self.get_fit_params(x, *args)
		# return self.error_func( self.fit_func(x, *fit_args) )
		return self.error_func( self.fit_func(x, *args) )


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

	# @staticmethod
	def print_equation(self):
		from IPython.display import display, Math, Latex
		display(Math(self.equation))
		return None

	def info(self):
		self.print_equation()
		print self.__doc__











def r_squared(y_data=None, y_fit=None):
	return 1. - np.sum((y_data - y_fit)**2) / np.sum((y_data - np.mean(y_data))**2)

def residuals(y_data=None, y_fit=None):
	return y_data - y_fit


def phi(ref_protein, mut_protein):
	""" Makes this easier to use! """
	from phi import phi
	return phi(ref_protein, cmp_protein)


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
	eq_sim = eq_raw + np.random.randn(100,)*0.01

	equilibrium.labels = ['[Denaturant] (M)', 'e1']
	equilibrium.data = {'x':{'e1':den}, 'y':{'e1':eq_sim}}
	equilibrium.fit_func = models.TwoStateEquilibrium


	# generate a two-state chevron curve, with Gaussian noise
	# kf, mf, ku, mu
	kin_model = models.TwoStateChevron()
	kin_raw = kin_model.fit_func(den, *truth['kin'])
	kin_sim = np.exp( np.log(kin_raw) + np.random.randn(100,)*0.001 )


	chevron.labels = ['[Denaturant] (M)', 'k1']
	chevron.data = {'x':{'k1':den}, 'y':{'k1':kin_sim}}
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


	print 'SUCCESS - Test completed!'

	# plot the output
	plot_figure(equilibrium, chevron, display=True)








if __name__ == "__main__":
	test()
