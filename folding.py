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

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy import optimize

import models

__author__ = "Alan R. Lowe"
__email__ = "a.lowe@ucl.ac.uk"



class Protein(object):
	""" Protein wrapper object
	"""
	def __init__(self, ID=None):
		self.ID = ID
		self.chevron = None
		self.equilibrium = None



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
			raise Exception("Fit function must be callable")

	def fit(self, p0=None, data=None):
		""" Fit the data to the defined model. Use p0 to 
		introduce the estimated start values.
		"""
		if self.__fit_func:

			# set the default fitting parameters
			if not p0: p0 = self.__fit_func.default_params

			if isinstance(self, Chevron):
				if data == 'unfolding':
					x,y = self.unfolding_limb()
				elif data == 'folding':
					x,y = self.refolding_limb()
				else:
					x,y = self.x, self.y
				y = np.log(y)
			else:
				x,y = self.x, self.y

			# perform the actual fitting
			out = optimize.curve_fit(self.__fit_func, x, y, p0=p0, maxfev=20000)

			self.fit_params = np.array( self.__fit_func.get_fit_params(x, *list(out[0])) )
			self.fit_covar = out[1]

			if hasattr(self.__fit_func, "components"):
				self.components = self.__fit_func.components(np.linspace(0.,10.,100), *list(self.fit_params))
			
		else:
			raise Exception("Fit function must be defined")

		print "--------------------"
		print " Fitting results "
		print "--------------------"
		print "Model: {0:s}".format(self.__fit_func.name)
		for e in self.errors:
			fit_arg, fit_val, fit_err = e
			print u"{0:s}: {1:2.5f} \u00B1 {2:2.5f}".format(fit_arg, fit_val, fit_err)
		print "--------------------"



	@property 
	def fitted(self):
		if isinstance(self.fit_params, np.ndarray):
			x = np.linspace(0.,10.,100)
			return np.copy( self.__fit_func.fit_func(x, *list(self.fit_params)) )

	@property 
	def residuals(self):
		if isinstance(self.fit_params, np.ndarray):
			if isinstance(self, Chevron):
				residuals = np.log(self.y) - self.__fit_func(self.x, *list(self.fit_params))
			else:
				residuals = self.y - self.__fit_func(self.x, *list(self.fit_params))
			return residuals
		else:
			return np.zeros(self.x.shape)

	@property 
	def error(self):
		""" Calculates the standard error based on the residuals of the fit.

		$S.E. = \frac{\sigma}{\sqrt{N}}$

		"""
		if isinstance(self.fit_params, np.ndarray):
			res = self.residuals
			return np.std(res) / np.sqrt(1.*len(res))

	@property
	def r_squared(self):
		""" Calculates the r-squared value of the fit.
		"""
		if isinstance(self.fit_params, np.ndarray):
			res = self.residuals
			

			if isinstance(self, Chevron): 
				y_bar = np.mean(np.log(self.y))
				y = np.log(self.y)
			else:
				y = self.y
				y_bar = np.mean(y)
			f_i = self.__fit_func(self.x, *list(self.fit_params))

			SS_tot = np.sum((y - y_bar)**2)
			SS_res = np.sum((y - f_i)**2)

			return 1.- SS_res / SS_tot


	@property
	def errors(self):
		""" Calculates the SEM of the fitted parameters, using
		the (hopefully well formed) covariance matrix from 
		the curve fitting.
		"""
		errors = []
		for p in xrange(len(self.fit_params)):
			fit_arg = self.__fit_func.fit_func_args[p]
			fit_val = self.fit_params[p]
			fit_variance = self.fit_covar[p][p]
			fit_err = np.sqrt(fit_variance) / np.sqrt(1.*len(self.x))
			errors.append([fit_arg, fit_val, fit_err])
		return errors


	def print_fit_params(self):
		if isinstance(self.fit_params, np.ndarray):
			print self.fit_params



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
	def y(self): return np.array(self.rates[self.phases[0]])


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
		for d,r in zip(self.denaturant[phase], self.rates[phase]):
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
		for d,r in zip(self.denaturant[phase], self.rates[phase]):
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
		return self.denaturant[phase], self.rates[phase]







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
	def normalised(self):
		""" TODO: Return a normalised equilbrium curve.
		"""
		if self.__fit_func:
			return
		else:
			return

	@property 
	def m_value(self):
		if isinstance(self.fit_params, np.ndarray):
			return self.fit_params[4]
		return None

	@property 
	def midpoint(self):
		if isinstance(self.fit_params, np.ndarray):
			return self.fit_params[5]
		else:
			return None

	def point(self, fraction_folded=0.5):
		""" Return the denaturant concentration for a particular
		fraction folded. Assumes a two-state transition since I
		had to derive this equation by hand.
		"""
		if self.m_value and self.midpoint:
			if fraction_folded<0. or fraction_folded>1.:
				raise Exception("Fraction folded must be in the range 0.<x<1.")
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






def read_kinetic_data(directory, filename):
	""" Read in kinetic data in the form of an Excel worksheet. It 
	should be arranged such that each file is a different protein,
	and columns represent the following:

	[den] k1 k2 ...

	This function then returns a chevron object with the data
	"""

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


def read_equilibrium_data(directory, filename):
	""" Read in an equilbrium denaturation curve from an Excel
	worksheet. It should be arranged such that each file is a 
	different protein, and columns represent the following: 

	[den] unfolding

	This function then returns an equilbrium curve object.
	"""
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




def plot_figure(equilibrium, chevron, pth, show=False, save=True):
	""" Plot a figure with the data and fits combined.
	"""

	dfifty = equilibrium.midpoint
	dfive = equilibrium.point(0.05)
	dninetyfive = equilibrium.point(0.95)

	fity = chevron.fitted
	res = chevron.residuals
	fiteq = equilibrium.fitted

	plt.figure()
	plt.subplot2grid((4,2),(0,0))
	plt.plot([dfifty,dfifty],[-.1,1.1],'b-',[dfive,dfive],[-.1,1.1],'b:', [dninetyfive,dninetyfive],[-.1,1.1],'b:')
	plt.plot(np.linspace(0., 10., 100), fiteq, 'k-', linewidth=2)
	plt.plot(equilibrium.x, equilibrium.y,'wo')
	plt.ylim([-.1,1.1])
	plt.ylabel('Signal (A.U.)')
	plt.title(chevron.ID)
	plt.subplot2grid((4,2),(1,0), rowspan=2)
	plt.semilogy()
	plt.plot([dfifty,dfifty],[0.01,100.],'b-',[dfive,dfive],[0.01,100.],'b:', [dninetyfive,dninetyfive],[0.01,100.],'b:')
	plt.plot(np.linspace(0., 10., 100), fity, 'k-', linewidth=2)
	if chevron.components:
		x = np.linspace(0., 10., 100)
		for c in chevron.components:
			plt.plot(x, chevron.components[c], 'r--', linewidth=2)

	if 'k2' in chevron.phases:
		plt.plot(chevron.denaturant['k2'], chevron.rates['k2'], 'w^')
	plt.plot(chevron.x, chevron.y, 'wo')		
	plt.ylabel(r'$k_{obs} (s^{-1})$')
	plt.xlim((0,10))
	plt.ylim((0.01,100.))
	plt.subplot2grid((4,2),(3,0))
	markerline, stemlines, baseline = plt.stem(chevron.x, res,'k:')
	plt.setp(markerline, 'markerfacecolor', 'w')
	plt.setp(stemlines, 'color', 'k')
	plt.setp(baseline, 'color', 'k')
	plt.plot([0.,10.],[0.,0.],'k-', linewidth=2)
	plt.xlabel(chevron.denaturant_label)
	plt.ylabel(r'$k_{fit}-k_{obs} (s^{-1})$')
	plt.xlim((0,10))
	plt.ylim((-0.5,0.5))
	# now plot some output
	t = u"Data-set: {0:s} \n".format(equilibrium.ID) 
	t+= u"Model: {0:s} \n".format(chevron.fit_func)
	t+= u"Folding midpoint {0:2.2f} M\n".format(equilibrium.midpoint)
	t+= u"Fit Standard Error: {0:2.2f} \n".format(chevron.error)
	t+= u"Fitting parameters: \n"
	for e in chevron.errors:
		fit_arg, fit_val, fit_err = e
		t+= u"{0:s}: {1:.2e} \u00B1 {2:.2e} \n".format(fit_arg, fit_val, fit_err)
	t+= u"$R^2$: {0:2.2f} \n".format(chevron.r_squared)

	ax = plt.gcf()
	ax.text(0.65, 0.95, t, horizontalalignment='left', verticalalignment='top', fontsize=10.)
	plt.tight_layout()
	if save: 
		plt.savefig(os.path.join(pth,"Fitting"+equilibrium.ID+"_{0:s}.pdf".format(chevron.fit_func)))
	if show: 
		plt.show()
	plt.close()



if __name__ == "__main__":
	pass


	

	

