#!/usr/bin/env python

"""
Python implementation of common model fitting to protein folding
data. Simply automates some fitting and value calculation with
common objects. Will be extended to include phi-value analysis 
and other common calculations.

Allows for quick model evaluation and plotting.

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

import fit

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

	@property
	def fit_func(self): return self.__fit_func.__name__
	@fit_func.setter
	def fit_func(self, fit_func=None):
		if hasattr(fit_func, "__call__"):
			self.__fit_func = fit_func
		else:
			raise Exception("Fit function must be callable")

	def fit(self, p0=None):
		""" Fit the data to the defined model. Use p0 to 
		introduce the estimated start values.
		"""
		if self.__fit_func:
			if isinstance(self, Chevron):
				# TODO: sort this out!
				objective_func = lambda p, x, y: np.ravel( y - np.log(self.__fit_func(x,p)) )
				y = np.log(self.y)
			else:
				objective_func = lambda p, x, y: np.ravel( y - self.__fit_func(x,p) )
				y = self.y

			out = optimize.leastsq(objective_func, p0, args=(self.x, y), maxfev=2000, xtol=1e-8, full_output=True)
			self.fit_covar = out[1]	
			self.fit_params = out[0]
		else:
			raise Exception("Fit function must be defined")

	@property 
	def fitted(self, x=None):
		if isinstance(self.fit_params, np.ndarray):
			if not x: x = np.linspace(0.,10.,100)
			return self.__fit_func(x, self.fit_params)

	@property 
	def residuals(self):
		if isinstance(self.fit_params, np.ndarray):
			if isinstance(self, Chevron):
				return np.log(self.y) - np.log(self.__fit_func(self.x, self.fit_params))
			else:
				return self.y - self.__fit_func(self.x, self.fit_params)

	@property 
	def error(self):
		""" Calculates the standard error based on the residuals of the fit.

		$S.E. = \frac{\sigma}{\sqrt{N}}$

		TODO: use the covariance matrix from the fit to estimate 
		fit parameter errors properly

		Notes:
			http://stackoverflow.com/questions/14581358/getting-standard-errors-on-fitted-parameters-using-the-optimize-leastsq-method-i

		"""
		if isinstance(self.fit_params, np.ndarray):
			res = self.residuals
			return np.std(res) / 1.*len(res)



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
				den.append(d)
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
				den.append(d)
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
		return None

	@property 
	def deltaG(self):
		""" Return the deltaG value based on the fit of the data """
		if self.m_value and self.midpoint:
			return self.m_value * self.midpoint
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

		

def plot_figure(equilibrium, chevron, pth):
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
	plt.ylabel(r'$F_{norm} (A.U.)$')
	plt.ylim([-.1,1.1])
	plt.title(chevron.ID)
	plt.subplot2grid((4,2),(1,0), rowspan=2)
	plt.semilogy()
	plt.plot([dfifty,dfifty],[0.01,100.],'b-',[dfive,dfive],[0.01,100.],'b:', [dninetyfive,dninetyfive],[0.01,100.],'b:')
	plt.plot(np.linspace(0., 10., 100), fity, 'k-', linewidth=2)
	plt.plot(chevron.x, chevron.y, 'wo')
	if 'k2' in chevron.phases:
		plt.plot(chevron.denaturant['k2'], chevron.rates['k2'], 'w^')
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
	plt.tight_layout()
	plt.savefig(os.path.join(pth,"Fitting"+equilibrium.ID+".pdf"))
	plt.show()








if __name__ == "__main__":

	all_proteins = ['WT','V10A','A14G','I26A','A43G','A47G','F58I','I79V','A80G','A81G','V89A','A91G','A113G','A122G','A147G','I155V','I157V','A179G','A188G','V192A','L209A','V211A']
	pth = "/Users/ubcg83a/Dropbox/Code/Folding/"
	fn = "WT.csv"

	start_2 = np.array([50., -1.3480, 5e-4, 1.])
	start_3 = np.array([4.5e-4, 1.1, 1.3e9, -6.9, -9.5e-1, 1.4e-8, -1.6])
	start_eq = np.array([1., 0.1, 0.0, 0.1, 1.5, 5.])
	start_oliveberg = np.array([50., -1.3480, 5e-4, 1., .1])

	# load the data
	chevron = read_kinetic_data(os.path.join(pth,"Chevrons"), fn)
	equilibrium = read_equilibrium_data(os.path.join(pth,"Equilibrium"), fn)

	# fit the data
	chevron.fit_func = fit.three_state_chevron
	chevron.fit(start_3)
	equilibrium.fit_func = fit.two_state_equilibrium
	equilibrium.fit(start_eq)

	# plot it and save a figure
	plot_figure(equilibrium, chevron, pth)

	

