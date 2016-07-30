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

import os
import csv
import folding


__author__ = "Alan R. Lowe"
__email__ = "a.lowe@ucl.ac.uk"

def display_fit_params(p_names, p_values, p_errs=None, fit_model=None):
	"""
	Write out the parameters to the display or to a log
	"""
	pass


def read_kinetic_data(directory, filename):
	""" Read in kinetic data in the form of an Excel worksheet. It 
	should be arranged such that each file is a different protein,
	and columns represent the following:

	[den] k1 k2 ...

	This function then returns a chevron object with the data
	"""

	protein_ID, ext = os.path.splitext(filename)
	chevron = folding.Chevron(ID=protein_ID)

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
	equilibrium = folding.EquilibriumDenaturationCurve(ID=protein_ID)

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


if __name__ == "__main__":
	pass