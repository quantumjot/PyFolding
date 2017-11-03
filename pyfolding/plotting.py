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
import sys
import time
import inspect
import numpy as np
import scipy as sp

import matplotlib.pyplot as plt

# use a global optimisation algorithm for the Ising model fitting
# from scipy.optimize import differential_evolution, minimize, leastsq, curve_fit

# import PyFolding specific libraries
import constants
# import models
import core
# from ising import IsingDomain

__author__ = "Alan R. Lowe"
__email__ = "a.lowe@ucl.ac.uk"


"""
CORE PLOTTING FUNCTIONS
"""

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
	plt.ylabel(r'$k_{obs}-k_{fit} (s^{-1})$', fontsize=constants.FONT_SIZE)
	plt.xlim(kin_x_range)
	# NOTE (ergm) editted out to make a prettier picture 4/9/2017
	# plt.ylim((-0.5,0.5))

	# now plot some output
	t = u"Data-set: {0:s} \n".format(equilibrium.ID)
	t+= u"\n"
	t+= u"Equilibrium Model: {0:s} \n".format(equilibrium.fit_func)
	for e in equilibrium.results.details:
		fit_arg, fit_val, fit_err = e.name, e.value, e.SE
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
		fit_arg, fit_val, fit_err = e.name, e.value, e.SE
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

	if not isinstance(protein, core.Chevron):
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

	if not isinstance(protein, core.EquilibriumDenaturationCurve):
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




if __name__ == '__main__':
	pass
