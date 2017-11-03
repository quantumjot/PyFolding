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
from ising import IsingDomain

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


"""
ISING PLOTTING FUNCTIONS
"""

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


	if not isinstance(topologies, list):
		raise TypeError("Topologies must be a list of topologies.")

	# first we need to check whether the topologies are instantiated classes or
	# not. If they are already instantiated, then we need to recover the
	# underlying classes...
	if isinstance(topologies[0][0], IsingDomain):
		tmp_topologies = []
		for t in topologies:
			new_topology = [d.__class__ for d in t]
			tmp_topologies.append(new_topology)
	else:
		tmp_topologies = topologies

	from matplotlib.patches import Patch

	if not labels:
		labels = ['Protein {0:d}'.format(i) for i in xrange(len(tmp_topologies))]

	if 'fold' in kwargs:
		raise DeprecationWarning('Fold keyword is being replaced with collapse.')
		collapse = kwargs['fold']


	# collapse the topologies
	compact = []
	domain_types = set()

	for l,t in zip(labels, tmp_topologies):

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
	plt.figure(figsize=(14, len(tmp_topologies)*1.5))
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



	ax.set_yticks(np.arange(len(tmp_topologies)), minor=False)
	ax.set_xticklabels([], minor=False, rotation='vertical')
	ax.set_yticklabels(labels, minor=False)

	# make the legend entries
	l = [ Patch(edgecolor=d_color(d), facecolor='w', label=d) for d in domain_types ]


	ax.set_xlim([-1.,11.])
	ax.set_ylim([-1., len(tmp_topologies)])
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


if __name__ == '__main__':
	pass
