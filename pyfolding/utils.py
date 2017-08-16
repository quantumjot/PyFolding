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


import sys
import os
import csv

import numpy as np
import constants



def disable_autoscroll(verbose=True):
	""" Disable autoscrolling in Ipython notebooks """

	# now test to see whether we're in an Ipython environment
	if not 'ipykernel' in sys.modules:
		return

	try:
		from IPython.display import display, Javascript
	except ImportError:
		return

	disable_js = """
	IPython.OutputArea.prototype._should_scroll = function(lines) {
	    return false;
	}
	"""

	# turn off the autoscrolling
	display(Javascript(disable_js))

	# give the user some feedback
	if verbose:
		print "PyFolding: Jupyter autoscrolling has been disabled"




def write_CSV(filename, data, verbose=True):
	"""
	Write out data in a convenient CSV format for import into other
	plotting and analysis packages.

	Args:
		filename
		data
		verbose

	Notes:

		WORK IN PROGRESS

	"""

	if not isinstance(filename, basestring):
		raise IOError("Filename must be a string")

	# if not filename.endswith('.csv', '.CSV'):
	# 	filename+='.csv'

	# data should be a dictionary
	csv_results = {k:iter(data[k]) for k in data.keys()}
	n_entries = len(data.values()[0])

	with open(filename, 'wb') as csvfile:
		print "Writing .csv file ({0:s})...".format(filename)
		r = csv.DictWriter(csvfile, fieldnames=data.keys(), dialect=csv.excel_tab,
				delimiter=',')
		r.writeheader()

		for i in xrange(n_entries):
			h = {k: csv_results[k].next() for k in csv_results.keys()}
			r.writerow(h)



if __name__ == "__main__":
  pass
