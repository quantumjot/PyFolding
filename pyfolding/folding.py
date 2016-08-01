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



