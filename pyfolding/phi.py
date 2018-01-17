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

import core
import constants

__author__ = "Alan R. Lowe"
__email__ = "a.lowe@ucl.ac.uk"




def phi(ref_protein, mut_protein):
    """ Calculate the Phi value corresponding to a particular mutation
    of the wild-type protein.


    Phi-values are calculated according to:
    Serrano L, Matouschek A, Fersht AR (1992) J Mol Biol 224:805–818.

    \Phi = \frac{RT\ln\biggl( \frac{k_f^{ref}}{k_f^{mut}} \biggr)} {\Delta\Delta G_{eq}}

    Args:
        ref_protein: A protein class of type pyfolding.Protein
        mut_protein: A protein class of type pyfolding.Protein

    Notes:




        TODO: Test this thoroughly!

    """

    if not isinstance(ref_protein, core.Protein):
        raise TypeError("Reference protein must be of type pyfolding.Protein")

    if not isinstance(mut_protein, core.Protein):
        raise TypeError("Mutant protein must be of type pyfolding.Protein")


    # Start by calculating the DDG value
    DDG = ref_protein.deltaG - mut_protein.deltaG

    # now calculate the differences in folding rate constants
    DTS = constants.RT * np.log( ref_protein.kf_H20 / mut_protein.kf_H20 )

    # now the phi value itself
    phi = DTS / DDG

    return phi



if __name__ == "__main__":
    pass
