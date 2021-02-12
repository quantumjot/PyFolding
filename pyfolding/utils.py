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

from . import core
import numpy as np
from . import constants

from collections import OrderedDict



class __Temperature(object):
    """ Maintain temperature information across all functions. Meant to be a
    wrapper for a global variable, temperature which is used across different
    models.

    Args:
        temperature: set the temperature in celsius

    Properties:
        temperature: return the temperature in celsius
        RT:    return the product of the ideal gas constant and the temperature

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
            raise ValueError("Temperature ({0:2.2f}) is not in valid range "
                "(0-100 degrees C)".format(value))
        self.__temperature = value

    @property
    def RT(self):
        return constants.IDEAL_GAS_CONSTANT_KCAL * (constants.ZERO_KELVIN + self.temperature)




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
        print("PyFolding: Jupyter autoscrolling has been disabled")



def check_filename(directory, filename):
    """ Check the filename for consistency.
    """

    if not isinstance(directory, str) or not isinstance(filename, str):
        raise TypeError('Pyfolding expects a filename as a string')

    if not filename.lower().endswith(('.csv', '.CSV')):
        raise IOError('PyFolding expects a .CSV file as input: {0:s}'.format(filename))

    if not os.path.exists(os.path.join(directory, filename)):
        raise IOError('PyFolding could not find the file: {0:s}'.format(os.path.join(directory, filename)))




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

    if not isinstance(filename, str):
        raise IOError("Filename must be a string")

    # data should be a dictionary
    csv_results = {k:iter(data[k]) for k in list(data.keys())}
    n_entries = len(list(data.values())[0])

    with open(filename, 'wb') as csvfile:
        if verbose:
            print("Writing .csv file ({0:s})...".format(filename))
        r = csv.DictWriter(csvfile, fieldnames=list(data.keys()), dialect=csv.excel_tab,
                delimiter=',')
        r.writeheader()

        for i in range(n_entries):
            h = {k: next(csv_results[k]) for k in list(csv_results.keys())}
            r.writerow(h)



class DataImporter(object):
    """
    Generic data importer class.

    This class will read in .CSV files containing data of the format:

    x    y_0     y_1 ...
    0    1.0        5.6 ...
    1    2.0        7.4 ...
    ...

    And output a data class of the type specified. This can be one of
        - EquilibriumDenaturationCurve
        - Chevron
        - Generic (this is meant to be a catch all for other, future types)

    Properties:
        type

    Members:
        load

    Notes:
        None

    """

    def __init__(self, datatype='GenericData'):
        self.type=datatype

    @property
    def type(self): return self.__type
    @type.setter
    def type(self, datatype):
        if not isinstance(datatype, str):
            raise TypeError('Data type must be specified as a string')
        if datatype not in ('EquilibriumDenaturationCurve','Chevron','GenericData'):
            raise ValueError('Data type {s} is not recognised'.format(datatype))

        self.__type = datatype

    def load(self, fullfilename):
        """ Load the data and instantiate the correct data model """

        directory, filename = os.path.split(fullfilename)

        # check the filename
        check_filename(directory, filename)

        protein_ID, ext = os.path.splitext(filename)

        # make a new data object of the specified type
        DataObject = getattr(core, self.type)
        data = DataObject(ID=protein_ID)

        # open the data and set up a new object
        with open(fullfilename, 'rU') as data_file:
            data_reader = csv.reader(data_file, delimiter=',', quotechar='|')
            header = next(data_reader)

            data.labels = header
            x_label = header[0]
            y_labels = header[1:]
            data.data = {'x':{y_lbl:[] for y_lbl in y_labels},
                        'y':{y_lbl:[] for y_lbl in y_labels}}

            # now group all of the data and insert into object
            for row in data_reader:
                x_val = float( row.pop(0) )
                for phase, y_val in zip(y_labels, row):
                    if y_val:
                        data.data['x'][phase].append(x_val)
                        data.data['y'][phase].append(float(y_val))

        return data




class FitExporter(object):
    """ A class to export Fit data from FitResult objects """
    def __init__(self):
        self.verbose = False

    def export(self, filename, results):
        """ Take a list of results or a single result and ouput files """
        if isinstance(results, list):
            for r in results:
                self.export(filename, r)
            return

        if not isinstance(results, core.FitResult):
            raise TypeError('FitExporter requires a FitResult object as input')

        # save out the data
        filepath, fn = os.path.split(filename)
        fn, ext = os.path.splitext(fn)

        # first save out the curves
        save_file = os.path.join(filepath, fn+'_'+results.ID+'_FITCURVE'+ext)
        data = OrderedDict([('x', results.x_fit), ('y',results.y_fit)])
        write_CSV(save_file, data, verbose=self.verbose)

        # now save out the fit parameters
        save_file = os.path.join(filepath, fn+'_'+results.ID+'_FITRESULT'+ext)
        data = OrderedDict([('Parameter', [f.name for f in results.fit_params]),
                            ('Type', [f.type for f in results.fit_params]),
                            ('Value',[f.value for f in results.fit_params]),
                            ('Error',[f.SE for f in results.fit_params])])
        write_CSV(save_file, data, verbose=self.verbose)




if __name__ == "__main__":
  pass
