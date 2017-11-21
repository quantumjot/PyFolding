# PyFolding


PyFolding is a simple Python based framework for fitting kinetic and
thermodynamic models to protein folding data.  It includes several basic models
and is extensible to enable fitting of more interesting models.

**We have now included video tutorials in the wiki to demonstrate installation
and usage of PyFolding**:
https://github.com/quantumjot/PyFolding/wiki

---
### Reference

**PyFolding: An open-source software package for graphing, analysis and simulation  
of thermodynamic and kinetic models of protein folding**  
Lowe AR, Perez-Riba A, Itzhaki L, Main E (2017) BioRXiv  
http://dx.doi.org/10.1101/191593

---

### Current models supported

+ Homopolymer Ising equilibrium
+ Heteropolymer Ising equilibrium
+ Parallel Two-State chevron
+ Parallel Two-State Unfolding chevron
+ Three State chevron
+ Three State Dimeric Intermediate equilibrium
+ Three State Fast-Phase chevron
+ Three State Monomeric Intermediate equilibrium
+ Three State Sequential Barriers chevron
+ Two State chevron
+ Two State chevron Moving Transition State
+ Two State Dimer equilibrium
+ Two State equilibrium
+ Two State equilibrium with sloping baselines

A 'TemplateModel' is provided for adding additional models. We encourage users
to generate their own models and contribute.

---
### Example scripts and Jupyter notebooks
PyFolding can be used within IPython Notebooks (Jupyter). Several example
notebooks are provided in the notebooks folder. To install Jupyter notebooks
please read the accompanying file "SETUP.md".

---

### Requirements

Pyfolding has been tested on Python 2.7, and requires the following additional packages:
+ Numpy
+ Scipy
+ Matplotlib
+ (Optional) Jupyter

The easiest way to get these packages on a Mac or PC is to install Anaconda:
https://www.anaconda.com/download/


### Installation

You can install PyFolding by cloning the repo and running the setup script:
```sh
$ git clone https://github.com/quantumjot/PyFolding.git
$ cd PyFolding
$ python setup.py install
```

+ For more detailed installation for setup on a Mac or PC please read the accompanying SETUP.md.
+ Watch the YouTube video demonstrating how to install and run PyFolding.

### Testing the installation

Once installed, a simple script will execute the test function. This tests
whether the package has installed correctly.  The following script is found in /notebooks.

```python
# import the pyfolding package
import pyfolding
# run a self test to check everything is installed
pyfolding.test()
```

Upon executing the script, the following output should be generated:

```sh
========================================================
Fitting results
========================================================
ID: Simulated protein
Model: TwoStateEquilibrium
Method: scipy.optimize.curve_fit

m:   1.48602 ± 0.00016   95% CI[1.48598, 1.48606]
d50: 4.99714 ± 0.00005   95% CI[4.99713, 4.99715]
--------------------------------------------------------
R^2: 0.99955
========================================================

========================================================
Fitting results
========================================================
ID: Simulated protein
Model: TwoStateChevron
Method: scipy.optimize.curve_fit

kf: 100.03261 ± 0.00112    95% CI[100.03232, 100.03289]
mf:   1.00006 ± 0.00000    95% CI[1.00006, 1.00006]
ku:   0.00499 ± 0.00000    95% CI[0.00499, 0.00499]
mu:   1.00015 ± 0.00000    95% CI[1.00015, 1.00015]
--------------------------------------------------------
R^2: 1.00000
========================================================
Test completed!
```

### Data format

Raw data for PyFolding should be provided in .csv files. Sample data for all of the notebooks is provided in the /examples folder.

```sh
GuHCL    CTPR2A MOPS
0        0.016262019
0.1068  -0.004045731
0.2136   0.005769455
0.3204   0.000484273
0.4272   0.001788867
0.534   -0.001676449
0.6408   0.005704187
0.7476   0.013683238
0.8544   0.018731417
0.9612   0.052727921
1.068    0.078283591
...
```

---
