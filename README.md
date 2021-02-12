# PyFolding


PyFolding is a simple Python based framework for fitting kinetic and
thermodynamic models to protein folding data.  It includes several basic models
and is extensible to enable fitting of more interesting models.

**We have now included video tutorials in the wiki to demonstrate installation
and usage of PyFolding**:
https://github.com/quantumjot/PyFolding/wiki

**Please note, we have moved the examples and Jupyter notebooks to a separate
repository**:
https://github.com/quantumjot/PyFolding-Notebooks

---
### Update for Python 3

PyFolding is not currently maintained. We have provided an updated version for Python 3, but it has not been thoroughly tested, so we cannot guarantee the results obtained. The original legacy version (for Python 2) can be found in the `py2-legacy` branch.

---
### Citation

**PyFolding: An open-source software package for graphing, analysis and simulation  
of thermodynamic and kinetic models of protein folding**  
Lowe AR, Perez-Riba A, Itzhaki L, Main E (2018) Biophys J.  
http://dx.doi.org/10.1016/j.bpj.2017.11.3779

```
@article{Lowe27112017,
  author = {Lowe, Alan R. and Perez-Riba, Albert and Itzhaki, Laura S. and Main, Ewan R.G.},
  title = {PyFolding: Open-Source Graphing, Simulation, and Analysis of the
    Biophysical Properties of Proteins},
  volume = {114},
  number = {3},
  pages = {511-521},
  year = {2018},
  doi = {10.1016/j.bpj.2017.11.3779},
  URL = {http://dx.doi.org/10.1016/j.bpj.2017.11.3779},
  eprint = {http://www.cell.com/biophysj/fulltext/S0006-3495(17)35041-5},
  journal = {Biophysical Journal}
}
```

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
### Example scripts, data and Jupyter notebooks
PyFolding can be used within IPython Notebooks (Jupyter). Several example
notebooks are provided in the notebooks folder of the PyFolding-Notebooks repo. To install Jupyter notebooks
please read the accompanying file "SETUP.md".

Example data can be found in the /examples folder.

---

### Requirements

Pyfolding requires the following additional packages:
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
================================================================================
Fitting results
================================================================================
ID: Simulated protein
Model: TwoStateEquilibrium
Optimiser: pyfolding.GlobalFit and scipy.optimize.curve_fit
Temperature: 25.00°C

(f) m        1.51913 ± 0.01683      	 95% CI[     1.51488,      1.52339]
(f) d50      5.00126 ± 0.00491      	 95% CI[     5.00002,      5.00250]
--------------------------------------------------------------------------------
R^2: 	0.99953
DOF: 	98
SS: 	9.98e-03
================================================================================


================================================================================
Fitting results
================================================================================
ID: Simulated protein
Model: TwoStateChevron
Optimiser: pyfolding.GlobalFit and scipy.optimize.curve_fit
Temperature: 25.00°C

(f) kf    100.04461 ± 0.03297      	 95% CI[   100.03628,    100.05295]
(f) mf      1.00005 ± 0.00013      	 95% CI[     1.00002,      1.00009]
(f) ku      0.00501 ± 0.00001      	 95% CI[     0.00500,      0.00501]
(f) mu      0.99988 ± 0.00013      	 95% CI[     0.99984,      0.99991]
--------------------------------------------------------------------------------
R^2: 	1.00000
DOF: 	96
SS: 	1.20e-04
================================================================================

SUCCESS - Test completed!
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
