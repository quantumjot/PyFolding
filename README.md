# PyFolding


PyFolding is a simple Python based framework for fitting kinetic and thermodynamic models to 
protein folding data.  It includes several basic models and is extensible to enable fitting
of more interesting models.

---

### Requirements

Pyfolding has been tested on Python 2.7, and requires the following additional packages:
+ Numpy
+ Scipy
+ Matplotlib
+ (Optional) Jupyter

### Installation

You can install Pyfolding by cloning the repo and running the setup script:
```sh
$ git clone https://github.com/quantumjot/PyFolding.git
$ cd PyFolding
$ python setup.py install
```

### Testing the installation

Once installed, a simple script will execute the test function. This tests
whether the package has installed correctly.  The following script is found in ./notebooks/

```python
# import the pyfolding package
import pyfolding
# run a self test to check everything is installed
pyfolding.test()
```

Upon executing the script, the following output should be generated:

```sh
==================================================
 Fitting results 
==================================================
ID: Simulated protein
Model: TwoStateEquilibrium
Method: scipy.optimize.curve_fit 

m: 1.49978 ± 0.00000
d50: 5.00013 ± 0.00000
--------------------------------------------------
==================================================
==================================================
 Fitting results 
==================================================
ID: Simulated protein
Model: TwoStateChevron
Method: scipy.optimize.curve_fit 

kf: 100.00017 ± 0.00008
mf: 1.00002 ± 0.00000
ku: 0.00500 ± 0.00000
mu: 0.99992 ± 0.00000
--------------------------------------------------
==================================================
Test completed!
```

### Data format

Raw data for PyFolding should be provided in .csv files.

---

### Current models supported

+ Homopolymer Ising Equilibrium
+ Heteropolymer Ising Equilibrium
+ Parallel Two-State Chevron
+ Parallel Two-State Unfolding Chevron
+ Three State Chevron
+ Three State Dimeric Intermediate Equilibrium
+ Three State Fast Phase Chevron
+ Three State Monomeric Intermediate Equilibrium
+ Three State Sequential Barriers Chevron
+ Two State Chevron
+ Two State Chevron Moving Transition State
+ Two State Dimer Equilibrium
+ Two State Equilibrium
+ Two State EquilibriumSloping

A 'TemplateModel' is provided for adding additional models.

---
### Example scripts

There are numerous example scripts in the examples folder. A simple script to
fit the equilbrium denaturation and kinetic chevron to two-state models is 
shown below:

```python
import pyfolding
from pyfolding import models

# load the kinetic and equilibrium folding data
chevron = pyfolding.read_kinetic_data("WT_chevron.csv")
equilibrium = pyfolding.read_equilibrium_data("WT_equilibrium.csv")

# fit the equilibrium data to a two-state model
equilibrium.fit_func = models.TwoStateEquilibrium
equilibrium.fit()

# print the m-value, transition midpoint and stability
print equilibrium.m_value
print equilibrium.midpoint
print equilibrium.deltaG

# use the midpoint (D_50) of the equilibrium curve as the kinetic midpoint
chevron.midpoint = equilibrium.midpoint

# now fit the chevron to a two-state model
chevron.fit_func = models.TwoStateChevron
chevron.fit()

# plot the output
folding.plot_figure(equilibrium, chevron)
```

---
### IPython Notebooks

PyFolding can be used within IPython Notebooks. Several example notebooks are provided in the notebooks folder.