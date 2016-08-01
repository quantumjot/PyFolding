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

### Installation

You can install Pyfolding by cloning the repo and running the setup script:
```sh
$ git clone https://github.com/quantumjot/PyFolding.git
$ python setup.py -install
```

### Testing the installation

Once installed, a simple script will execute the test function. This tests
whether the package has installed correctly.

```python
# import the pyfolding package
import pyfolding
# run a self test to check everything is installed
pyfolding.test()
```

---
### Example scripts

There are numerous example scripts in the examples folder. A simple script to
fit the equilbrium denaturation and kinetic chevron to two-state models is 
shown below:

```python
import pyfolding
import models

# load the kinetic and equilibrium folding data
chevron = pyfolding.read_kinetic_data("WT_chevron.csv")
equilibrium = folding.read_equilibrium_data("WT_equilibrium.csv")

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

### Current models supported

+ Two-state equilibrium unfolding
+ Three-state equilibrium unfolding
+ Homopoylmer Ising model for equilibrium unfolding
+ Heteropolymer Ising model for equilibrium unfolding
+ Two-state chevron
+ Two-state with moving transition state chevron
+ Three-state chevron with fast pre-equilibrium chevron 
+ Three-state chevron with fast phase chevron
+ Three-state sequential barriers chevron
+ Parallel two-state chevrons
